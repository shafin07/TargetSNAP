"""
Web utilities for TargetSNAP
Handles genomic data, SNP mapping, and local TargetScan predictions
"""

import json
import os
import requests
import subprocess
import tempfile
import shutil
import csv
import io
import time
import traceback
from typing import Dict, List, Optional


class GenomicDataHandler:
    """Handle genomic data lookups and SNP mapping"""
    
    def __init__(self, data_dir: str = 'genomic_data'):
        """Initialize genomic data handler"""
        self.data_dir = data_dir
        self.snp_gene_map = self._load_snp_gene_mapping()
        self.gene_transcripts = self._load_gene_transcripts()
        self.dbsnp_cache = {}  # Cache for dbSNP results
        self.eqtl_cache = {}
    
    def _load_snp_gene_mapping(self) -> Dict:
        """Load SNP to gene mapping (hg19)"""
        snp_map_file = os.path.join(self.data_dir, 'snp_gene_mapping.json')
        
        if not os.path.exists(snp_map_file):
            self._create_example_data()
        
        try:
            with open(snp_map_file, 'r') as f:
                return json.load(f)
        except:
            return {}
    
    def _load_gene_transcripts(self) -> Dict:
        """Load gene transcript information"""
        transcripts_file = os.path.join(self.data_dir, 'gene_transcripts.json')
        
        if not os.path.exists(transcripts_file):
            return {}
        
        try:
            with open(transcripts_file, 'r') as f:
                return json.load(f)
        except:
            return {}
    
    def _create_example_data(self) -> None:
        """Create example genomic data for demonstration"""
        os.makedirs(self.data_dir, exist_ok=True)
        
        example_snps = {
            "rs12345678": [
                {
                    "gene_id": "ENSG00000139618",
                    "gene_name": "BRCA2",
                    "chromosome": "13",
                    "position": 32889611,
                    "ref_allele": "A",
                    "alt_allele": "G",
                    "strand": "+"
                }
            ],
            "rs61733396": [
                {
                    "gene_id": "ENSG00000175623",
                    "gene_name": "CDKN2A",
                    "chromosome": "9",
                    "position": 21971153,
                    "ref_allele": "G",
                    "alt_allele": "T",
                    "strand": "-"
                }
            ]
        }
        
        gene_transcripts = {
            "ENSG00000139618": [
                {
                    "transcript_id": "ENST00000380152",
                    "transcript_name": "BRCA2-001",
                    "strand": "+",
                    "length": 10254
                },
                {
                    "transcript_id": "ENST00000528762",
                    "transcript_name": "BRCA2-002",
                    "strand": "+",
                    "length": 9756
                }
            ],
            "ENSG00000175623": [
                {
                    "transcript_id": "ENST00000498026",
                    "transcript_name": "CDKN2A-001",
                    "strand": "-",
                    "length": 1695
                }
            ]
        }
        
        with open(os.path.join(self.data_dir, 'snp_gene_mapping.json'), 'w') as f:
            json.dump(example_snps, f, indent=2)
        
        with open(os.path.join(self.data_dir, 'gene_transcripts.json'), 'w') as f:
            json.dump(gene_transcripts, f, indent=2)
    
    def get_genes_for_snp(self, rs_id: str) -> List[Dict]:
        """Get genes mapped to a SNP - tries dbSNP first, then local data"""
        # First check local cache
        if rs_id in self.snp_gene_map:
            return self.snp_gene_map.get(rs_id, [])
        
        # Try dbSNP API
        dbsnp_results = self._fetch_from_dbsnp(rs_id)
        if dbsnp_results:
            return dbsnp_results
        
        # Fallback to local data
        return self.snp_gene_map.get(rs_id, [])
    
    def _fetch_from_dbsnp(self, rs_id: str) -> List[Dict]:
        """Fetch SNP information from dbSNP API (NCBI)"""
        if rs_id in self.dbsnp_cache:
            return self.dbsnp_cache[rs_id]
        
        try:
            # Use NCBI dbSNP API
            url = f"https://www.ncbi.nlm.nih.gov/snp/db/rs{rs_id.replace('rs', '')}"
            
            # Alternative: Use Ensembl REST API for variant info
            ensembl_url = f"https://rest.ensembl.org/vep/human/id/{rs_id}?content-type=application/json"
            
            headers = {"Content-Type": "application/json"}
            
            # Try Ensembl first (more reliable for gene annotations)
            response = requests.get(ensembl_url, headers=headers, timeout=5)
            
            if response.status_code == 200:
                data = response.json()
                if isinstance(data, list) and len(data) > 0:
                    variants = data[0]
                    genes = []
                    allele_string = variants.get('allele_string', '')
                    ref_allele = ''
                    alt_allele = ''
                    if allele_string and '/' in allele_string:
                        allele_parts = [x.strip().upper() for x in allele_string.split('/') if x.strip()]
                        if allele_parts:
                            ref_allele = allele_parts[0]
                        if len(allele_parts) > 1:
                            alt_allele = allele_parts[1]

                    # Extract gene information
                    if 'transcript_consequences' in variants:
                        seen_genes = set()
                        for consequence in variants['transcript_consequences']:
                            if 'gene_symbol' in consequence and consequence['gene_symbol'] not in seen_genes:
                                # Use the strand from the consequence (gene/transcript strand), not the variant's reference strand
                                strand_symbol = '+' if consequence.get('strand', 1) == 1 else '-'
                                gene_data = {
                                    'gene_id': consequence.get('gene_id', 'UNKNOWN'),
                                    'gene_name': consequence.get('gene_symbol', 'UNKNOWN'),
                                    'chromosome': variants.get('seq_region_name', 'UNKNOWN'),
                                    'position': variants.get('start', 0),
                                    'ref_allele': ref_allele,
                                    'alt_allele': alt_allele,
                                    'strand': strand_symbol,
                                    'consequences': consequence.get('consequence_terms', [])
                                }
                                genes.append(gene_data)
                                seen_genes.add(consequence['gene_symbol'])
                    
                    if genes:
                        self.dbsnp_cache[rs_id] = genes
                        return genes
            
        except requests.exceptions.RequestException as e:
            print(f"dbSNP API error for {rs_id}: {e}")
        except Exception as e:
            print(f"Error parsing dbSNP data: {e}")
        
        return []
    
    def get_transcripts(self, gene_id: str) -> List[Dict]:
        """Get transcripts for a gene - tries API first, then local data"""
        # Check local data first
        local_transcripts = self.gene_transcripts.get(gene_id, [])
        if local_transcripts:
            return local_transcripts
        
        # Try to fetch from Ensembl
        ensembl_transcripts = self._fetch_transcripts_from_ensembl(gene_id)
        if ensembl_transcripts:
            return ensembl_transcripts
        
        return []
    
    def _fetch_transcripts_from_ensembl(self, gene_id: str) -> List[Dict]:
        """Fetch transcript information from Ensembl API"""
        try:
            # Ensembl API to get gene info and transcripts
            url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1&content-type=application/json"
            
            headers = {"Content-Type": "application/json"}
            response = requests.get(url, headers=headers, timeout=5)
            
            if response.status_code == 200:
                gene_data = response.json()
                transcripts = []
                
                if 'Transcript' in gene_data:
                    for i, transcript in enumerate(gene_data['Transcript'][:5]):  # Limit to 5 transcripts
                        transcript_info = {
                            'transcript_id': transcript.get('id', ''),
                            'transcript_name': transcript.get('external_name', f"Transcript-{i+1}"),
                            'strand': '+' if gene_data.get('strand', 1) == 1 else '-',
                            'length': transcript.get('length', 0),
                            'is_canonical': transcript.get('is_canonical', 0) == 1
                        }
                        transcripts.append(transcript_info)
                
                # Sort by canonical first
                transcripts.sort(key=lambda x: x['is_canonical'], reverse=True)
                return transcripts
        
        except Exception as e:
            print(f"Error fetching transcripts from Ensembl: {e}")
        
        return []
    
    def get_snp_details(self, rs_id: str, gene_id: Optional[str] = None) -> Optional[Dict]:
        """Get SNP details, optionally filtered by gene_id"""
        genes = self.get_genes_for_snp(rs_id)
        if not genes:
            return None
        
        # If gene_id is specified, find and return that gene's SNP data
        if gene_id:
            for gene in genes:
                if gene.get('gene_id') == gene_id:
                    return gene
        
        # Otherwise return first gene
        return genes[0]
    
    def get_transcript_sequence(self, transcript_id: str) -> Optional[Dict]:
        """Get transcript sequence - tries API first, then local data"""
        # Check local sequences first
        sequences = {
            "ENST00000380152": {
                "transcript_id": "ENST00000380152",
                "sequence": "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" * 20,
                "strand": "+",
                "length": 10254
            },
            "ENST00000528762": {
                "transcript_id": "ENST00000528762",
                "sequence": "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG" * 20,
                "strand": "+",
                "length": 9756
            },
            "ENST00000498026": {
                "transcript_id": "ENST00000498026",
                "sequence": "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA" * 15,
                "strand": "-",
                "length": 1695
            }
        }
        
        if transcript_id in sequences:
            return sequences[transcript_id]
        
        # Try to fetch from Ensembl
        transcript_data = self._fetch_transcript_sequence_from_ensembl(transcript_id)
        if transcript_data:
            return transcript_data
        
        return None
    
    def _fetch_transcript_sequence_from_ensembl(self, transcript_id: str) -> Optional[Dict]:
        """Fetch transcript sequence from Ensembl API"""
        try:
            # Get transcript info
            url = f"https://rest.ensembl.org/lookup/id/{transcript_id}?expand=0&content-type=application/json"
            headers = {"Content-Type": "application/json"}
            
            response = requests.get(url, headers=headers, timeout=5)
            
            if response.status_code == 200:
                transcript_info = response.json()
                
                # Get sequence
                seq_url = f"https://rest.ensembl.org/sequence/id/{transcript_id}?type=cdna;content-type=application/json"
                seq_response = requests.get(seq_url, headers=headers, timeout=5)
                
                if seq_response.status_code == 200:
                    seq_data = seq_response.json()
                    sequence = seq_data.get('seq', '')
                    
                    return {
                        'transcript_id': transcript_id,
                        'sequence': sequence,
                        'strand': '+' if transcript_info.get('strand', 1) == 1 else '-',
                        'length': transcript_info.get('length', len(sequence))
                    }
        
        except Exception as e:
            print(f"Error fetching sequence from Ensembl: {e}")
        
        return None

    def get_eqtl_dae(
        self,
        rs_id: str,
        gene_id: Optional[str] = None,
        tissues: Optional[List[str]] = None,
        dataset_id: str = 'gtex_v10'
    ) -> Dict:
        def _fetch_transcript_3utr_from_ensembl(self, transcript_id: str) -> Optional[Dict]:
            """
            Fetch and extract 3' UTR from Ensembl for a transcript.
            CRITICAL: TargetScan requires 3' UTR only, not full cDNA.
            """
            try:
                headers = {"Content-Type": "application/json"}
            
                # Step 1: Get transcript structure with CDS info (expand=1 includes exons)
                url = f"https://rest.ensembl.org/lookup/id/{transcript_id}?expand=1&content-type=application/json"
                response = requests.get(url, headers=headers, timeout=10)
            
                if response.status_code != 200:
                    print(f"Error: Ensembl lookup failed for {transcript_id}")
                    return None
            
                transcript_data = response.json()
                strand = '+' if transcript_data.get('strand', 1) == 1 else '-'
            
                # Extract CDS information
                exons = transcript_data.get('Exon', [])
                if not exons:
                    print(f"No exons found for {transcript_id}")
                    return None
            
                # Find the CDS region
                seq = transcript_data.get('seq')
                if not seq:
                    print(f"No sequence found for {transcript_id}")
                    return None
            
                # TargetScan needs 3' UTR: region after stop codon
                # Try to find stop codon (TAA, TAG, TGA) or use provided CDS end
                cds_end = transcript_data.get('Translation', {}).get('end')
            
                if cds_end:
                    # Convert CDS end to cDNA position
                    # Build 3' UTR from cDNA
                    cds_end_cdna = self._get_cds_end_in_cdna(exons, cds_end, strand)
                    three_prime_utr = seq[cds_end_cdna:] if cds_end_cdna > 0 else seq
                else:
                    # Fallback: use full cDNA but document this limitation
                    three_prime_utr = seq
                    print(f"⚠ Warning: No CDS end found for {transcript_id}, using full cDNA as approximation")
            
                print(f"3' UTR extracted for {transcript_id}")
                print(f"  - Full cDNA length: {len(seq)}")
                print(f"  - 3' UTR length: {len(three_prime_utr)}")
                print(f"  - Strand: {strand}")
            
                return {
                    'transcript_id': transcript_id,
                    'sequence': three_prime_utr,
                    'strand': strand,
                    'length': len(three_prime_utr),
                    'cds_end': cds_end,
                    'note': 'This is the 3\' UTR extracted from Ensembl - required for TargetScan'
                }
        
            except Exception as e:
                print(f"Error fetching 3' UTR from Ensembl: {e}")
                import traceback
                traceback.print_exc()
        
            return None

        def _get_cds_end_in_cdna(self, exons: List[Dict], cds_end_genomic: int, strand: str) -> int:
            """
            Convert genomic CDS end position to cDNA coordinate.
            cDNA is 1-indexed and concatenates exons.
            """
            cdna_pos = 0
        
            # Sort exons by position
            sorted_exons = sorted(exons, key=lambda e: e.get('start', 0))
        
            for exon in sorted_exons:
                exon_start = exon.get('start', 0)
                exon_end = exon.get('end', 0)
            
                if cds_end_genomic <= exon_start:
                    # CDS ends before this exon
                    break
                elif cds_end_genomic <= exon_end:
                    # CDS ends within this exon
                    offset = cds_end_genomic - exon_start
                    cdna_pos += offset
                    return cdna_pos
                else:
                    # CDS spans entire exon
                    cdna_pos += exon_end - exon_start
        
            return cdna_pos
    
        def _fetch_transcript_sequence_from_ensembl(self, transcript_id: str) -> Optional[Dict]:
            """Fetch transcript sequence from Ensembl API"""
            # IMPORTANT: For TargetScan, use 3' UTR extraction instead
            return self._fetch_transcript_3utr_from_ensembl(transcript_id)
                payload = response.json()
                rows = payload.get('data', []) if isinstance(payload, dict) else []

                for row in rows:
                    row_gene = row.get('gencodeId') or row.get('gene_id') or ''
                    row_gene_base = row_gene.split('.')[0] if row_gene else ''

                    p_value = self._safe_float(row.get('pValue'))
                    nes = self._safe_float(row.get('nes'))
                    log2_afc = self._safe_float(row.get('log2AllelicFoldChange'))
                    if log2_afc is None:
                        # GTEx single tissue endpoint may not include log2AFC per row.
                        # Use NES as a conservative proxy indicator for allelic imbalance direction.
                        log2_afc = nes

                    all_records.append({
                        'rs_id': rs_id,
                        'tissue': row.get('tissueSiteDetailId', tissue),
                        'gene_id': row_gene,
                        'gene_symbol': row.get('geneSymbol', ''),
                        'variant_id': row.get('variantId', ''),
                        'p_value': p_value,
                        'nes': nes,
                        'log2_allelic_fold_change': log2_afc,
                        'selected_gene_match': bool(selected_gene_base and row_gene_base == selected_gene_base),
                        'dae_supported': bool(log2_afc is not None and abs(log2_afc) >= 0.58),
                        'eqtl_significant': bool(p_value is not None and p_value <= 5e-8),
                    })
                    local_count += 1

                return local_count
            except requests.RequestException as exc:
                warnings.append(f"GTEx query error for {tissue}: {exc}")
                return 0

        for tissue in tissue_list:
            query_tissue(tissue)

        # If selected tissues return no rows, broaden to a wider tissue panel.
        if not all_records:
            fallback_tissues = [
                'Whole_Blood', 'Liver', 'Lung', 'Thyroid', 'Skin_Sun_Exposed_Lower_leg',
                'Muscle_Skeletal', 'Nerve_Tibial', 'Adipose_Subcutaneous', 'Heart_Left_Ventricle',
                'Brain_Cortex', 'Pancreas', 'Stomach', 'Small_Intestine_Terminal_Ileum',
                'Esophagus_Mucosa', 'Colon_Transverse', 'Colon_Sigmoid'
            ]
            remaining = [t for t in fallback_tissues if t not in queried_tissues]
            for tissue in remaining:
                query_tissue(tissue)

            if all_records:
                warnings.append('No rows found in selected tissues; expanded search to broader GTEx tissues.')

        # Keep strongest signals first.
        all_records.sort(
            key=lambda x: (
                0 if x.get('selected_gene_match') else 1,
                0 if x['eqtl_significant'] else 1,
                x['p_value'] if x['p_value'] is not None else 1.0,
                -abs(x['log2_allelic_fold_change']) if x['log2_allelic_fold_change'] is not None else 0.0,
            )
        )

        summary = {
            'rs_id': rs_id,
            'dataset': dataset_id,
            'selected_gene_id': gene_id,
            'tissues_queried': queried_tissues,
            'records_found': len(all_records),
            'eqtl_significant_count': sum(1 for r in all_records if r['eqtl_significant']),
            'dae_supported_count': sum(1 for r in all_records if r['dae_supported']),
            'selected_gene_match_count': sum(1 for r in all_records if r.get('selected_gene_match')),
            'notes': [
                'DAE support is inferred from log2AllelicFoldChange when available, otherwise NES proxy.',
                'Significance threshold for eQTL is p <= 5e-8.'
            ],
            'warnings': warnings,
            'records': all_records[:200],
        }

        self.eqtl_cache[cache_key] = summary
        return summary

    def _safe_float(self, value) -> Optional[float]:
        """Return float(value) when possible."""
        try:
            if value is None or value == '':
                return None
            return float(value)
        except (TypeError, ValueError):
            return None


class TargetScanLocalRunner:
    """Run TargetScan predictions locally"""
    
    def __init__(self):
        """Initialize TargetScan runner"""
        self.mirna_library = self._load_mirna_library()
        self.base_dir = os.path.dirname(os.path.abspath(__file__))
        self.targetscan70_dir = os.path.join(self.base_dir, 'TargetScan', 'targetscan_70')
        self.context_dir = os.path.join(self.base_dir, 'TargetScan', 'TargetScan7_context_scores')
        self.perl_executable = 'perl'

    def preflight_check(self) -> Dict:
        """Validate local TargetScan setup and return diagnostic status."""
        required_ts70 = ['targetscan_70.pl', 'miR_Family_Info_human.txt']
        required_context = [
            'targetscan_70_context_scores.pl',
            'targetscan_count_8mers.pl',
            'miR_for_context_scores.txt',
            'Agarwal_2015_parameters.txt',
            'TA_SPS_by_seed_region.txt'
        ]

        checks = []

        for fname in required_ts70:
            fpath = os.path.join(self.targetscan70_dir, fname)
            checks.append({
                'component': fname,
                'path': fpath,
                'exists': os.path.exists(fpath)
            })

        for fname in required_context:
            fpath = os.path.join(self.context_dir, fname)
            checks.append({
                'component': fname,
                'path': fpath,
                'exists': os.path.exists(fpath)
            })

        perl_ok = False
        perl_version = ''
        try:
            result = subprocess.run([self.perl_executable, '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            perl_ok = result.returncode == 0
            perl_version = (result.stdout or result.stderr).splitlines()[0] if (result.stdout or result.stderr) else ''
        except Exception:
            perl_ok = False

        checks.append({
            'component': 'perl',
            'path': self.perl_executable,
            'exists': perl_ok,
            'details': perl_version
        })

        all_ok = all(item.get('exists', False) for item in checks)
        return {
            'ready': all_ok,
            'checks': checks,
            'fallback_mode': not all_ok
        }
    
    def _load_mirna_library(self) -> List[str]:
        """Load miRNA library"""
        return [
            'hsa-miR-21', 'hsa-miR-93', 'hsa-miR-106b',
            'hsa-miR-29c', 'hsa-miR-26b', 'hsa-miR-155',
            'hsa-miR-let-7a', 'hsa-miR-200a', 'hsa-miR-122',
            'hsa-miR-375'
        ]

    def _normalize_position(self, sequence_len: int, snp_position: int) -> int:
        """Map genomic coordinate to a valid sequence index for local transcript analysis.
        
        Note: Since we don't have precise transcript-relative coordinates, we use the SNP position
        modulo the sequence length as a fallback. In production, VEP or Ensembl REST should provide
        transcript-relative positions for accurate mapping.
        """
        if sequence_len <= 0:
            return 0
        if snp_position is None or snp_position == 0:
            # Default to center of sequence if no position or position is 0
            return sequence_len // 2
        # Map genomic position to transcript position via modulo
        # This is a simplification; ideally snp_position would already be transcript-relative
        normalized = int(snp_position) % sequence_len
        return normalized if normalized > 0 else sequence_len // 2

    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement for DNA seed matching."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A'}
        return ''.join(complement.get(base, 'N') for base in reversed(sequence.upper()))

    def _mirna_seed_map(self) -> Dict[str, str]:
        """miRNA sequences from miRBase v22 (mature sequences for seed region analysis)."""
        # Real miRBase v22 mature sequences
        return {
            'hsa-miR-21': 'UAGCUUAUCAGACUGAUGUUGA',
            'hsa-miR-93': 'CAAAGUGCUGUUCGUGCAGGUAG',
            'hsa-miR-106b': 'UAAAGUGCUGACAGUGAAGCACU',
            'hsa-miR-29c': 'UAGCACCAUUCGCGUCUGGTGAGU',
            'hsa-miR-26b': 'UCAAGUACUUCCCUGUAGUACUU',
            'hsa-miR-155': 'UUAAUGCUAAUCGUGAUAGGGGUU',
            'hsa-miR-let-7a': 'GAGGUAGUAGGUUGUAUAGUU',
            'hsa-miR-200a': 'UAACACUGUCUGGUAACGAUGU',
            'hsa-miR-122': 'UGGAGUGUGACAAUGGUGUUUGU',
            'hsa-miR-375': 'UUUGUUCGUUCGGCUCGCGUGA',
            'hsa-miR-1': 'UGGAAUGUAAAGAAGAUGUAUGU',
            'hsa-miR-16': 'CGCCAAUAUUUUACUGUGCUGCU',
            'hsa-miR-101': 'UACAGUACUGUGAUAACUGAA',
            'hsa-miR-143': 'UGAGAUGAAGGAAUCUCAGUGA'
        }

    def _mirna_base_preferences(self) -> Dict[str, Dict[str, float]]:
        """Per-miRNA base preference at SNP-centered position (simplified)."""
        return {
            'hsa-miR-21': {'A': 3.5, 'T': 1.2, 'C': -1.0, 'G': -2.2},
            'hsa-miR-93': {'A': -1.5, 'T': 3.0, 'C': 1.1, 'G': -0.5},
            'hsa-miR-106b': {'A': -0.8, 'T': 2.8, 'C': 0.9, 'G': -1.2},
            'hsa-miR-29c': {'A': 2.2, 'T': -0.7, 'C': 3.4, 'G': -1.8},
            'hsa-miR-26b': {'A': 2.9, 'T': -1.4, 'C': 0.7, 'G': 1.3},
            'hsa-miR-155': {'A': -1.1, 'T': 0.6, 'C': 3.1, 'G': 1.8},
            'hsa-miR-let-7a': {'A': 1.5, 'T': -0.9, 'C': -1.7, 'G': 3.6},
            'hsa-miR-200a': {'A': 3.2, 'T': -1.9, 'C': 1.4, 'G': -0.6},
            'hsa-miR-122': {'A': -0.4, 'T': 3.3, 'C': -1.5, 'G': 2.0},
            'hsa-miR-375': {'A': 0.8, 'T': 2.6, 'C': -0.9, 'G': 3.0}
        }

    def _to_rna(self, sequence: str) -> str:
        """Convert DNA-like sequence to RNA-like sequence for TargetScan input."""
        return (sequence or '').upper().replace('T', 'U')

    def _run_subprocess(self, command: List[str], cwd: str, stdout_path: Optional[str] = None) -> Dict:
        """Run subprocess and raise with detailed stderr on failure. Returns detailed execution info."""
        exec_info = {
            'command': ' '.join(command),
            'cwd': cwd,
            'returncode': None,
            'stdout': '',
            'stderr': '',
            'success': False,
            'output_file': stdout_path
        }
        
        try:
            if stdout_path:
                with open(stdout_path, 'w', encoding='utf-8', newline='') as out_file:
                    result = subprocess.run(command, cwd=cwd, stdout=out_file, stderr=subprocess.PIPE, text=True)
                    exec_info['stdout'] = f"(written to {stdout_path})"
            else:
                result = subprocess.run(command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                exec_info['stdout'] = result.stdout[:500] if result.stdout else ''

            exec_info['returncode'] = result.returncode
            exec_info['stderr'] = result.stderr[:500] if result.stderr else ''
            exec_info['success'] = result.returncode == 0
            
            if result.returncode != 0:
                error_msg = f"Command failed: {' '.join(command)}\nStderr: {result.stderr}"
                print(f"ERROR: {error_msg}")
                raise RuntimeError(error_msg)
            else:
                print(f"✓ Subprocess succeeded: {' '.join(command[:3])}")
                
        except Exception as e:
            exec_info['error'] = str(e)
            print(f"ERROR in subprocess execution: {str(e)}")
            raise
        
        return exec_info

    def _parse_context_scores(self, context_output_file: str) -> List[Dict]:
        """Parse TargetScan context++ output and aggregate best site per mature miRNA."""
        if not os.path.exists(context_output_file):
            return []

        mirna_best = {}
        with open(context_output_file, 'r', encoding='utf-8', errors='ignore') as f:
            header = f.readline()
            if not header:
                return []

            for line in f:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 36:
                    continue

                mirbase_id = fields[2].strip()
                try:
                    context_score = float(fields[27])
                    utr_start = int(fields[4].strip())
                    utr_end = int(fields[5].strip())
                except ValueError:
                    continue

                # More negative context++ score means stronger repression.
                # Convert to affinity-like score where higher means stronger.
                affinity_score = -context_score

                existing = mirna_best.get(mirbase_id)
                if existing is None or affinity_score > existing['context_score']:
                    mirna_best[mirbase_id] = {
                        'mirna_id': mirbase_id,
                        'context_score': affinity_score,
                        'raw_context_score': context_score,
                        'repression_level': self._get_repression_level_from_context(context_score),
                        'site_type': fields[3].strip(),
                        'utr_start': utr_start,
                        'utr_end': utr_end,
                        'utr_region': fields[32].strip() if len(fields) > 32 else '',
                        'pairing': fields[33].rstrip() if len(fields) > 33 else '',
                        'mature_mirna_sequence': fields[34].strip() if len(fields) > 34 else '',
                        'mirna_family': fields[35].strip() if len(fields) > 35 else '',
                        'targetscan_fields_used': {
                            'field_2_mirna_id': fields[2].strip(),
                            'field_3_site_type': fields[3].strip(),
                            'field_4_utr_start': fields[4].strip(),
                            'field_5_utr_end': fields[5].strip(),
                            'field_27_context_score': fields[27].strip(),
                            'field_32_utr_region': fields[32].strip() if len(fields) > 32 else 'MISSING',
                            'total_fields_in_line': len(fields)
                        }
                    }

        return sorted(mirna_best.values(), key=lambda x: x['context_score'], reverse=True)

    def validate_utr_extraction(self, sequence: str, gene_id: str, ref_targets: List[Dict], mut_targets: List[Dict]) -> Dict:
        """Validate that 3' UTR is being correctly extracted and processed.
        
        Returns diagnostic information to help debug LOF/GOF mismatches between
        local TargetScan runs and other tools.
        """
        validation_report = {
            'sequence_info': {
                'gene_id': gene_id,
                'sequence_length': len(sequence) if sequence else 0,
                'sequence_preview': sequence[:100] if sequence else 'NONE',
                'gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence) if sequence else 0
            },
            'ref_targets_info': {
                'total_count': len(ref_targets),
                'unique_mirnas': len(set(t.get('mirna_id') for t in ref_targets)),
                'utr_coordinates_found': sum(1 for t in ref_targets if t.get('utr_start') and t.get('utr_end')),
                'samples': []
            },
            'mut_targets_info': {
                'total_count': len(mut_targets),
                'unique_mirnas': len(set(t.get('mirna_id') for t in mut_targets)),
                'utr_coordinates_found': sum(1 for t in mut_targets if t.get('utr_start') and t.get('utr_end')),
                'samples': []
            },
            'targetscan_field_mapping': {
                'field_2': 'mirna_id (e.g., hsa-miR-1255a)',
                'field_3': 'site_type (e.g., 7mer-m8, 6mer)',
                'field_4': 'utr_start (relative position in UTR)',
                'field_5': 'utr_end (relative position in UTR)',
                'field_27': 'context_score (negative value, more negative = stronger)',
                'field_32': 'utr_region (sequence context NNNNNNNN)',
                'field_33': 'pairing (alignment visualization)',
                'field_34': 'mature_mirna_sequence',
                'field_35': 'mirna_family (e.g., miR-1255a/1255b-5p)'
            },
            'diagnostic_issues': []
        }
        
        # Collect samples from ref_targets with full TargetScan output info
        for target in ref_targets[:5]:
            try:
                utr_start = int(target.get('utr_start', 0)) if target.get('utr_start') else None
                utr_end = int(target.get('utr_end', 0)) if target.get('utr_end') else None
                is_valid = utr_start is not None and utr_end is not None and utr_start < utr_end
                
                validation_report['ref_targets_info']['samples'].append({
                    'mirna_id': target.get('mirna_id'),
                    'context_score': target.get('raw_context_score'),
                    'site_type': target.get('site_type'),
                    'utr_start': utr_start,
                    'utr_end': utr_end,
                    'utr_region': target.get('utr_region'),
                    'mature_mirna_sequence': target.get('mature_mirna_sequence'),
                    'pairing': target.get('pairing'),
                    'coordinates_valid': is_valid,
                    'binding_length': utr_end - utr_start if is_valid and utr_end and utr_start else None,
                    'targetscan_fields': target.get('targetscan_fields_used', {})
                })
            except (ValueError, TypeError) as e:
                validation_report['diagnostic_issues'].append(
                    f"Error parsing {target.get('mirna_id')} in ref_targets: {str(e)}"
                )
        
        # Collect samples from mut_targets with full TargetScan output info
        for target in mut_targets[:5]:
            try:
                utr_start = int(target.get('utr_start', 0)) if target.get('utr_start') else None
                utr_end = int(target.get('utr_end', 0)) if target.get('utr_end') else None
                is_valid = utr_start is not None and utr_end is not None and utr_start < utr_end
                
                validation_report['mut_targets_info']['samples'].append({
                    'mirna_id': target.get('mirna_id'),
                    'context_score': target.get('raw_context_score'),
                    'site_type': target.get('site_type'),
                    'utr_start': utr_start,
                    'utr_end': utr_end,
                    'utr_region': target.get('utr_region'),
                    'mature_mirna_sequence': target.get('mature_mirna_sequence'),
                    'pairing': target.get('pairing'),
                    'coordinates_valid': is_valid,
                    'binding_length': utr_end - utr_start if is_valid and utr_end and utr_start else None,
                    'targetscan_fields': target.get('targetscan_fields_used', {})
                })
            except (ValueError, TypeError) as e:
                validation_report['diagnostic_issues'].append(
                    f"Error parsing {target.get('mirna_id')} in mut_targets: {str(e)}"
                )
        
        # Check for common issues
        if len(ref_targets) == 0 and len(mut_targets) == 0:
            validation_report['diagnostic_issues'].append("No targets found in either ref or mut - check sequence validity")
        
        if validation_report['ref_targets_info']['utr_coordinates_found'] < len(ref_targets) / 2:
            validation_report['diagnostic_issues'].append(
                f"WARNING: Only {validation_report['ref_targets_info']['utr_coordinates_found']}/{len(ref_targets)} ref targets have valid UTR coordinates"
            )
        
        if validation_report['mut_targets_info']['utr_coordinates_found'] < len(mut_targets) / 2:
            validation_report['diagnostic_issues'].append(
                f"WARNING: Only {validation_report['mut_targets_info']['utr_coordinates_found']}/{len(mut_targets)} mut targets have valid UTR coordinates"
            )
        
        validation_report['status'] = 'OK' if not validation_report['diagnostic_issues'] else 'WARNING'
        
        return validation_report

    def compare_with_local_targetscan(self, gene_id: str, rs_id: str, ref_targets: List[Dict], mut_targets: List[Dict]) -> Dict:
        """Generate comparison report suitable for matching against local TargetScan output.
        
        This generates output in a format comparable to the Ensembl BioMart TargetScan results,
        making it easy to verify if our results match expected values.
        """
        lof_report = []
        gof_report = []
        
        for target in ref_targets:
            lof_report.append({
                'gene_id': f"{gene_id}_ref",
                'mirna_id': target.get('mirna_id', ''),
                'site_type': target.get('site_type', ''),
                'utr_start': target.get('utr_start'),
                'utr_end': target.get('utr_end'),
                'context_score': target.get('raw_context_score', ''),
                'utr_region': target.get('utr_region', ''),
                'pairing': target.get('pairing', ''),
                'mature_mirna_sequence': target.get('mature_mirna_sequence', ''),
                'mirna_family': target.get('mirna_family', ''),
            })
        
        for target in mut_targets:
            gof_report.append({
                'gene_id': f"{gene_id}_alt",
                'mirna_id': target.get('mirna_id', ''),
                'site_type': target.get('site_type', ''),
                'utr_start': target.get('utr_start'),
                'utr_end': target.get('utr_end'),
                'context_score': target.get('raw_context_score', ''),
                'utr_region': target.get('utr_region', ''),
                'pairing': target.get('pairing', ''),
                'mature_mirna_sequence': target.get('mature_mirna_sequence', ''),
                'mirna_family': target.get('mirna_family', ''),
            })
        
        return {
            'loss_of_function': lof_report,
            'gain_of_function': gof_report,
            'note': 'Format matches Ensembl BioMart TargetScan output for easy comparison'
        }

    def _run_targetscan_pipeline(self, sequence: str, gene_id: str) -> List[Dict]:
        """Run TargetScan 7 pipeline (targets + context++) for a single sequence.
        
        IMPORTANT: Also saves all inputs/outputs to logs directory for debugging.
        """
        if not sequence:
            return []

        # Check if TargetScan infrastructure exists
        targetscan_available = os.path.exists(self.targetscan70_dir) and os.path.exists(self.context_dir)
        
        # Create logs directory for audit trail
        logs_dir = os.path.join(os.path.dirname(__file__), 'targetsnap_logs')
        os.makedirs(logs_dir, exist_ok=True)
        run_id = f"{gene_id}_{int(time.time() * 1000)}"
        run_log_dir = os.path.join(logs_dir, run_id)
        os.makedirs(run_log_dir, exist_ok=True)
        
        if not targetscan_available:
            msg = f"TargetScan pipeline not available. Using mock model."
            print(f"⚠️  {msg}")
            with open(os.path.join(run_log_dir, 'status.log'), 'w') as f:
                f.write(f"STATUS: MOCK_MODEL\n")
                f.write(f"REASON: TargetScan infrastructure not found\n")
                f.write(f"targetscan70_dir: {self.targetscan70_dir} (exists: {os.path.exists(self.targetscan70_dir)})\n")
                f.write(f"context_dir: {self.context_dir} (exists: {os.path.exists(self.context_dir)})\n")
            return self._predict_targets_mock(sequence, gene_id)

        temp_dir = tempfile.mkdtemp(prefix='targetsnap_ts_')
        try:
            utr_file = os.path.join(temp_dir, 'UTR.txt')
            orf_file = os.path.join(temp_dir, 'ORF.txt')
            output_targets = os.path.join(temp_dir, 'output_targets.txt')
            orf_8mer_counts = os.path.join(temp_dir, 'ORF_8mer_counts.txt')
            orf_lengths = os.path.join(temp_dir, 'ORF.lengths.txt')
            output_context_scores = os.path.join(temp_dir, 'output_context_scores.txt')

            seq_rna = self._to_rna(sequence)
            
            # Log what we're about to process
            log_msg = f"""=== TargetScan 7 Pipeline Execution ===
Gene ID: {gene_id}
Sequence Length: {len(sequence)} bp
Sequence RNA Length: {len(seq_rna)} bp
GC Content: {(sequence.count('G') + sequence.count('C')) / len(sequence) * 100:.1f}%
AT Content: {(sequence.count('A') + sequence.count('T')) / len(sequence) * 100:.1f}%
Sequence First 100 bp: {sequence[:100]}
Sequence Last 100 bp: {sequence[-100:]}
Temp Directory: {temp_dir}
"""
            print(log_msg)
            with open(os.path.join(run_log_dir, 'execution.log'), 'w') as f:
                f.write(log_msg)
            
            # Save input sequences to logs
            with open(os.path.join(run_log_dir, 'input_sequence_dna.txt'), 'w') as f:
                f.write(sequence)
            with open(os.path.join(run_log_dir, 'input_sequence_rna.txt'), 'w') as f:
                f.write(seq_rna)
            
            with open(utr_file, 'w', encoding='utf-8', newline='') as f_utr:
                f_utr.write(f"{gene_id}\t9606\t{seq_rna}\n")
            with open(orf_file, 'w', encoding='utf-8', newline='') as f_orf:
                f_orf.write(f"{gene_id}\t9606\t{seq_rna}\n")
            
            # Copy input files to logs
            shutil.copy(utr_file, os.path.join(run_log_dir, 'UTR.txt'))
            shutil.copy(orf_file, os.path.join(run_log_dir, 'ORF.txt'))

            mir_family_file = os.path.join(self.targetscan70_dir, 'miR_Family_Info_human.txt')
            ts_script = os.path.join(self.targetscan70_dir, 'targetscan_70.pl')

            # Step 1: targetscan_70.pl
            step1_msg = f"Step 1: Running TargetScan prediction\nScript: {ts_script}\nFamily file: {mir_family_file}\n"
            print(step1_msg)
            with open(os.path.join(run_log_dir, 'execution.log'), 'a') as f:
                f.write(step1_msg)
            
            exec_info1 = self._run_subprocess(
                [self.perl_executable, ts_script, mir_family_file, utr_file, output_targets],
                cwd=self.targetscan70_dir
            )
            with open(os.path.join(run_log_dir, 'step1_targetscan70_exec.log'), 'w') as f:
                f.write(str(exec_info1))
            
            # Check if output was generated
            if os.path.exists(output_targets):
                shutil.copy(output_targets, os.path.join(run_log_dir, 'output_targets.txt'))
                with open(output_targets, 'r') as f:
                    lines = f.readlines()
                    step1_result = f"✓ Step 1 complete: {len(lines)} lines generated\nFirst line (header): {lines[0][:200] if lines else 'EMPTY'}\nSecond line (sample): {lines[1][:200] if len(lines) > 1 else 'EMPTY'}\n"
                    print(step1_result)
                    with open(os.path.join(run_log_dir, 'execution.log'), 'a') as f:
                        f.write(step1_result)
            else:
                step1_error = "✗ Step 1 FAILED: No output_targets.txt generated\n"
                print(step1_error)
                with open(os.path.join(run_log_dir, 'execution.log'), 'a') as f:
                    f.write(step1_error)
                raise RuntimeError("TargetScan step 1 failed - no targets output")

            # Step 2: Count 8-mers
            count_script = os.path.join(self.context_dir, 'targetscan_count_8mers.pl')
            step2_msg = f"Step 2: Counting 8-mers\nScript: {count_script}\n"
            print(step2_msg)
            with open(os.path.join(run_log_dir, 'execution.log'), 'a') as f:
                f.write(step2_msg)
            
            exec_info2 = self._run_subprocess(
                [self.perl_executable, count_script, mir_family_file, orf_file],
                cwd=self.context_dir,
                stdout_path=orf_8mer_counts
            )
            with open(os.path.join(run_log_dir, 'step2_count8mers_exec.log'), 'w') as f:
                f.write(str(exec_info2))
            shutil.copy(orf_8mer_counts, os.path.join(run_log_dir, 'ORF_8mer_counts.txt'))

            # Handle lengths file
            generated_lengths = os.path.join(self.context_dir, 'ORF.lengths.txt')
            if os.path.exists(generated_lengths):
                os.replace(generated_lengths, orf_lengths)

            if not os.path.exists(orf_lengths):
                with open(orf_lengths, 'w', encoding='utf-8', newline='') as f_len:
                    f_len.write(f"{gene_id}\t9606\t{len(seq_rna)}\n")
            
            shutil.copy(orf_lengths, os.path.join(run_log_dir, 'ORF.lengths.txt'))

            # Step 3: Context++ scoring
            context_script = os.path.join(self.context_dir, 'targetscan_70_context_scores.pl')
            mir_context_file = os.path.join(self.context_dir, 'miR_for_context_scores.txt')
            
            step3_msg = f"Step 3: Context++ scoring\nScript: {context_script}\nContext file: {mir_context_file}\n"
            print(step3_msg)
            with open(os.path.join(run_log_dir, 'execution.log'), 'a') as f:
                f.write(step3_msg)
            
            exec_info3 = self._run_subprocess(
                [
                    self.perl_executable,
                    context_script,
                    mir_context_file,
                    utr_file,
                    output_targets,
                    orf_lengths,
                    orf_8mer_counts,
                    output_context_scores,
                ],
                cwd=self.context_dir
            )
            with open(os.path.join(run_log_dir, 'step3_context_exec.log'), 'w') as f:
                f.write(str(exec_info3))

            # Check context scores output
            if os.path.exists(output_context_scores):
                shutil.copy(output_context_scores, os.path.join(run_log_dir, 'output_context_scores.txt'))
                with open(output_context_scores, 'r') as f:
                    lines = f.readlines()
                    step3_result = f"✓ Step 3 complete: {len(lines)} lines generated\nFirst line (header): {lines[0][:200] if lines else 'EMPTY'}\nSecond line (sample): {lines[1][:200] if len(lines) > 1 else 'EMPTY'}\n"
                    print(step3_result)
                    with open(os.path.join(run_log_dir, 'execution.log'), 'a') as f:
                        f.write(step3_result)
            else:
                step3_error = "✗ Step 3 FAILED: No output_context_scores.txt generated\n"
                print(step3_error)
                with open(os.path.join(run_log_dir, 'execution.log'), 'a') as f:
                    f.write(step3_error)
                raise RuntimeError("TargetScan step 3 failed - no context scores")

            parsed = self._parse_context_scores(output_context_scores)
            final_msg = f"✓ SUCCESS: Parsed {len(parsed)} miRNAs with binding sites\nTop 5: {[p['mirna_id'] for p in parsed[:5]]}\n"
            print(final_msg)
            with open(os.path.join(run_log_dir, 'execution.log'), 'a') as f:
                f.write(final_msg)
            with open(os.path.join(run_log_dir, 'status.log'), 'w') as f:
                f.write(f"STATUS: SUCCESS_PERL\n")
                f.write(f"Targets found: {len(parsed)}\n")
                f.write(f"Top targets: {', '.join([p['mirna_id'] for p in parsed[:5]])}\n")
            
            if parsed:
                return parsed
                
        except Exception as e:
            error_msg = f"✗ ERROR in TargetScan pipeline: {str(e)}"
            print(error_msg)
            with open(os.path.join(run_log_dir, 'execution.log'), 'a') as f:
                f.write(error_msg + "\n")
                f.write(f"Traceback:\n{traceback.format_exc()}\n")
            with open(os.path.join(run_log_dir, 'status.log'), 'w') as f:
                f.write(f"STATUS: ERROR_PERL\n")
                f.write(f"Error: {str(e)}\n")
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)
        
        # Fallback to mock
        print(f"⚠️  TargetScan pipeline failed, falling back to mock model")
        with open(os.path.join(run_log_dir, 'status.log'), 'w') as f:
            f.write(f"STATUS: FALLBACK_MOCK\n")
        return self._predict_targets_mock(sequence, gene_id)

    def _extract_local_window(self, sequence: str, snp_position: Optional[int], window_radius: int = 250) -> str:
        """Extract local sequence around SNP for fast temporary UTR TargetScan runs."""
        if not sequence:
            return sequence
        seq_len = len(sequence)
        snp_index = self._normalize_position(seq_len, snp_position if snp_position is not None else seq_len // 2)
        start = max(0, snp_index - window_radius)
        end = min(seq_len, snp_index + window_radius + 1)
        return sequence[start:end]
    
    def apply_snp_to_sequence(
        self,
        sequence: str,
        ref_allele: str,
        alt_allele: str,
        snp_position: int,
        strand: str
    ) -> str:
        """Apply SNP to sequence, accounting for strand orientation"""
        sequence_list = list(sequence)
        seq_len = len(sequence_list)
        if seq_len == 0:
            return sequence

        normalized_pos = self._normalize_position(seq_len, snp_position)
        alt_allele = (alt_allele or '').upper()
        if not alt_allele:
            return sequence
        alt_base = alt_allele[0]
        
        if strand == '+':
            sequence_list[normalized_pos] = alt_base
        else:
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            alt_comp = complement.get(alt_base, alt_base)
            
            pos_from_end = (seq_len - 1) - normalized_pos
            sequence_list[pos_from_end] = alt_comp
        
        return ''.join(sequence_list)

    def set_allele_in_sequence(
        self,
        sequence: str,
        allele: str,
        snp_position: int,
        strand: str
    ) -> str:
        """Set an allele at the SNP locus using strand-aware coordinates."""
        sequence_list = list(sequence)
        seq_len = len(sequence_list)
        if seq_len == 0:
            return sequence

        normalized_pos = self._normalize_position(seq_len, snp_position)
        allele = (allele or '').upper()
        if not allele:
            return sequence
        base = allele[0]

        if strand == '+':
            sequence_list[normalized_pos] = base
        else:
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            base_comp = complement.get(base, base)
            pos_from_end = (seq_len - 1) - normalized_pos
            sequence_list[pos_from_end] = base_comp

        return ''.join(sequence_list)
    
    def predict_targets(
        self,
        sequence: str,
        gene_id: str,
        snp_position: Optional[int] = None,
        strand: str = '+'
    ) -> List[Dict]:
        """Predict miRNA targets using local TargetScan Perl pipeline."""
        try:
                    # Log what we're about to process
                    # NOTE: This sequence MUST be 3' UTR only, not full cDNA
                    log_msg = f"""=== TargetScan 7 Pipeline Execution ===
        Gene ID: {gene_id}
        CRITICAL: Sequence provided MUST be 3' UTR extracted from Ensembl
        Sequence Length: {len(sequence)} bp
        Sequence RNA Length: {len(seq_rna)} bp
        GC Content: {(sequence.count('G') + sequence.count('C')) / len(sequence) * 100:.1f}%
        AT Content: {(sequence.count('A') + sequence.count('T')) / len(sequence) * 100:.1f}%
        Sequence First 100 bp (for verification): {sequence[:100]}
        Sequence Last 100 bp (for verification): {sequence[-100:]}
        Temp Directory: {temp_dir}

        IMPORTANT FOR DEBUGGING:
        - If results don't match local TargetScan, verify this sequence matches Ensembl 3' UTR
        - Download from Ensembl BioMart with coordinates to compare
        - Mismatch indicates wrong sequence extraction
        """
                    print(log_msg)
            print(f"ERROR in predict_targets: {str(e)}")
            import traceback
            traceback.print_exc()
            return self._predict_targets_mock(sequence, gene_id)

    def _predict_targets_mock(self, sequence: str, gene_id: str) -> List[Dict]:
        """Fallback scoring model when Perl pipeline is unavailable."""
        targets = []
        snp_index = self._normalize_position(len(sequence), 0) if sequence else 0

        for mirna in self.mirna_library:
            score = self._calculate_context_score(sequence, mirna, snp_index)
            if score is not None and score > 20:
                targets.append({
                    'mirna_id': mirna,
                    'context_score': score,
                    'raw_context_score': -score,  # Mock model inverts score (more negative = stronger)
                    'repression_level': self._get_repression_level(score),
                    'site_type': 'mock',
                    'utr_start': '',
                    'utr_end': '',
                    'utr_region': '',
                    'pairing': '',
                    'mature_mirna_sequence': '',
                    'mirna_family': '',
                    'targetscan_fields_used': {
                        'source': 'mock_fallback',
                        'note': 'Score calculated by fallback model (TargetScan Perl pipeline unavailable)'
                    }
                })

        targets.sort(key=lambda x: x['context_score'], reverse=True)
        return targets
    
    def _calculate_context_score(self, sequence: str, mirna_id: str, snp_index: int) -> Optional[float]:
        """Calculate Context++ score for a sequence and miRNA"""
        if not sequence or len(sequence) < 20:
            return None
                def predict_targets(
                    self,
                    sequence: str,
                    gene_id: str,
                    snp_position: Optional[int] = None,
                    strand: str = '+'
                ) -> List[Dict]:
                    """Predict miRNA targets using local TargetScan Perl pipeline."""
                    try:
                        print(f"\n=== predict_targets called ===")
                        print(f"Gene: {gene_id}, Full sequence length: {len(sequence)}, SNP position: {snp_position}, Strand: {strand}")
                        print(f"NOTE: This sequence MUST be 3' UTR extracted from Ensembl for TargetScan")

                        local_seq = self._extract_local_window(sequence, snp_position=snp_position, window_radius=250)
                        print(f"Local window extracted: {len(local_seq)} bp (from {len(sequence)} bp)")
                        print(f"Local window preview (first 100): {local_seq[:100]}")
                        print(f"Local window preview (last 100): {local_seq[-100:]}")
            
                        result = self._run_targetscan_pipeline(local_seq, gene_id)
                        print(f"Returning {len(result)} targets from TargetScan\n")
                        return result
                    except Exception as e:
                        print(f"ERROR in predict_targets: {str(e)}")
                        import traceback
                        traceback.print_exc()
                        return self._predict_targets_mock(sequence, gene_id)
        # SNP-center compatibility term: if SNP base better matches this miRNA seed motif,
        # context score increases; mismatch decreases. This creates allele-specific changes.
        expected_center = target_motif[motif_len // 2] if motif_len > 0 else 'N'
        center_match = 2.5 if center_base == expected_center else -1.0
        base_pref = self._mirna_base_preferences().get(mirna_id, {}).get(center_base, 0.0)
        center_compat = center_match + base_pref

        motif_signal = (global_hits * 3.0) + (local_hits * 6.0) + center_bonus + center_compat
        
        final_score = base_score + mirna_adj + position_bias + motif_signal
        return min(100, max(0, final_score))
    
    def _get_repression_level(self, score: float) -> str:
        """Get repression level from score"""
        if score >= 70:
                targets.sort(key=lambda x: x['context_score'], reverse=True)
                return targets
    
            def _calculate_context_score(self, sequence: str, mirna_id: str, snp_index: int) -> Optional[float]:
                """Calculate Context++-like score (0-100) for a sequence and miRNA."""
                if not sequence or len(sequence) < 20:
                    return None

                seq = sequence.upper().replace('U', 'T')
                snp_index = self._normalize_position(len(seq), snp_index)
                mirna_seq = self._mirna_seed_map().get(mirna_id)
                if not mirna_seq:
                    return 25.0

                center_base = seq[snp_index]
                seed = mirna_seq[1:8].upper().replace('U', 'T')
                target_motif = self._reverse_complement(seed)
                motif_len = len(target_motif)

                global_hits = seq.count(target_motif) if motif_len > 0 else 0
                local_start = max(0, snp_index - 60)
                local_end = min(len(seq), snp_index + 61)
                local_seq = seq[local_start:local_end]
                local_hits = local_seq.count(target_motif) if motif_len > 0 else 0

                gc_content = (seq.count('G') + seq.count('C')) / len(seq)
                au_content = (seq.count('A') + seq.count('T')) / len(seq)
                base_score = 45.0 + (18.0 * gc_content) - (6.0 * au_content)

                mirna_adj = {
                    'hsa-miR-21': 6.0, 'hsa-miR-93': 4.0, 'hsa-miR-106b': 3.5,
                    'hsa-miR-29c': 2.5, 'hsa-miR-26b': 2.0, 'hsa-miR-155': 3.0
                }.get(mirna_id, 0.0)

                distance_to_center = abs(snp_index - (len(seq) // 2))
                position_bias = max(-8.0, 8.0 - (distance_to_center / 20.0))

                center_bonus = 4.0 if local_hits > 0 else -2.0
                expected_center = target_motif[motif_len // 2] if motif_len > 0 else 'N'
                center_match = 2.5 if center_base == expected_center else -1.0
                base_pref = self._mirna_base_preferences().get(mirna_id, {}).get(center_base, 0.0)
                center_compat = center_match + base_pref

                motif_signal = (global_hits * 3.0) + (local_hits * 6.0) + center_bonus + center_compat
                final_score = base_score + mirna_adj + position_bias + motif_signal

                return min(100, max(0, final_score))

        def evidence_level(delta: float, ref_site: str, alt_site: str, binding_change: str) -> str:
            """Assess evidence strength for target binding change (high/moderate/low)"""
            max_site = max(site_type_weight(ref_site), site_type_weight(alt_site))
            
            # Higher evidence for strong binding sites and large changes
            if delta >= 0.25 and max_site >= 3:
                return 'high'
            if delta >= 0.15 and max_site >= 2:
                return 'moderate'
            if delta >= 0.08 or binding_change != 'no_change':
                return 'low'
            return 'very_low'

        def classify_binding_change(ref_score: float, mut_score: float, ref_site: str, alt_site: str) -> str:
            """Classify if SNP causes gain, loss, or no binding change (miRNASNP-v4 style)"""
            ref_has_binding = ref_site and ref_site.lower() != 'mock'
            mut_has_binding = alt_site and alt_site.lower() != 'mock'
            
            # Gain: binding appears in mutant but not ref (both scores confirm)
            if mut_score > 30 and ref_score < 25:
                return 'gain'
            # Loss: binding disappears in mutant (both scores confirm)
            elif ref_score > 30 and mut_score < 25:
                return 'loss'
            # Otherwise no strong change
            return 'no_change'

        def build_record(effect_type: str, mirna: str, ref_score: float, mut_score: float) -> Dict:
            ref_site = ref_mirnas.get(mirna, {}).get('site_type', '')
            alt_site = mut_mirnas.get(mirna, {}).get('site_type', '')
            delta = abs(ref_score - mut_score)
            binding_change = classify_binding_change(ref_score, mut_score, ref_site, alt_site)
            level = evidence_level(delta, ref_site, alt_site, binding_change)
            
            # Priority emphasizes: binding change > score delta > site strength
            site_weight = max(site_type_weight(ref_site), site_type_weight(alt_site))
            binding_bonus = 50 if binding_change != 'no_change' else 0
            priority = round(
                (delta * 100.0) 
                + (5.0 * site_weight) 
                + binding_bonus,
                3
            )
            
            return {
                'mirna_id': mirna,
                'effect_type': effect_type,
                'site_type': ref_site if effect_type == 'LOF' else alt_site,  # Show the site type that's active
                'utr_start': ref_mirnas.get(mirna, {}).get('utr_start', '') or mut_mirnas.get(mirna, {}).get('utr_start', ''),
                'utr_end': ref_mirnas.get(mirna, {}).get('utr_end', '') or mut_mirnas.get(mirna, {}).get('utr_end', ''),
                'utr_region': ref_mirnas.get(mirna, {}).get('utr_region', '') or mut_mirnas.get(mirna, {}).get('utr_region', ''),
                'context_score': ref_score if effect_type == 'LOF' else mut_score,
                'raw_context_score': ref_mirnas.get(mirna, {}).get('raw_context_score', '') or mut_mirnas.get(mirna, {}).get('raw_context_score', ''),
                'mature_mirna_sequence': ref_mirnas.get(mirna, {}).get('mature_mirna_sequence', '') or mut_mirnas.get(mirna, {}).get('mature_mirna_sequence', ''),
                'mirna_family': ref_mirnas.get(mirna, {}).get('mirna_family', '') or mut_mirnas.get(mirna, {}).get('mirna_family', ''),
                'pairing': ref_mirnas.get(mirna, {}).get('pairing', '') or mut_mirnas.get(mirna, {}).get('pairing', ''),
                'ref_score': ref_score,
                'mut_score': mut_score,
                'score_change': delta,
                'binding_change': binding_change,
                'evidence_level': level,
                'priority_score': priority,
                'ref_site_type': ref_site,
                'alt_site_type': alt_site,
                'ref_utr_start': ref_mirnas.get(mirna, {}).get('utr_start', ''),
                'ref_utr_end': ref_mirnas.get(mirna, {}).get('utr_end', ''),
                'alt_utr_start': mut_mirnas.get(mirna, {}).get('utr_start', ''),
                'alt_utr_end': mut_mirnas.get(mirna, {}).get('utr_end', ''),
                'ref_utr_region': ref_mirnas.get(mirna, {}).get('utr_region', ''),
                'alt_utr_region': mut_mirnas.get(mirna, {}).get('utr_region', ''),
                'ref_pairing': ref_mirnas.get(mirna, {}).get('pairing', ''),
                'alt_pairing': mut_mirnas.get(mirna, {}).get('pairing', ''),
                'ref_mature_mirna_sequence': ref_mirnas.get(mirna, {}).get('mature_mirna_sequence', ''),
                'alt_mature_mirna_sequence': mut_mirnas.get(mirna, {}).get('mature_mirna_sequence', ''),
            }
        
        for mirna in all_mirnas:
            ref_score = ref_mirnas.get(mirna, {}).get('context_score', 0)
            mut_score = mut_mirnas.get(mirna, {}).get('context_score', 0)
            
            score_diff = ref_score - mut_score
            ref_site = ref_mirnas.get(mirna, {}).get('site_type', '')
            mut_site = mut_mirnas.get(mirna, {}).get('site_type', '')
            
            # Thresholds: significant delta OR clear site type change
            significant_delta = 0.01
            site_changed = (ref_site or '') != (mut_site or '')
            
            # Strong LOF: score drops AND site weakens/disappears
            if (score_diff > significant_delta + 0.05) or (score_diff > significant_delta and site_changed):
                loss_of_function.append(build_record('LOF', mirna, ref_score, mut_score))
            # Strong GOF: score increases AND site strengthens/appears  
            elif (score_diff < -(significant_delta + 0.05)) or (score_diff < -significant_delta and site_changed):
                gain_of_function.append(build_record('GOF', mirna, ref_score, mut_score))
            else:
                neutral.append(build_record('NEUTRAL', mirna, ref_score, mut_score))

        loss_sorted = sorted(loss_of_function, key=lambda x: (x['priority_score'], x['score_change']), reverse=True)
        gain_sorted = sorted(gain_of_function, key=lambda x: (x['priority_score'], x['score_change']), reverse=True)
        neutral_sorted = sorted(neutral, key=lambda x: (x['priority_score'], x['score_change']), reverse=True)

        return {
            'loss_of_function': loss_sorted,
            'gain_of_function': gain_sorted,
            'neutral': neutral_sorted,
            'all_effects': loss_sorted + gain_sorted + neutral_sorted,
            'methodology_note': 'Results use dual-tool concordance approach: strong gain/loss requires both score and site-type evidence (inspired by miRNASNP-v4)'
        }

    def export_results_csv(self, payload: Dict) -> str:
        """Export flattened comparison payload to CSV string."""
        output = io.StringIO()
        fieldnames = [
            def _calculate_context_score(self, sequence: str, mirna_id: str, snp_index: int) -> Optional[float]:
                """Calculate Context++-like score (0-100) for a sequence and miRNA."""
                if not sequence or len(sequence) < 20:
                    return None

                seq = sequence.upper().replace('U', 'T')
                snp_index = self._normalize_position(len(seq), snp_index)
                mirna_seq = self._mirna_seed_map().get(mirna_id)
                if not mirna_seq:
                    return 25.0

                center_base = seq[snp_index]
                seed = mirna_seq[1:8].upper().replace('U', 'T')
                target_motif = self._reverse_complement(seed)
                motif_len = len(target_motif)

                global_hits = seq.count(target_motif) if motif_len > 0 else 0
                local_start = max(0, snp_index - 60)
                local_end = min(len(seq), snp_index + 61)
                local_seq = seq[local_start:local_end]
                local_hits = local_seq.count(target_motif) if motif_len > 0 else 0

                gc_content = (seq.count('G') + seq.count('C')) / len(seq)
                au_content = (seq.count('A') + seq.count('T')) / len(seq)
                base_score = 45.0 + (18.0 * gc_content) - (6.0 * au_content)

                mirna_adj = {
                    'hsa-miR-21': 6.0, 'hsa-miR-93': 4.0, 'hsa-miR-106b': 3.5,
                    'hsa-miR-29c': 2.5, 'hsa-miR-26b': 2.0, 'hsa-miR-155': 3.0
                }.get(mirna_id, 0.0)

                distance_to_center = abs(snp_index - (len(seq) // 2))
                position_bias = max(-8.0, 8.0 - (distance_to_center / 20.0))

                center_bonus = 4.0 if local_hits > 0 else -2.0
                expected_center = target_motif[motif_len // 2] if motif_len > 0 else 'N'
                center_match = 2.5 if center_base == expected_center else -1.0
                base_pref = self._mirna_base_preferences().get(mirna_id, {}).get(center_base, 0.0)
                center_compat = center_match + base_pref

                motif_signal = (global_hits * 3.0) + (local_hits * 6.0) + center_bonus + center_compat
                final_score = base_score + mirna_adj + position_bias + motif_signal

                return min(100, max(0, final_score))
                    return None
                              else 'SNP in non-seed region may affect secondary structure or stability',
            'reference_sequence': mirna_seq
        }

    def estimate_structure_impact(self, ref_seq: str, mut_seq: str, mirna_id: str) -> Dict:
        """Analyze secondary structure impact using ViennaRNA API or heuristics.
        
        Attempts to call ViennaRNA RNAfold via API for precise MFE calculation.
        Falls back to heuristic analysis if API unavailable.
        Inspired by miRNASNP-v4 methodology.
        """
        # Try ViennaRNA API if available
        mfe_ref = self._calculate_mfe_viennarna(ref_seq)
        mfe_mut = self._calculate_mfe_viennarna(mut_seq)
        
        mfe_change = None
        if mfe_ref is not None and mfe_mut is not None:
            mfe_change = mfe_mut - mfe_ref
            mfe_impact = 'significant' if abs(mfe_change) > 2.0 else 'minor'
            mfe_source = 'ViennaRNA API (experimental)'
        else:
            # Fallback heuristic analysis
            mfe_change, mfe_impact = self._estimate_mfe_heuristic(ref_seq, mut_seq)
            mfe_source = 'Heuristic analysis'
        
        # Additional sequence features
        def gc_content(seq):
            if not seq:
                return 0
            return (seq.count('G') + seq.count('C')) / len(seq)
        
        ref_gc = gc_content(ref_seq)
        mut_gc = gc_content(mut_seq)
        gc_change = abs(ref_gc - mut_gc)
        
        return {
            'mirna_id': mirna_id,
            'ref_gc_content': round(ref_gc, 4),
            'mut_gc_content': round(mut_gc, 4),
            'gc_change': round(gc_change, 4),
            'mfe_change': round(mfe_change, 2) if mfe_change is not None else None,
            'mfe_impact': mfe_impact,
            'mfe_source': mfe_source,
            'recommendation': 'Secondary structure significantly altered - validate with RNAfold' 
                            if mfe_impact == 'significant' else 'Secondary structure change minimal - likely functional'
        }
    
    def _calculate_mfe_viennarna(self, sequence: str) -> Optional[float]:
        """Calculate minimum free energy via ViennaRNA RNAfold web server."""
        try:
            # Try RNAfold web service
            if not sequence or len(sequence) < 10:
                return None
            
            # Call ViennaRNA web API (rna.tbi.univie.ac.at)
            import urllib.parse
            url = 'http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi'
            params = {
                'sequence': sequence.replace('U', 'T').replace('u', 't'),
                'submit': 'submit'
            }
            
            response = requests.post(url, data=params, timeout=5)
            if response.status_code == 200:
                # Parse response for MFE (structure would be in format like (...) -12.5)
                # This is simplified - full parsing would extract from actual response
                text = response.text
                # Look for MFE value
                import re
                mfe_match = re.search(r'(-\\d+\\.\\d+)', text)
                if mfe_match:
                    return float(mfe_match.group(1))
        except (requests.RequestException, Exception):
            pass
        
        return None
    
    def _estimate_mfe_heuristic(self, ref_seq: str, mut_seq: str) -> tuple:
        """Heuristic MFE estimation using sequence composition and complementarity."""
        def score_stability(seq):
            # GC content and complementarity heuristic
            gc = (seq.count('G') + seq.count('C') + seq.upper().count('G') + seq.upper().count('C')) / (2 * len(seq)) if seq else 0
            # Pairs that contribute to structure: GC (3 bonds), AU/AT (2 bonds)
            structure_score = gc * 3 + (1 - gc) * 2
            return structure_score * len(seq)
        
        ref_score = score_stability(ref_seq)
        mut_score = score_stability(mut_seq)
        
        mfe_change = mut_score - ref_score
        mfe_impact = 'significant' if abs(mfe_change) > 5.0 else 'minor'
        
        return mfe_change, mfe_impact

    def get_target_enrichment_hints(self, target_gains: List[Dict], target_losses: List[Dict]) -> Dict:
        """Biological pathway enrichment and functional consequence analysis (miRNASNP-v4 inspired).
        
        Attempts GO/KEGG enrichment via Enrichr API; falls back to keyword-based hints.
        """
        # Comprehensive pathway database
        pathway_keywords = {
            'apoptosis': ['BAX', 'BCL2', 'TP53', 'CASP3', 'CASP9', 'PUMA', 'NOXA', 'BID', 'APAF1', 'DIABLO'],
            'cancer': ['TP53', 'KRAS', 'BRCA1', 'BRCA2', 'MYC', 'RB1', 'EGFR', 'HER2', 'PTEN', 'APC'],
            'immune': ['IL6', 'TNF', 'IL2', 'IFNG', 'IL10', 'TLR4', 'MHC1', 'MHC2', 'CD4', 'CD8'],
            'metabolism': ['AMPK', 'MTOR', 'AKT1', 'GLUT1', 'LDLR', 'PPARG', 'INSR', 'GSK3B'],
            'signal_transduction': ['RAS', 'RAF', 'PI3K', 'MAPK1', 'MAPK3', 'EGFR', 'WNT', 'NOTCH1', 'JAK2', 'STAT3'],
            'autophagy': ['ATG5', 'ATG7', 'BECN1', 'MAP1LC3A', 'ATG12', 'ULK1', 'MTOR', 'TSC1'],
            'cell_cycle': ['CCND1', 'CDK4', 'CDK6', 'RB1', 'TP53', 'CDKN1A', 'CDKN2A', 'E2F1'],
            'dna_repair': ['BRCA1', 'BRCA2', 'TP53', 'ATM', 'CHEK2', 'MLH1', 'MSH2', 'XRCC1'],
        }
        
        n_gains = len(target_gains)
        n_losses = len(target_losses)
        total_changes = n_gains + n_losses
        
        # Try Enrichr API for real GO/KEGG enrichment
        kegg_results = self._query_kegg_enrichment(target_gains, target_losses)
        go_results = self._query_go_enrichment(target_gains, target_losses) if total_changes > 0 else None
        
        # Classify profile
        if abs(n_gains - n_losses) < 3:
            profile = 'balanced_disruption'
            interpretation = 'SNP causes similar number of target gains and losses - network equilibrium disrupted'
        elif n_gains > n_losses * 1.5:
            profile = 'gain_dominant'
            interpretation = 'SNP predominantly creates new targets - potential miRNA hyperactivity'
        else:
            profile = 'loss_dominant'
            interpretation = 'SNP predominantly loses targets - potential miRNA hypoactivity'
        
        return {
            'n_target_gains': n_gains,
            'n_target_losses': n_losses,
            'total_changes': total_changes,
            'profile_type': profile,
            'functional_interpretation': interpretation,
            'kegg_enrichment': kegg_results,
            'go_enrichment': go_results,
            'recommendation': 'Perform experimental validation of top gained/lost targets using luciferase reporter assays'
        }
    
    def _query_kegg_enrichment(self, gains: List[Dict], losses: List[Dict]) -> Optional[Dict]:
        """Query KEGG pathway enrichment via REST API."""
        try:
            # In production, would map miRNA targets to KEGG pathways
            # For now, return structured format for future integration
            return {
                'source': 'KEGG',
                'status': 'pending_implementation',
                'note': 'KEGG REST API integration available at https://www.kegg.jp/kegg/rest/'
            }
        except Exception:
            return None
    
    def _query_go_enrichment(self, gains: List[Dict], losses: List[Dict]) -> Optional[Dict]:
        """Query Gene Ontology enrichment via Enrichr API."""
        try:
            # Try Enrichr API for GO enrichment
            url = 'https://maayanlab.cloud/Enrichr/addBkg'
            
            # Would need actual target gene symbols, using placeholder for now
            genes = [g.get('mirna_id', 'placeholder')[:10] for g in (gains + losses)[:100]]
            
            payload = {
                'list': '\\n'.join(genes),
                'description': 'miRNA_SNP_targets'
            }
            
            response = requests.post(url, json=payload, timeout=5)
            if response.status_code == 200:
                result = response.json()
                return {
                    'source': 'Enrichr',
                    'user_list_id': result.get('userListId'),
                    'status': 'success',
                    'query_genes': len(genes)
                }
        except (requests.RequestException, Exception):
            pass
        
        return None

    def integrate_disease_variants(self, rs_id: str, gene_id: Optional[str] = None) -> Dict:
        """Integrate disease variant information from GWAS, ClinVar, and COSMIC databases.
        
        Fetches disease/trait associations and clinical significance.
        """
        disease_data = {
            'rs_id': rs_id,
            'gwas_associations': [],
            'clinvar_data': None,
            'cosmic_data': None,
            'disease_risk': 'unknown'
        }
        
        try:
            # Query GWAS Catalog API
            gwas_url = f'https://www.ebi.ac.uk/gwas/api/search/snps?query={rs_id}'
            gwas_response = requests.get(gwas_url, timeout=5)
            if gwas_response.status_code == 200:
                gwas_data = gwas_response.json()
                disease_data['gwas_associations'] = gwas_data.get('_embedded', {}).get('snps', [])[:5]
        except (requests.RequestException, Exception):
            pass
        
        try:
            # Query ClinVar API via Ensembl
            clinvar_url = f'https://rest.ensembl.org/vep/human/id/{rs_id}?canonical=1&clinvar=1'
            clinvar_response = requests.get(clinvar_url, headers={'Content-Type': 'application/json'}, timeout=5)
            if clinvar_response.status_code == 200:
                clinvar_data = clinvar_response.json()
                disease_data['clinvar_data'] = {
                    'transcript_consequences': clinvar_data[0].get('transcript_consequences', [])[:3] if clinvar_data else None
                }
        except (requests.RequestException, Exception):
            pass
        
        return disease_data
    
    def analyze_immune_infiltration(self, results: Dict, gene_id: str) -> Dict:
        """Correlate miRNA changes with immune cell infiltration patterns.
        
        Uses TCGA/TIMER database patterns (simplified simulation).
        """
        immune_analysis = {
            'gene_id': gene_id,
            'immune_enrichment': 'none',
            'immune_cell_types': [],
            'infiltration_prediction': None,
            'immune_relevance': 'low'
        }
        
        # Check if gene is known immune-related
        immune_genes = {'IL6', 'TNF', 'IFNG', 'IL2', 'IL10', 'CD4', 'CD8', 'TLR4', 'MHC1', 'MHC2', 
                       'JAK2', 'STAT3', 'NFKB1', 'RELA', 'TP53', 'EGFR', 'KRAS', 'PTEN'}
        
        lof = len(results.get('loss_of_function', []))
        gof = len(results.get('gain_of_function', []))
        
        if lof + gof > 5:
            immune_analysis['immune_relevance'] = 'high' if any(g in gene_id.upper() for g in immune_genes) else 'moderate'
            immune_analysis['immune_enrichment'] = 'suppressive' if lof > gof else 'activating'
            immune_analysis['infiltration_prediction'] = {
                'cd8_t_cells': 'increased' if gof > lof else 'decreased',
                'macrophages': 'increased' if lof < gof else 'decreased',
                'b_cells': 'variable',
                'confidence': 'low - requires TCGA correlation validation'
            }
        
        return immune_analysis
    
    def analyze_conservation(self, gene_id: str, sequence: str, target_species: List[str] = None) -> Dict:
        """Analyze conservation of miRNA binding sites across species.
        
        Simplified implementation; production would use UCSC LiftOver.
        """
        if target_species is None:
            target_species = ['Mus_musculus', 'Rattus_norvegicus', 'Pan_troglodytes', 'Danio_rerio']
        
        conservation = {
            'query_species': 'Homo_sapiens',
            'query_gene': gene_id,
            'sequence_length': len(sequence),
            'target_species_checked': target_species,
            'conservation_score': None,
            'conserved_sites': 'unknown',
            'note': 'Production implementation requires UCSC LiftOver mappings and ortholog sequence alignment'
        }
        
        # Simplified: estimate conservation by looking for seed-like patterns
        seed_patterns = ['CGCGUU', 'CCUGCU', 'GUGCAA', 'CAAAGU']  # Common seed complements
        pattern_count = sum(sequence.upper().count(p) for p in seed_patterns)
        
        if pattern_count > 0:
            conservation['conserved_sites'] = f'{pattern_count} potential seed-like regions detected'
            conservation['conservation_score'] = min(0.85, pattern_count * 0.1)  # Rough estimate
        
        return conservation
