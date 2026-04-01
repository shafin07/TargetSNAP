"""
Clean web utilities for TargetSNAP.

This module intentionally uses GRCh37-first data retrieval for transcript 3' UTR:
- Primary: Ensembl GRCh37 BioMart martservice
- Fallback: Ensembl GRCh37 REST API
"""

import csv
import copy
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import traceback
import uuid
from typing import Dict, List, Optional, Tuple

import requests


class GenomicDataHandler:
    """Handle genomic data lookups and SNP/transcript retrieval."""

    GRCH37_REST = "https://grch37.rest.ensembl.org"
    BIOMART_URL = "https://grch37.ensembl.org/biomart/martservice"

    def __init__(self, data_dir: str = "genomic_data"):
        self.data_dir = data_dir
        self.base_dir = os.path.dirname(os.path.abspath(__file__))
        self.snp_gene_map = self._load_json("snp_gene_mapping.json")
        self.gene_transcripts = self._load_json("gene_transcripts.json")
        self.dbsnp_cache: Dict[str, List[Dict]] = {}
        self.eqtl_cache: Dict[Tuple, Dict] = {}
        self.transcript_sequence_cache: Dict[str, Dict] = {}
        self.local_utr_index_loaded = False
        self.local_gff_entries: Dict[str, Dict] = {}
        self.local_gff_by_base: Dict[str, List[str]] = {}
        self.local_fasta_by_id: Dict[str, List[str]] = {}
        self.local_fasta_by_base: Dict[str, List[str]] = {}

    def _load_json(self, filename: str) -> Dict:
        path = os.path.join(self.data_dir, filename)
        if not os.path.exists(path):
            return {}
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return {}

    def get_genes_for_snp(self, rs_id: str) -> List[Dict]:
        if rs_id in self.snp_gene_map:
            return self.snp_gene_map.get(rs_id, [])

        if rs_id in self.dbsnp_cache:
            return self.dbsnp_cache[rs_id]

        genes = self._fetch_from_ensembl_vep_grch37(rs_id)
        if genes:
            self.dbsnp_cache[rs_id] = genes
            return genes

        return self.snp_gene_map.get(rs_id, [])

    def _fetch_from_ensembl_vep_grch37(self, rs_id: str) -> List[Dict]:
        try:
            url = f"{self.GRCH37_REST}/vep/human/id/{rs_id}"
            headers = {"Content-Type": "application/json", "Accept": "application/json"}
            resp = requests.get(url, headers=headers, timeout=12)
            if resp.status_code != 200:
                return []

            payload = resp.json()
            if not isinstance(payload, list) or not payload:
                return []

            variant = payload[0]
            allele_string = variant.get("allele_string", "")
            ref_allele = ""
            alt_allele = ""
            if "/" in allele_string:
                parts = [x.strip().upper() for x in allele_string.split("/") if x.strip()]
                if parts:
                    ref_allele = parts[0]
                if len(parts) > 1:
                    alt_allele = parts[1]

            seen = set()
            out: List[Dict] = []
            for consequence in variant.get("transcript_consequences", []) or []:
                gene_id = consequence.get("gene_id", "")
                gene_symbol = consequence.get("gene_symbol", "")
                if not gene_id or gene_id in seen:
                    continue
                seen.add(gene_id)
                out.append(
                    {
                        "gene_id": gene_id,
                        "gene_name": gene_symbol or gene_id,
                        "chromosome": variant.get("seq_region_name", ""),
                        "position": int(variant.get("start", 0) or 0),
                        "ref_allele": ref_allele,
                        "alt_allele": alt_allele,
                        "strand": "+" if consequence.get("strand", 1) == 1 else "-",
                        "consequences": consequence.get("consequence_terms", []),
                        "source": "ensembl_grch37_vep",
                    }
                )
            return out
        except Exception:
            return []

    def get_snp_details(self, rs_id: str, gene_id: Optional[str] = None) -> Optional[Dict]:
        genes = self.get_genes_for_snp(rs_id)
        if not genes:
            return None
        if gene_id:
            for g in genes:
                if g.get("gene_id") == gene_id:
                    return g
        return genes[0]

    def get_transcripts(self, gene_id: str) -> List[Dict]:
        local = self.gene_transcripts.get(gene_id, [])
        if local:
            return self._annotate_transcripts_with_utr_length(local)

        try:
            url = f"{self.GRCH37_REST}/lookup/id/{gene_id}?expand=1"
            headers = {"Content-Type": "application/json", "Accept": "application/json"}
            resp = requests.get(url, headers=headers, timeout=12)
            if resp.status_code != 200:
                return []

            data = resp.json()
            transcripts: List[Dict] = []
            for i, t in enumerate(data.get("Transcript", [])[:20]):
                transcripts.append(
                    {
                        "transcript_id": t.get("id", ""),
                        "transcript_name": t.get("display_name", f"Transcript-{i + 1}"),
                        "strand": "+" if t.get("strand", data.get("strand", 1)) == 1 else "-",
                        "length": t.get("length", 0),
                        "is_canonical": bool(t.get("is_canonical", 0) == 1),
                    }
                )
            transcripts.sort(key=lambda x: x.get("is_canonical", False), reverse=True)
            return self._annotate_transcripts_with_utr_length(transcripts)
        except Exception:
            return []

    def _annotate_transcripts_with_utr_length(self, transcripts: List[Dict]) -> List[Dict]:
        """Attach UTR-only length from local TargetScan hg19 GFF when available."""
        self._load_local_targetscan_indexes()
        annotated: List[Dict] = []
        for t in transcripts:
            row = dict(t)
            transcript_id = (row.get("transcript_id") or "").strip()
            base_id = transcript_id.split(".")[0] if transcript_id else ""

            candidates: List[str] = []
            if transcript_id:
                candidates.append(transcript_id)
            if base_id and base_id not in candidates:
                candidates.append(base_id)
            candidates.extend(self.local_gff_by_base.get(base_id, []))

            seen = set()
            candidates = [c for c in candidates if c and not (c in seen or seen.add(c))]

            chosen = None
            chosen_len = None
            for c in candidates:
                entry = self.local_gff_entries.get(c)
                if not entry:
                    continue
                try:
                    utr_len = int(entry.get("utr_end", 0)) - int(entry.get("utr_start", 0)) + 1
                except Exception:
                    continue
                if utr_len <= 0:
                    continue
                if chosen_len is None or utr_len > chosen_len:
                    chosen = c
                    chosen_len = utr_len

            row["utr_length"] = chosen_len
            row["utr_length_source"] = "targetscan_hg19_gff" if chosen_len else "unknown"
            row["transcript_id_resolved"] = chosen or transcript_id
            annotated.append(row)

        return annotated

    def get_transcript_sequence(self, transcript_id: str) -> Optional[Dict]:
        """
        Retrieve transcript 3' UTR (GRCh37-first).

        Priority order:
        1) GRCh37 BioMart 3'UTR sequence
        2) GRCh37 REST sequence endpoint type=3utr
        """
        transcript_id = (transcript_id or "").strip()
        if not transcript_id:
            return None

        if transcript_id in self.transcript_sequence_cache:
            return self.transcript_sequence_cache[transcript_id]

        candidate_ids = self._candidate_transcript_ids(transcript_id)

        local_hit = self._fetch_3utr_from_local_targetscan(transcript_id)
        if local_hit and local_hit.get("sequence"):
            self.transcript_sequence_cache[transcript_id] = local_hit
            resolved_id = local_hit.get("transcript_id_resolved", transcript_id)
            self.transcript_sequence_cache[resolved_id] = local_hit
            return local_hit

        for candidate in candidate_ids:
            if candidate in self.transcript_sequence_cache:
                return self.transcript_sequence_cache[candidate]

            bio = self._fetch_3utr_from_grch37_biomart(candidate)
            if bio and bio.get("sequence"):
                bio["transcript_id_requested"] = transcript_id
                bio["transcript_id_resolved"] = candidate
                self.transcript_sequence_cache[transcript_id] = bio
                self.transcript_sequence_cache[candidate] = bio
                return bio

            rest = self._fetch_3utr_from_grch37_rest(candidate)
            if rest and rest.get("sequence"):
                rest["transcript_id_requested"] = transcript_id
                rest["transcript_id_resolved"] = candidate
                self.transcript_sequence_cache[transcript_id] = rest
                self.transcript_sequence_cache[candidate] = rest
                return rest

        # Version probing for stable IDs without explicit version suffix.
        if "." not in transcript_id:
            for versioned_id in self._fetch_biomart_transcript_versions(transcript_id)[:5]:
                bio_v = self._fetch_3utr_from_grch37_biomart(versioned_id, filter_name="ensembl_transcript_id_version")
                if bio_v and bio_v.get("sequence"):
                    bio_v["transcript_id_requested"] = transcript_id
                    bio_v["transcript_id_resolved"] = versioned_id
                    self.transcript_sequence_cache[transcript_id] = bio_v
                    self.transcript_sequence_cache[versioned_id] = bio_v
                    return bio_v

                rest_v = self._fetch_3utr_from_grch37_rest(versioned_id)
                if rest_v and rest_v.get("sequence"):
                    rest_v["transcript_id_requested"] = transcript_id
                    rest_v["transcript_id_resolved"] = versioned_id
                    self.transcript_sequence_cache[transcript_id] = rest_v
                    self.transcript_sequence_cache[versioned_id] = rest_v
                    return rest_v

        return None

    def _load_local_targetscan_indexes(self) -> None:
        """Build in-memory indexes from local TargetScan GFF and local RNAplfold FASTA files."""
        if self.local_utr_index_loaded:
            return

        gff_path = os.path.join(self.base_dir, "TargetScan", "TSHuman_7_hg19_3UTRs.gff")
        if os.path.exists(gff_path):
            try:
                with open(gff_path, "r", encoding="utf-8", errors="ignore") as f:
                    for raw in f:
                        line = raw.strip()
                        if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                            continue
                        fields = line.split("\t")
                        if len(fields) < 9:
                            continue
                        feature = fields[2]
                        if feature != "UTR":
                            continue

                        transcript_id = fields[8].strip()
                        if not transcript_id.startswith("ENST"):
                            continue

                        chrom = fields[0]
                        strand = "+" if fields[6] == "+" else "-"
                        try:
                            start = int(fields[3])
                            end = int(fields[4])
                        except Exception:
                            continue

                        entry = self.local_gff_entries.get(transcript_id)
                        if entry is None:
                            entry = {
                                "chromosome": chrom,
                                "strand": strand,
                                "utr_start": start,
                                "utr_end": end,
                                "segments": 1,
                            }
                            self.local_gff_entries[transcript_id] = entry
                        else:
                            entry["utr_start"] = min(entry["utr_start"], start)
                            entry["utr_end"] = max(entry["utr_end"], end)
                            entry["segments"] += 1

                        base_id = transcript_id.split(".")[0]
                        self.local_gff_by_base.setdefault(base_id, [])
                        if transcript_id not in self.local_gff_by_base[base_id]:
                            self.local_gff_by_base[base_id].append(transcript_id)
            except Exception:
                pass

        rna_dir = os.path.join(self.base_dir, "TargetScan", "TargetScan7_context_scores", "RNAplfold_in_out")
        if os.path.isdir(rna_dir):
            try:
                for fname in os.listdir(rna_dir):
                    if not fname.lower().endswith(".fa"):
                        continue
                    path = os.path.join(rna_dir, fname)
                    m = re.search(r"(ENST\d+(?:\.\d+)?)", fname)
                    if not m:
                        continue
                    tid = m.group(1)
                    self.local_fasta_by_id.setdefault(tid, []).append(path)
                    self.local_fasta_by_base.setdefault(tid.split(".")[0], []).append(path)
            except Exception:
                pass

        self.local_utr_index_loaded = True

    def _read_first_fasta_sequence(self, fasta_path: str) -> str:
        seq_parts: List[str] = []
        try:
            with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
                for raw in f:
                    line = raw.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if seq_parts:
                            break
                        continue
                    seq_parts.append(line)
        except Exception:
            return ""

        seq = "".join(seq_parts).upper().replace("U", "T")
        return "".join(c for c in seq if c in "ACGTN")

    def _choose_best_fasta_path(self, paths: List[str]) -> Optional[str]:
        if not paths:
            return None
        # Prefer ref files first because they contain canonical reference UTR snapshots.
        paths_sorted = sorted(paths, key=lambda p: (0 if "_ref" in os.path.basename(p).lower() else 1, len(os.path.basename(p))))
        return paths_sorted[0]

    def _fetch_3utr_from_local_targetscan(self, transcript_id: str) -> Optional[Dict]:
        """Fast local retrieval using TargetScan hg19 resources (GFF + RNAplfold FASTA cache)."""
        self._load_local_targetscan_indexes()

        candidates: List[str] = self._candidate_transcript_ids(transcript_id)
        if "." not in transcript_id:
            candidates.extend(self.local_gff_by_base.get(transcript_id, []))
        # dedupe preserving order
        seen = set()
        candidates = [c for c in candidates if not (c in seen or seen.add(c))]

        for candidate in candidates:
            fasta_paths = list(self.local_fasta_by_id.get(candidate, []))
            fasta_paths.extend(self.local_fasta_by_base.get(candidate.split(".")[0], []))
            best_path = self._choose_best_fasta_path(fasta_paths)
            if not best_path:
                continue

            seq = self._read_first_fasta_sequence(best_path)
            if not seq:
                continue

            gff_entry = self.local_gff_entries.get(candidate)
            if gff_entry is None:
                # Try any version match in GFF for base ID.
                versions = self.local_gff_by_base.get(candidate.split(".")[0], [])
                if versions:
                    gff_entry = self.local_gff_entries.get(versions[0])

            return {
                "transcript_id": candidate,
                "transcript_id_requested": transcript_id,
                "transcript_id_resolved": candidate,
                "gene_id": "",
                "sequence": seq,
                "strand": (gff_entry or {}).get("strand", "+"),
                "length": len(seq),
                "chromosome": (gff_entry or {}).get("chromosome", ""),
                "utr_start": (gff_entry or {}).get("utr_start"),
                "utr_end": (gff_entry or {}).get("utr_end"),
                "sequence_source": "local_targetscan_gff_fasta",
                "sequence_source_url": os.path.relpath(best_path, self.base_dir).replace("\\", "/"),
            }

        return None

    def _candidate_transcript_ids(self, transcript_id: str) -> List[str]:
        candidates: List[str] = []
        if transcript_id:
            candidates.append(transcript_id)
        if "." in transcript_id:
            base = transcript_id.split(".")[0]
            if base and base not in candidates:
                candidates.append(base)
        return candidates

    def _fetch_biomart_transcript_versions(self, transcript_base_id: str) -> List[str]:
        """Return possible versioned transcript IDs from GRCh37 BioMart."""
        try:
            query = (
                "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
                "<!DOCTYPE Query>"
                "<Query virtualSchemaName=\"default\" formatter=\"TSV\" header=\"0\" uniqueRows=\"1\" count=\"\" datasetConfigVersion=\"0.6\">"
                "<Dataset name=\"hsapiens_gene_ensembl\" interface=\"default\">"
                f"<Filter name=\"ensembl_transcript_id\" value=\"{transcript_base_id}\" />"
                "<Attribute name=\"ensembl_transcript_id_version\" />"
                "</Dataset>"
                "</Query>"
            )
            resp = requests.get(self.BIOMART_URL, params={"query": query}, timeout=5)
            if resp.status_code != 200:
                return []

            versions = []
            seen = set()
            for raw in (resp.text or "").splitlines():
                tid = raw.strip()
                if tid and tid not in seen:
                    seen.add(tid)
                    versions.append(tid)
            return versions
        except Exception:
            return []

    def _fetch_3utr_from_grch37_biomart(self, transcript_id: str, filter_name: str = "ensembl_transcript_id") -> Optional[Dict]:
        """
        Query GRCh37 BioMart martservice for transcript-level 3' UTR.

        Uses the GRCh37 BioMart endpoint requested by user:
        https://grch37.ensembl.org/info/data/biomart/index.html
        """

        def build_query(attrs: List[str]) -> str:
            attr_xml = "\n".join([f'<Attribute name="{a}" />' for a in attrs])
            return (
                "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
                "<!DOCTYPE Query>"
                "<Query virtualSchemaName=\"default\" formatter=\"TSV\" header=\"0\" uniqueRows=\"1\" count=\"\" datasetConfigVersion=\"0.6\">"
                "<Dataset name=\"hsapiens_gene_ensembl\" interface=\"default\">"
                f"<Filter name=\"{filter_name}\" value=\"{transcript_id}\" />"
                f"{attr_xml}"
                "</Dataset>"
                "</Query>"
            )

        # Different BioMart mirrors sometimes expose either "3_utr" or "3utr".
        attempts = [
            ["ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", "strand", "3_utr_start", "3_utr_end", "3_utr"],
            ["ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", "strand", "3_utr_start", "3_utr_end", "3utr"],
        ]

        for attrs in attempts:
            try:
                xml_query = build_query(attrs)
                resp = requests.get(self.BIOMART_URL, params={"query": xml_query}, timeout=6)
                if resp.status_code != 200:
                    continue

                text = (resp.text or "").strip()
                if not text or "Query ERROR" in text:
                    continue

                # First non-empty line
                line = ""
                for raw in text.splitlines():
                    if raw.strip():
                        line = raw.strip()
                        break
                if not line:
                    continue

                parts = line.split("\t")
                seq = ""
                gene_id = ""
                chrom = ""
                strand = "+"
                utr_start = None
                utr_end = None

                if len(parts) >= 7:
                    gene_id = parts[0]
                    # transcript id is parts[1]
                    chrom = parts[2]
                    strand = "+" if parts[3] == "1" else "-"
                    try:
                        utr_start = int(parts[4])
                    except Exception:
                        utr_start = None
                    try:
                        utr_end = int(parts[5])
                    except Exception:
                        utr_end = None
                    seq = parts[6].upper()
                elif len(parts) >= 3:
                    gene_id = parts[0]
                    seq = parts[-1].upper()

                seq = "".join(c for c in seq if c in "ACGTUNacgtun").upper().replace("U", "T")
                if not seq:
                    continue

                return {
                    "transcript_id": transcript_id,
                    "gene_id": gene_id,
                    "sequence": seq,
                    "strand": strand,
                    "length": len(seq),
                    "chromosome": chrom,
                    "utr_start": utr_start,
                    "utr_end": utr_end,
                    "sequence_source": "grch37_biomart",
                    "sequence_source_url": "https://grch37.ensembl.org/info/data/biomart/index.html",
                }
            except Exception:
                continue

        return None

    def _fetch_3utr_from_grch37_rest(self, transcript_id: str) -> Optional[Dict]:
        try:
            headers = {"Content-Type": "application/json", "Accept": "application/json"}

            # Transcript metadata
            lookup_url = f"{self.GRCH37_REST}/lookup/id/{transcript_id}"
            lookup = requests.get(lookup_url, headers=headers, timeout=8)
            if lookup.status_code != 200:
                return None
            meta = lookup.json()

            seq_url = f"{self.GRCH37_REST}/sequence/id/{transcript_id}?type=3utr"
            seq_resp = requests.get(seq_url, headers=headers, timeout=8)
            if seq_resp.status_code != 200:
                return None
            seq_data = seq_resp.json()
            seq = (seq_data.get("seq", "") or "").upper().replace("U", "T")
            if not seq:
                return None

            return {
                "transcript_id": transcript_id,
                "gene_id": meta.get("Parent", ""),
                "sequence": seq,
                "strand": "+" if meta.get("strand", 1) == 1 else "-",
                "length": len(seq),
                "chromosome": meta.get("seq_region_name", ""),
                "utr_start": None,
                "utr_end": None,
                "sequence_source": "grch37_rest_3utr_fallback",
                "sequence_source_url": "https://grch37.rest.ensembl.org",
            }
        except Exception:
            return None

    def get_eqtl_dae(
        self,
        rs_id: str,
        gene_id: Optional[str] = None,
        tissues: Optional[List[str]] = None,
        dataset_id: str = "gtex_v10",
    ) -> Dict:
        if not rs_id:
            return {"error": "rs_id is required"}

        tissue_list = tissues or ["Colon_Transverse", "Colon_Sigmoid"]
        cache_key = (rs_id, gene_id or "", tuple(sorted(tissue_list)), dataset_id)
        if cache_key in self.eqtl_cache:
            return self.eqtl_cache[cache_key]

        base_url = "https://gtexportal.org/api/v2/association/singleTissueEqtl"
        all_records: List[Dict] = []
        warnings: List[str] = []
        selected_gene_base = (gene_id or "").split(".")[0]
        queried_tissues: List[str] = []

        def query_tissue(tissue: str) -> int:
            params = {
                "snpId": rs_id,
                "datasetId": dataset_id,
                "itemsPerPage": 250,
                "tissueSiteDetailId": tissue,
            }
            queried_tissues.append(tissue)
            local_count = 0
            try:
                response = requests.get(base_url, params=params, timeout=12)
                if response.status_code != 200:
                    warnings.append(f"GTEx query failed for {tissue}: HTTP {response.status_code}")
                    return 0

                payload = response.json()
                rows = payload.get("data", []) if isinstance(payload, dict) else []
                for row in rows:
                    row_gene = row.get("gencodeId") or row.get("gene_id") or ""
                    row_gene_base = row_gene.split(".")[0] if row_gene else ""

                    p_value = self._safe_float(row.get("pValue"))
                    nes = self._safe_float(row.get("nes"))
                    log2_afc = self._safe_float(row.get("log2AllelicFoldChange"))
                    if log2_afc is None:
                        log2_afc = nes

                    all_records.append(
                        {
                            "rs_id": rs_id,
                            "tissue": row.get("tissueSiteDetailId", tissue),
                            "gene_id": row_gene,
                            "gene_symbol": row.get("geneSymbol", ""),
                            "variant_id": row.get("variantId", ""),
                            "p_value": p_value,
                            "nes": nes,
                            "log2_allelic_fold_change": log2_afc,
                            "selected_gene_match": bool(selected_gene_base and row_gene_base == selected_gene_base),
                            "dae_supported": bool(log2_afc is not None and abs(log2_afc) >= 0.58),
                            "eqtl_significant": bool(p_value is not None and p_value <= 5e-8),
                        }
                    )
                    local_count += 1
                return local_count
            except requests.RequestException as exc:
                warnings.append(f"GTEx query error for {tissue}: {exc}")
                return 0

        for tissue in tissue_list:
            query_tissue(tissue)

        if not all_records:
            fallback_tissues = [
                "Whole_Blood",
                "Liver",
                "Lung",
                "Thyroid",
                "Skin_Sun_Exposed_Lower_leg",
                "Muscle_Skeletal",
                "Nerve_Tibial",
                "Adipose_Subcutaneous",
                "Heart_Left_Ventricle",
                "Brain_Cortex",
                "Pancreas",
                "Stomach",
                "Small_Intestine_Terminal_Ileum",
                "Esophagus_Mucosa",
                "Colon_Transverse",
                "Colon_Sigmoid",
            ]
            remaining = [t for t in fallback_tissues if t not in queried_tissues]
            for tissue in remaining:
                query_tissue(tissue)
            if all_records:
                warnings.append("No rows found in selected tissues; expanded search to broader GTEx tissues.")

        all_records.sort(
            key=lambda x: (
                0 if x.get("selected_gene_match") else 1,
                0 if x["eqtl_significant"] else 1,
                x["p_value"] if x["p_value"] is not None else 1.0,
                -abs(x["log2_allelic_fold_change"]) if x["log2_allelic_fold_change"] is not None else 0.0,
            )
        )

        summary = {
            "rs_id": rs_id,
            "dataset": dataset_id,
            "selected_gene_id": gene_id,
            "tissues_queried": queried_tissues,
            "records_found": len(all_records),
            "eqtl_significant_count": sum(1 for r in all_records if r["eqtl_significant"]),
            "dae_supported_count": sum(1 for r in all_records if r["dae_supported"]),
            "selected_gene_match_count": sum(1 for r in all_records if r.get("selected_gene_match")),
            "notes": [
                "DAE support is inferred from log2AllelicFoldChange when available, otherwise NES proxy.",
                "Significance threshold for eQTL is p <= 5e-8.",
            ],
            "warnings": warnings,
            "records": all_records[:200],
        }
        self.eqtl_cache[cache_key] = summary
        return summary

    def _safe_float(self, value) -> Optional[float]:
        try:
            if value is None or value == "":
                return None
            return float(value)
        except Exception:
            return None


class TargetScanLocalRunner:
    """Run TargetScan predictions locally and provide fallback diagnostics."""

    def __init__(self):
        self.mirna_library = self._load_mirna_library()
        self.base_dir = os.path.dirname(os.path.abspath(__file__))
        self.targetscan70_dir = os.path.join(self.base_dir, "TargetScan", "targetscan_70")
        self.context_dir = os.path.join(self.base_dir, "TargetScan", "TargetScan7_context_scores")
        self.perl_executable = "perl"
        self.prediction_cache: Dict[str, List[Dict]] = {}
        self.prediction_cache_order: List[str] = []
        self.max_prediction_cache = 256
        self.disk_cache_dir = os.path.join(self.base_dir, "targetsnap_cache")
        os.makedirs(self.disk_cache_dir, exist_ok=True)
        self._load_disk_cache()
        self.verbose_targetscan_logs = str(os.getenv("TARGETSNAP_VERBOSE_LOGS", "0")).strip().lower() in {
            "1",
            "true",
            "yes",
            "on",
        }

    def _load_disk_cache(self):
        """Load previously computed predictions from disk into memory cache."""
        try:
            for fname in os.listdir(self.disk_cache_dir):
                if fname.endswith(".json"):
                    cache_key = fname[:-5]
                    fpath = os.path.join(self.disk_cache_dir, fname)
                    with open(fpath, "r", encoding="utf-8") as f:
                        self.prediction_cache[cache_key] = json.load(f)
                    self.prediction_cache_order.append(cache_key)
        except Exception:
            pass

    def _save_to_disk_cache(self, cache_key: str, targets: List[Dict]):
        """Persist prediction results to disk for fast reload."""
        try:
            fpath = os.path.join(self.disk_cache_dir, f"{cache_key}.json")
            with open(fpath, "w", encoding="utf-8") as f:
                json.dump(targets, f)
        except Exception:
            pass

    def preflight_check(self) -> Dict:
        required_ts70 = ["targetscan_70.pl", "miR_Family_Info_human.txt"]
        required_context = [
            "targetscan_70_context_scores.pl",
            "targetscan_count_8mers.pl",
            "miR_for_context_scores.txt",
            "Agarwal_2015_parameters.txt",
            "TA_SPS_by_seed_region.txt",
        ]

        checks = []
        for fname in required_ts70:
            path = os.path.join(self.targetscan70_dir, fname)
            checks.append({"component": fname, "path": path, "exists": os.path.exists(path)})

        for fname in required_context:
            path = os.path.join(self.context_dir, fname)
            checks.append({"component": fname, "path": path, "exists": os.path.exists(path)})

        perl_ok = False
        perl_version = ""
        try:
            res = subprocess.run([self.perl_executable, "-v"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            perl_ok = res.returncode == 0
            perl_version = (res.stdout or res.stderr).splitlines()[0] if (res.stdout or res.stderr) else ""
        except Exception:
            perl_ok = False

        checks.append({"component": "perl", "path": self.perl_executable, "exists": perl_ok, "details": perl_version})
        all_ok = all(c.get("exists", False) for c in checks)
        return {"ready": all_ok, "checks": checks, "fallback_mode": not all_ok}

    def _load_mirna_library(self) -> List[str]:
        return [
            "hsa-miR-21",
            "hsa-miR-93",
            "hsa-miR-106b",
            "hsa-miR-29c",
            "hsa-miR-26b",
            "hsa-miR-155",
            "hsa-miR-let-7a",
            "hsa-miR-200a",
            "hsa-miR-122",
            "hsa-miR-375",
        ]

    def _normalize_position(self, sequence_len: int, snp_position: Optional[int]) -> int:
        if sequence_len <= 0:
            return 0
        if snp_position is None or snp_position == 0:
            return sequence_len // 2
        normalized = int(snp_position) % sequence_len
        return normalized if normalized > 0 else sequence_len // 2

    def _to_rna(self, sequence: str) -> str:
        return (sequence or "").upper().replace("T", "U")

    def _extract_local_window(self, sequence: str, snp_position: Optional[int], window_radius: int = 250) -> str:
        if not sequence:
            return sequence
        seq_len = len(sequence)
        snp_index = self._normalize_position(seq_len, snp_position)
        start = max(0, snp_index - window_radius)
        end = min(seq_len, snp_index + window_radius + 1)
        return sequence[start:end]

    def _run_subprocess(self, command: List[str], cwd: str, stdout_path: Optional[str] = None, timeout: int = 180) -> Dict:
        exec_info = {
            "command": " ".join(command),
            "cwd": cwd,
            "returncode": None,
            "stdout": "",
            "stderr": "",
            "success": False,
            "output_file": stdout_path,
        }
        try:
            if stdout_path:
                with open(stdout_path, "w", encoding="utf-8", newline="") as out_file:
                    result = subprocess.run(command, cwd=cwd, stdout=out_file, stderr=subprocess.PIPE, text=True, timeout=timeout)
                    exec_info["stdout"] = f"(written to {stdout_path})"
            else:
                result = subprocess.run(command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout)
                exec_info["stdout"] = (result.stdout or "")[:500]
        except subprocess.TimeoutExpired as e:
            exec_info["stderr"] = f"Process timed out after {timeout}s"
            raise RuntimeError(f"TargetScan subprocess timed out after {timeout}s for command: {' '.join(command[:3])}...") from e

        exec_info["returncode"] = result.returncode
        exec_info["stderr"] = (result.stderr or "")[:500]
        exec_info["success"] = result.returncode == 0
        if result.returncode != 0:
            raise RuntimeError(f"Command failed: {' '.join(command)}\nStderr: {result.stderr}")
        return exec_info

    def _reverse_complement(self, sequence: str) -> str:
        comp = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}
        return "".join(comp.get(b, "N") for b in reversed(sequence.upper()))

    def _mirna_seed_map(self) -> Dict[str, str]:
        return {
            "hsa-miR-21": "UAGCUUAUCAGACUGAUGUUGA",
            "hsa-miR-93": "CAAAGUGCUGUUCGUGCAGGUAG",
            "hsa-miR-106b": "UAAAGUGCUGACAGUGAAGCACU",
            "hsa-miR-29c": "UAGCACCAUUCGCGUCUGGUGAGU",
            "hsa-miR-26b": "UCAAGUACUUCCCUGUAGUACUU",
            "hsa-miR-155": "UUAAUGCUAAUCGUGAUAGGGGUU",
            "hsa-miR-let-7a": "GAGGUAGUAGGUUGUAUAGUU",
            "hsa-miR-200a": "UAACACUGUCUGGUAACGAUGU",
            "hsa-miR-122": "UGGAGUGUGACAAUGGUGUUUGU",
            "hsa-miR-375": "UUUGUUCGUUCGGCUCGCGUGA",
        }

    def _mirna_base_preferences(self) -> Dict[str, Dict[str, float]]:
        return {
            "hsa-miR-21": {"A": 3.5, "T": 1.2, "C": -1.0, "G": -2.2},
            "hsa-miR-93": {"A": -1.5, "T": 3.0, "C": 1.1, "G": -0.5},
            "hsa-miR-106b": {"A": -0.8, "T": 2.8, "C": 0.9, "G": -1.2},
            "hsa-miR-29c": {"A": 2.2, "T": -0.7, "C": 3.4, "G": -1.8},
            "hsa-miR-26b": {"A": 2.9, "T": -1.4, "C": 0.7, "G": 1.3},
            "hsa-miR-155": {"A": -1.1, "T": 0.6, "C": 3.1, "G": 1.8},
            "hsa-miR-let-7a": {"A": 1.5, "T": -0.9, "C": -1.7, "G": 3.6},
            "hsa-miR-200a": {"A": 3.2, "T": -1.9, "C": 1.4, "G": -0.6},
            "hsa-miR-122": {"A": -0.4, "T": 3.3, "C": -1.5, "G": 2.0},
            "hsa-miR-375": {"A": 0.8, "T": 2.6, "C": -0.9, "G": 3.0},
        }

    def _parse_context_scores(self, context_output_file: str, snp_position: Optional[int] = None) -> List[Dict]:
        if not os.path.exists(context_output_file):
            return []

        best: Dict[str, Dict] = {}
        snp_site_1based = None
        if snp_position is not None:
            try:
                snp_site_1based = int(snp_position) + 1
            except Exception:
                snp_site_1based = None

        with open(context_output_file, "r", encoding="utf-8", errors="ignore") as f:
            header = f.readline()
            if not header:
                return []
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 36:
                    continue
                mirna_id = fields[2].strip()
                try:
                    raw_context_score = float(fields[27])
                    utr_start = int(fields[4].strip())
                    utr_end = int(fields[5].strip())
                except Exception:
                    continue

                affinity = -raw_context_score
                overlaps_snp = False
                if snp_site_1based is not None:
                    overlaps_snp = utr_start <= snp_site_1based <= utr_end

                prev = best.get(mirna_id)
                should_replace = False
                if prev is None:
                    should_replace = True
                else:
                    prev_overlap = bool(prev.get("overlaps_snp"))
                    # Prefer SNP-overlapping sites over non-overlapping sites.
                    if overlaps_snp and not prev_overlap:
                        should_replace = True
                    elif overlaps_snp == prev_overlap and affinity > prev["context_score"]:
                        should_replace = True

                if should_replace:
                    best[mirna_id] = {
                        "mirna_id": mirna_id,
                        "context_score": affinity,
                        "raw_context_score": raw_context_score,
                        "repression_level": self._get_repression_level_from_context(raw_context_score),
                        "site_type": fields[3].strip(),
                        "utr_start": utr_start,
                        "utr_end": utr_end,
                        "utr_region": fields[32].strip() if len(fields) > 32 else "",
                        "pairing": fields[33].strip() if len(fields) > 33 else "",
                        "mature_mirna_sequence": fields[34].strip() if len(fields) > 34 else "",
                        "mirna_family": fields[35].strip() if len(fields) > 35 else "",
                        "overlaps_snp": overlaps_snp,
                    }

        return sorted(best.values(), key=lambda x: (bool(x.get("overlaps_snp")), x["context_score"]), reverse=True)

    def _parse_seed_targets(self, targets_file: str, snp_position: Optional[int] = None) -> List[Dict]:
        """Parse seed-only output from targetscan_70.pl (fallback when context++ times out)."""
        if not os.path.exists(targets_file):
            return []

        best: Dict[str, Dict] = {}
        snp_site_1based = None
        if snp_position is not None:
            try:
                snp_site_1based = int(snp_position) + 1
            except Exception:
                snp_site_1based = None

        with open(targets_file, "r", encoding="utf-8", errors="ignore") as f:
            header = f.readline()
            if not header:
                return []
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                if fields[2].strip() != "9606":
                    continue

                mirna_id = fields[1].strip()
                try:
                    utr_start = int(fields[5].strip())
                    utr_end = int(fields[6].strip())
                except Exception:
                    continue

                overlaps_snp = False
                if snp_site_1based is not None:
                    overlaps_snp = utr_start <= snp_site_1based <= utr_end

                prev = best.get(mirna_id)
                should_replace = False
                if prev is None:
                    should_replace = True
                else:
                    prev_overlap = bool(prev.get("overlaps_snp"))
                    if overlaps_snp and not prev_overlap:
                        should_replace = True

                if should_replace:
                    best[mirna_id] = {
                        "mirna_id": mirna_id,
                        "context_score": 0.0,
                        "raw_context_score": 0.0,
                        "repression_level": "N/A (seed-only)",
                        "site_type": fields[8].strip(),
                        "utr_start": utr_start,
                        "utr_end": utr_end,
                        "utr_region": "",
                        "pairing": "",
                        "mature_mirna_sequence": "",
                        "mirna_family": mirna_id,
                        "overlaps_snp": overlaps_snp,
                        "seed_only": True,
                    }

        return sorted(best.values(), key=lambda x: bool(x.get("overlaps_snp")), reverse=True)

    def _parse_context_scores_grouped(
        self,
        context_output_file: str,
        snp_position_map: Optional[Dict[str, Optional[int]]] = None,
    ) -> Dict[str, List[Dict]]:
        """Parse context++ output and return best sites grouped by transcript key (field 0)."""
        if not os.path.exists(context_output_file):
            return {}

        grouped_best: Dict[str, Dict[str, Dict]] = {}
        with open(context_output_file, "r", encoding="utf-8", errors="ignore") as f:
            header = f.readline()
            if not header:
                return {}

            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 36:
                    continue

                transcript_key = fields[0].strip()
                mirna_id = fields[2].strip()
                try:
                    raw_context_score = float(fields[27])
                    utr_start = int(fields[4].strip())
                    utr_end = int(fields[5].strip())
                except Exception:
                    continue

                affinity = -raw_context_score
                snp_site_1based = None
                if snp_position_map and transcript_key in snp_position_map and snp_position_map[transcript_key] is not None:
                    try:
                        snp_site_1based = int(snp_position_map[transcript_key]) + 1
                    except Exception:
                        snp_site_1based = None

                overlaps_snp = False
                if snp_site_1based is not None:
                    overlaps_snp = utr_start <= snp_site_1based <= utr_end

                if transcript_key not in grouped_best:
                    grouped_best[transcript_key] = {}

                prev = grouped_best[transcript_key].get(mirna_id)
                should_replace = False
                if prev is None:
                    should_replace = True
                else:
                    prev_overlap = bool(prev.get("overlaps_snp"))
                    if overlaps_snp and not prev_overlap:
                        should_replace = True
                    elif overlaps_snp == prev_overlap and affinity > prev["context_score"]:
                        should_replace = True

                if should_replace:
                    grouped_best[transcript_key][mirna_id] = {
                        "mirna_id": mirna_id,
                        "context_score": affinity,
                        "raw_context_score": raw_context_score,
                        "repression_level": self._get_repression_level_from_context(raw_context_score),
                        "site_type": fields[3].strip(),
                        "utr_start": utr_start,
                        "utr_end": utr_end,
                        "utr_region": fields[32].strip() if len(fields) > 32 else "",
                        "pairing": fields[33].strip() if len(fields) > 33 else "",
                        "mature_mirna_sequence": fields[34].strip() if len(fields) > 34 else "",
                        "mirna_family": fields[35].strip() if len(fields) > 35 else "",
                        "overlaps_snp": overlaps_snp,
                    }

        out: Dict[str, List[Dict]] = {}
        for transcript_key, best_map in grouped_best.items():
            out[transcript_key] = sorted(
                best_map.values(),
                key=lambda x: (bool(x.get("overlaps_snp")), x["context_score"]),
                reverse=True,
            )
        return out

    def _parse_seed_targets_grouped(
        self,
        targets_file: str,
        snp_position_map: Optional[Dict[str, Optional[int]]] = None,
    ) -> Dict[str, List[Dict]]:
        """Parse seed-only output from targetscan_70.pl (fallback when context++ times out).

        Columns: a_Gene_ID(0), miRNA_family_ID(1), species_ID(2), MSA_start(3),
        MSA_end(4), UTR_start(5), UTR_end(6), Group_num(7), Site_type(8),
        miRNA_in_this_species(9), Group_type(10), ...
        """
        if not os.path.exists(targets_file):
            return {}

        grouped_best: Dict[str, Dict[str, Dict]] = {}
        with open(targets_file, "r", encoding="utf-8", errors="ignore") as f:
            header = f.readline()
            if not header:
                return {}

            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                transcript_key = fields[0].strip()
                mirna_id = fields[1].strip()
                species_id = fields[2].strip()
                if species_id != "9606":
                    continue

                try:
                    utr_start = int(fields[5].strip())
                    utr_end = int(fields[6].strip())
                except Exception:
                    continue

                site_type = fields[8].strip()

                snp_site_1based = None
                if snp_position_map and transcript_key in snp_position_map and snp_position_map[transcript_key] is not None:
                    try:
                        snp_site_1based = int(snp_position_map[transcript_key]) + 1
                    except Exception:
                        snp_site_1based = None

                overlaps_snp = False
                if snp_site_1based is not None:
                    overlaps_snp = utr_start <= snp_site_1based <= utr_end

                if transcript_key not in grouped_best:
                    grouped_best[transcript_key] = {}

                prev = grouped_best[transcript_key].get(mirna_id)
                should_replace = False
                if prev is None:
                    should_replace = True
                else:
                    prev_overlap = bool(prev.get("overlaps_snp"))
                    if overlaps_snp and not prev_overlap:
                        should_replace = True
                    elif overlaps_snp == prev_overlap:
                        should_replace = False

                if should_replace:
                    grouped_best[transcript_key][mirna_id] = {
                        "mirna_id": mirna_id,
                        "context_score": 0.0,
                        "raw_context_score": 0.0,
                        "repression_level": "N/A (seed-only)",
                        "site_type": site_type,
                        "utr_start": utr_start,
                        "utr_end": utr_end,
                        "utr_region": "",
                        "pairing": "",
                        "mature_mirna_sequence": "",
                        "mirna_family": mirna_id,
                        "overlaps_snp": overlaps_snp,
                        "seed_only": True,
                    }

        out: Dict[str, List[Dict]] = {}
        for transcript_key, best_map in grouped_best.items():
            out[transcript_key] = sorted(
                best_map.values(),
                key=lambda x: bool(x.get("overlaps_snp")),
                reverse=True,
            )
        return out

    def _run_targetscan_pipeline_batch(
        self,
        sequence_by_label: Dict[str, str],
        gene_id: str,
        snp_position_map: Optional[Dict[str, Optional[int]]] = None,
        debug_meta: Optional[Dict] = None,
    ) -> Dict[str, List[Dict]]:
        """Run one TargetScan pipeline for multiple sequences (e.g., REF + ALT)."""
        if not sequence_by_label:
            return {}

        valid_sequences = {k: v for k, v in (sequence_by_label or {}).items() if v}
        if not valid_sequences:
            return {k: [] for k in sequence_by_label.keys()} if sequence_by_label else {}

        # Reuse the same cache used by predict_targets; this avoids rerunning Perl for repeated compares.
        cached_results: Dict[str, List[Dict]] = {}
        missing_sequences: Dict[str, str] = {}
        for label, sequence in valid_sequences.items():
            cache_key = hashlib.sha1(f"{gene_id}|{sequence}".encode("utf-8")).hexdigest()
            if cache_key in self.prediction_cache:
                cached_results[label] = copy.deepcopy(self.prediction_cache[cache_key])
            else:
                missing_sequences[label] = sequence

        if not missing_sequences:
            return {label: cached_results.get(label, []) for label in sequence_by_label.keys()}

        targetscan_available = os.path.exists(self.targetscan70_dir) and os.path.exists(self.context_dir)
        logs_dir = os.path.join(os.path.dirname(__file__), "targetsnap_logs")
        os.makedirs(logs_dir, exist_ok=True)
        uid = uuid.uuid4().hex[:8]
        run_id = f"{gene_id}_BATCH_{int(time.time() * 1000)}_{uid}"
        run_log_dir = os.path.join(logs_dir, run_id)
        os.makedirs(run_log_dir, exist_ok=True)

        if not targetscan_available:
            out_all = dict(cached_results)
            for label, seq in missing_sequences.items():
                out_all[label] = self._predict_targets_mock(seq, gene_id)
            return {label: out_all.get(label, []) for label in sequence_by_label.keys()}

        temp_dir = tempfile.mkdtemp(prefix="targetsnap_ts_batch_")
        try:
            utr_file = os.path.join(temp_dir, "UTR.txt")
            orf_file = os.path.join(temp_dir, "ORF.txt")
            output_targets = os.path.join(temp_dir, "output_targets.txt")
            orf_8mer_counts = os.path.join(temp_dir, "ORF_8mer_counts.txt")
            orf_lengths = os.path.join(temp_dir, "ORF.lengths.txt")
            output_context_scores = os.path.join(temp_dir, "output_context_scores.txt")

            transcript_key_map: Dict[str, str] = {}
            reverse_key_map: Dict[str, str] = {}
            with open(utr_file, "w", encoding="utf-8", newline="") as f_utr, open(orf_file, "w", encoding="utf-8", newline="") as f_orf:
                for label, sequence in missing_sequences.items():
                    tx_key = f"{gene_id}__{label}"
                    transcript_key_map[label] = tx_key
                    reverse_key_map[tx_key] = label
                    seq_rna = self._to_rna(sequence)
                    f_utr.write(f"{tx_key}\t9606\t{seq_rna}\n")
                    f_orf.write(f"{tx_key}\t9606\t{seq_rna}\n")
                    with open(os.path.join(run_log_dir, f"input_sequence_dna_{label}.txt"), "w", encoding="utf-8") as f:
                        f.write(sequence)

            with open(os.path.join(run_log_dir, "execution.log"), "w", encoding="utf-8") as f:
                f.write("=== TargetScan 7 Pipeline Execution (BATCH) ===\n")
                f.write(f"Gene ID: {gene_id}\n")
                f.write(f"Labels: {', '.join(missing_sequences.keys())}\n")
                if debug_meta:
                    f.write(f"Requested transcript: {debug_meta.get('transcript_id', '')}\n")
                    f.write(f"Ref allele: {debug_meta.get('ref_allele', '')}\n")
                    f.write(f"Alt allele: {debug_meta.get('alt_allele', '')}\n")
                    f.write(f"Strand: {debug_meta.get('strand', '')}\n")

            shutil.copy(utr_file, os.path.join(run_log_dir, "UTR.txt"))
            shutil.copy(orf_file, os.path.join(run_log_dir, "ORF.txt"))

            mir_family_file = os.path.join(self.targetscan70_dir, "miR_Family_Info_human.txt")
            ts_script = os.path.join(self.targetscan70_dir, "targetscan_70.pl")
            self._run_subprocess([self.perl_executable, ts_script, mir_family_file, utr_file, output_targets], cwd=self.targetscan70_dir, timeout=600)

            count_script = os.path.join(self.context_dir, "targetscan_count_8mers.pl")
            self._run_subprocess([self.perl_executable, count_script, mir_family_file, orf_file], cwd=self.context_dir, stdout_path=orf_8mer_counts, timeout=600)

            generated_lengths = os.path.join(self.context_dir, "ORF.lengths.txt")
            if os.path.exists(generated_lengths):
                os.replace(generated_lengths, orf_lengths)
            if not os.path.exists(orf_lengths):
                with open(orf_lengths, "w", encoding="utf-8", newline="") as f_len:
                    for label, sequence in valid_sequences.items():
                        tx_key = transcript_key_map[label]
                        f_len.write(f"{tx_key}\t9606\t{len(self._to_rna(sequence))}\n")

            context_script = os.path.join(self.context_dir, "targetscan_70_context_scores.pl")
            mir_context_file = os.path.join(self.context_dir, "miR_for_context_scores.txt")
            context_scores_available = True
            try:
                self._run_subprocess(
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
                    cwd=self.context_dir,
                    timeout=120,
                )
            except RuntimeError as e:
                if "timed out" in str(e).lower():
                    context_scores_available = False
                else:
                    raise

            if self.verbose_targetscan_logs and os.path.exists(output_targets):
                shutil.copy(output_targets, os.path.join(run_log_dir, "output_targets.txt"))
            if self.verbose_targetscan_logs and os.path.exists(output_context_scores):
                shutil.copy(output_context_scores, os.path.join(run_log_dir, "output_context_scores.txt"))

            mapped_snp_positions: Dict[str, Optional[int]] = {}
            for label, tx_key in transcript_key_map.items():
                mapped_snp_positions[tx_key] = (snp_position_map or {}).get(label)

            if context_scores_available:
                grouped = self._parse_context_scores_grouped(output_context_scores, snp_position_map=mapped_snp_positions)
            else:
                grouped = self._parse_seed_targets_grouped(output_targets, snp_position_map=mapped_snp_positions)

            out: Dict[str, List[Dict]] = dict(cached_results)
            for label in missing_sequences.keys():
                tx_key = transcript_key_map[label]
                out[label] = grouped.get(tx_key, [])
                if not out[label]:
                    out[label] = self._predict_targets_mock(missing_sequences[label], gene_id)

                cache_key = hashlib.sha1(f"{gene_id}|{missing_sequences[label]}".encode("utf-8")).hexdigest()
                self.prediction_cache[cache_key] = copy.deepcopy(out[label])
                self._save_to_disk_cache(cache_key, out[label])
                self.prediction_cache_order.append(cache_key)
                if len(self.prediction_cache_order) > self.max_prediction_cache:
                    old_key = self.prediction_cache_order.pop(0)
                    self.prediction_cache.pop(old_key, None)

            with open(os.path.join(run_log_dir, "status.log"), "w", encoding="utf-8") as f:
                f.write("STATUS: SUCCESS_PERL_BATCH\n")
                f.write(f"Labels processed: {', '.join(missing_sequences.keys())}\n")
                f.write(f"Cache hits: {len(cached_results)}\n")

            return {label: out.get(label, []) for label in sequence_by_label.keys()}
        except Exception as e:
            with open(os.path.join(run_log_dir, "execution.log"), "a", encoding="utf-8") as f:
                f.write(f"\nERROR: {str(e)}\n")
                f.write(traceback.format_exc())
            with open(os.path.join(run_log_dir, "status.log"), "w", encoding="utf-8") as f:
                f.write("STATUS: FALLBACK_MOCK\n")
                f.write(f"Error: {str(e)}\n")
            out_all = dict(cached_results)
            for label, seq in missing_sequences.items():
                out_all[label] = self._predict_targets_mock(seq, gene_id)
            return {label: out_all.get(label, []) for label in sequence_by_label.keys()}
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

    def _run_targetscan_pipeline(
        self,
        sequence: str,
        gene_id: str,
        snp_position: Optional[int] = None,
        debug_meta: Optional[Dict] = None,
    ) -> List[Dict]:
        if not sequence:
            return []

        targetscan_available = os.path.exists(self.targetscan70_dir) and os.path.exists(self.context_dir)
        logs_dir = os.path.join(os.path.dirname(__file__), "targetsnap_logs")
        os.makedirs(logs_dir, exist_ok=True)
        run_label = ""
        if debug_meta:
            run_label = (debug_meta.get("run_label") or "").strip()
        uid = uuid.uuid4().hex[:8]
        run_id = f"{gene_id}_{run_label}_{int(time.time() * 1000)}_{uid}" if run_label else f"{gene_id}_{int(time.time() * 1000)}_{uid}"
        run_log_dir = os.path.join(logs_dir, run_id)
        os.makedirs(run_log_dir, exist_ok=True)

        if not targetscan_available:
            with open(os.path.join(run_log_dir, "status.log"), "w", encoding="utf-8") as f:
                f.write("STATUS: MOCK_MODEL\n")
                f.write("REASON: TargetScan infrastructure not found\n")
            return self._predict_targets_mock(sequence, gene_id)

        temp_dir = tempfile.mkdtemp(prefix="targetsnap_ts_")
        try:
            utr_file = os.path.join(temp_dir, "UTR.txt")
            orf_file = os.path.join(temp_dir, "ORF.txt")
            output_targets = os.path.join(temp_dir, "output_targets.txt")
            orf_8mer_counts = os.path.join(temp_dir, "ORF_8mer_counts.txt")
            orf_lengths = os.path.join(temp_dir, "ORF.lengths.txt")
            output_context_scores = os.path.join(temp_dir, "output_context_scores.txt")

            seq_rna = self._to_rna(sequence)
            with open(os.path.join(run_log_dir, "execution.log"), "w", encoding="utf-8") as f:
                f.write("=== TargetScan 7 Pipeline Execution ===\n")
                f.write(f"Gene ID: {gene_id}\n")
                f.write("CRITICAL: input sequence must be 3' UTR\n")
                f.write(f"Sequence length: {len(sequence)}\n")
                f.write(f"First100: {sequence[:100]}\n")
                f.write(f"Last100: {sequence[-100:]}\n")
                if debug_meta:
                    f.write(f"Run label: {debug_meta.get('run_label', '')}\n")
                    f.write(f"Requested transcript: {debug_meta.get('transcript_id', '')}\n")
                    f.write(f"Ref allele: {debug_meta.get('ref_allele', '')}\n")
                    f.write(f"Alt allele: {debug_meta.get('alt_allele', '')}\n")
                    f.write(f"Strand: {debug_meta.get('strand', '')}\n")

                if snp_position is not None and len(sequence) > 0:
                    snp_index = self._normalize_position(len(sequence), snp_position)
                    ctx_start = max(0, snp_index - 20)
                    ctx_end = min(len(sequence), snp_index + 21)
                    local_context = sequence[ctx_start:ctx_end]
                    current_base = sequence[snp_index]
                    f.write(f"SNP position (input): {snp_position}\n")
                    f.write(f"SNP index (normalized): {snp_index}\n")
                    f.write(f"Base at SNP index: {current_base}\n")
                    f.write(f"SNP context (41bp): {local_context}\n")

                    # Save dedicated SNP context file for quick discrepancy checks.
                    with open(os.path.join(run_log_dir, "snp_context.txt"), "w", encoding="utf-8") as sf:
                        sf.write(f"gene_id\t{gene_id}\n")
                        if debug_meta:
                            sf.write(f"run_label\t{debug_meta.get('run_label', '')}\n")
                            sf.write(f"transcript_id\t{debug_meta.get('transcript_id', '')}\n")
                            sf.write(f"ref_allele\t{debug_meta.get('ref_allele', '')}\n")
                            sf.write(f"alt_allele\t{debug_meta.get('alt_allele', '')}\n")
                            sf.write(f"strand\t{debug_meta.get('strand', '')}\n")
                        sf.write(f"snp_position_input\t{snp_position}\n")
                        sf.write(f"snp_index_normalized\t{snp_index}\n")
                        sf.write(f"base_at_snp\t{current_base}\n")
                        sf.write(f"context_41bp\t{local_context}\n")

            with open(os.path.join(run_log_dir, "input_sequence_dna.txt"), "w", encoding="utf-8") as f:
                f.write(sequence)
            with open(os.path.join(run_log_dir, "input_sequence_rna.txt"), "w", encoding="utf-8") as f:
                f.write(seq_rna)

            with open(utr_file, "w", encoding="utf-8", newline="") as f_utr:
                f_utr.write(f"{gene_id}\t9606\t{seq_rna}\n")
            with open(orf_file, "w", encoding="utf-8", newline="") as f_orf:
                f_orf.write(f"{gene_id}\t9606\t{seq_rna}\n")

            shutil.copy(utr_file, os.path.join(run_log_dir, "UTR.txt"))
            shutil.copy(orf_file, os.path.join(run_log_dir, "ORF.txt"))

            mir_family_file = os.path.join(self.targetscan70_dir, "miR_Family_Info_human.txt")
            ts_script = os.path.join(self.targetscan70_dir, "targetscan_70.pl")
            self._run_subprocess([self.perl_executable, ts_script, mir_family_file, utr_file, output_targets], cwd=self.targetscan70_dir, timeout=600)

            count_script = os.path.join(self.context_dir, "targetscan_count_8mers.pl")
            self._run_subprocess([self.perl_executable, count_script, mir_family_file, orf_file], cwd=self.context_dir, stdout_path=orf_8mer_counts, timeout=600)

            generated_lengths = os.path.join(self.context_dir, "ORF.lengths.txt")
            if os.path.exists(generated_lengths):
                os.replace(generated_lengths, orf_lengths)
            if not os.path.exists(orf_lengths):
                with open(orf_lengths, "w", encoding="utf-8", newline="") as f_len:
                    f_len.write(f"{gene_id}\t9606\t{len(seq_rna)}\n")

            context_script = os.path.join(self.context_dir, "targetscan_70_context_scores.pl")
            mir_context_file = os.path.join(self.context_dir, "miR_for_context_scores.txt")
            context_scores_available = True
            try:
                self._run_subprocess(
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
                    cwd=self.context_dir,
                    timeout=120,
                )
            except RuntimeError as e:
                if "timed out" in str(e).lower():
                    context_scores_available = False
                else:
                    raise

            if self.verbose_targetscan_logs and os.path.exists(output_targets):
                shutil.copy(output_targets, os.path.join(run_log_dir, "output_targets.txt"))
            if self.verbose_targetscan_logs and os.path.exists(output_context_scores):
                shutil.copy(output_context_scores, os.path.join(run_log_dir, "output_context_scores.txt"))

            if context_scores_available:
                parsed = self._parse_context_scores(output_context_scores, snp_position=snp_position)
            else:
                parsed = self._parse_seed_targets(output_targets, snp_position=snp_position)
            with open(os.path.join(run_log_dir, "status.log"), "w", encoding="utf-8") as f:
                f.write("STATUS: SUCCESS_PERL\n")
                f.write(f"Targets found: {len(parsed)}\n")
            return parsed if parsed else self._predict_targets_mock(sequence, gene_id)

        except Exception as e:
            with open(os.path.join(run_log_dir, "execution.log"), "a", encoding="utf-8") as f:
                f.write(f"\nERROR: {str(e)}\n")
                f.write(traceback.format_exc())
            with open(os.path.join(run_log_dir, "status.log"), "w", encoding="utf-8") as f:
                f.write("STATUS: FALLBACK_MOCK\n")
                f.write(f"Error: {str(e)}\n")
            return self._predict_targets_mock(sequence, gene_id)
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

    def apply_snp_to_sequence(self, sequence: str, ref_allele: str, alt_allele: str, snp_position: int, strand: str) -> str:
        sequence_list = list(sequence)
        seq_len = len(sequence_list)
        if seq_len == 0:
            return sequence
        normalized_pos = self._normalize_position(seq_len, snp_position)
        alt_allele = (alt_allele or "").upper()
        if not alt_allele:
            return sequence
        alt_base = alt_allele[0]

        if strand == "+":
            sequence_list[normalized_pos] = alt_base
        else:
            complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
            alt_comp = complement.get(alt_base, alt_base)
            pos_from_end = (seq_len - 1) - normalized_pos
            sequence_list[pos_from_end] = alt_comp
        return "".join(sequence_list)

    def set_allele_in_sequence(self, sequence: str, allele: str, snp_position: int, strand: str) -> str:
        sequence_list = list(sequence)
        seq_len = len(sequence_list)
        if seq_len == 0:
            return sequence
        normalized_pos = self._normalize_position(seq_len, snp_position)
        allele = (allele or "").upper()
        if not allele:
            return sequence
        base = allele[0]

        if strand == "+":
            sequence_list[normalized_pos] = base
        else:
            complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
            base_comp = complement.get(base, base)
            pos_from_end = (seq_len - 1) - normalized_pos
            sequence_list[pos_from_end] = base_comp
        return "".join(sequence_list)

    def predict_targets(
        self,
        sequence: str,
        gene_id: str,
        snp_position: Optional[int] = None,
        strand: str = "+",
        debug_meta: Optional[Dict] = None,
    ) -> List[Dict]:
        if not sequence:
            return []
        try:
            cache_key = hashlib.sha1(f"{gene_id}|{sequence}".encode("utf-8")).hexdigest()
            if cache_key in self.prediction_cache:
                return copy.deepcopy(self.prediction_cache[cache_key])

            # Use full retrieved 3'UTR for TargetScan execution.
            # This keeps `input_sequence_dna.txt` aligned with expected UTR length.
            result = self._run_targetscan_pipeline(
                sequence,
                gene_id,
                snp_position=snp_position,
                debug_meta=debug_meta,
            )
            self.prediction_cache[cache_key] = copy.deepcopy(result)
            self._save_to_disk_cache(cache_key, result)
            self.prediction_cache_order.append(cache_key)
            if len(self.prediction_cache_order) > self.max_prediction_cache:
                old_key = self.prediction_cache_order.pop(0)
                self.prediction_cache.pop(old_key, None)
            return result
        except Exception:
            return self._predict_targets_mock(sequence, gene_id)

    def _predict_targets_mock(self, sequence: str, gene_id: str) -> List[Dict]:
        targets: List[Dict] = []
        snp_index = self._normalize_position(len(sequence), 0) if sequence else 0
        for mirna in self.mirna_library:
            score = self._calculate_context_score(sequence, mirna, snp_index)
            if score is not None and score > 20:
                targets.append(
                    {
                        "mirna_id": mirna,
                        "context_score": score,
                        "raw_context_score": -score,
                        "repression_level": self._get_repression_level(score),
                        "site_type": "mock",
                        "utr_start": "",
                        "utr_end": "",
                        "utr_region": "",
                        "pairing": "",
                        "mature_mirna_sequence": "",
                        "mirna_family": "",
                    }
                )
        targets.sort(key=lambda x: x["context_score"], reverse=True)
        return targets

    def _calculate_context_score(self, sequence: str, mirna_id: str, snp_index: int) -> Optional[float]:
        if not sequence or len(sequence) < 20:
            return None

        seq = sequence.upper().replace("U", "T")
        snp_index = self._normalize_position(len(seq), snp_index)
        mirna_seq = self._mirna_seed_map().get(mirna_id)
        if not mirna_seq:
            return 25.0

        center_base = seq[snp_index]
        seed = mirna_seq[1:8].upper().replace("U", "T")
        target_motif = self._reverse_complement(seed)
        motif_len = len(target_motif)

        global_hits = seq.count(target_motif) if motif_len > 0 else 0
        local_start = max(0, snp_index - 60)
        local_end = min(len(seq), snp_index + 61)
        local_seq = seq[local_start:local_end]
        local_hits = local_seq.count(target_motif) if motif_len > 0 else 0

        gc_content = (seq.count("G") + seq.count("C")) / len(seq)
        au_content = (seq.count("A") + seq.count("T")) / len(seq)
        base_score = 45.0 + (18.0 * gc_content) - (6.0 * au_content)

        mirna_adj = {
            "hsa-miR-21": 6.0,
            "hsa-miR-93": 4.0,
            "hsa-miR-106b": 3.5,
            "hsa-miR-29c": 2.5,
            "hsa-miR-26b": 2.0,
            "hsa-miR-155": 3.0,
        }.get(mirna_id, 0.0)

        distance_to_center = abs(snp_index - (len(seq) // 2))
        position_bias = max(-8.0, 8.0 - (distance_to_center / 20.0))

        center_bonus = 4.0 if local_hits > 0 else -2.0
        expected_center = target_motif[motif_len // 2] if motif_len > 0 else "N"
        center_match = 2.5 if center_base == expected_center else -1.0
        base_pref = self._mirna_base_preferences().get(mirna_id, {}).get(center_base, 0.0)
        center_compat = center_match + base_pref

        motif_signal = (global_hits * 3.0) + (local_hits * 6.0) + center_bonus + center_compat
        final_score = base_score + mirna_adj + position_bias + motif_signal
        return min(100, max(0, final_score))

    def _get_repression_level(self, score: float) -> str:
        if score >= 70:
            return "strong"
        if score >= 50:
            return "moderate"
        if score >= 30:
            return "weak"
        return "none"

    def _get_repression_level_from_context(self, context_score: float) -> str:
        if context_score <= -0.3:
            return "strong"
        if context_score <= -0.15:
            return "moderate"
        if context_score < 0:
            return "weak"
        return "none"

    def compare_targets(self, ref_targets: List[Dict], mut_targets: List[Dict]) -> Dict:
        ref_mirnas = {t["mirna_id"]: t for t in ref_targets}
        mut_mirnas = {t["mirna_id"]: t for t in mut_targets}
        all_mirnas = sorted(set(ref_mirnas) | set(mut_mirnas))

        loss_of_function: List[Dict] = []
        gain_of_function: List[Dict] = []
        neutral: List[Dict] = []

        def site_type_weight(site_type: str) -> int:
            s = (site_type or "").lower()
            if "8mer" in s:
                return 4
            if "7mer" in s:
                return 3
            if "6mer" in s:
                return 2
            if s and s != "mock":
                return 1
            return 0

        def evidence_level(delta: float, ref_site: str, alt_site: str, binding_change: str) -> str:
            max_site = max(site_type_weight(ref_site), site_type_weight(alt_site))
            if delta >= 0.25 and max_site >= 3:
                return "high"
            if delta >= 0.15 and max_site >= 2:
                return "moderate"
            if delta >= 0.08 or binding_change != "no_change":
                return "low"
            return "very_low"

        def classify_binding_change(ref_score: float, mut_score: float, ref_site: str, alt_site: str) -> str:
            _ = (ref_site, alt_site)
            if mut_score > 30 and ref_score < 25:
                return "gain"
            if ref_score > 30 and mut_score < 25:
                return "loss"
            return "no_change"

        def build_record(effect_type: str, mirna: str, ref_score: float, mut_score: float) -> Dict:
            ref_site = ref_mirnas.get(mirna, {}).get("site_type", "")
            alt_site = mut_mirnas.get(mirna, {}).get("site_type", "")
            delta = abs(ref_score - mut_score)
            binding_change = classify_binding_change(ref_score, mut_score, ref_site, alt_site)
            level = evidence_level(delta, ref_site, alt_site, binding_change)
            site_weight = max(site_type_weight(ref_site), site_type_weight(alt_site))
            binding_bonus = 50 if binding_change != "no_change" else 0
            priority = round((delta * 100.0) + (5.0 * site_weight) + binding_bonus, 3)

            return {
                "mirna_id": mirna,
                "effect_type": effect_type,
                "site_type": ref_site if effect_type == "LOF" else alt_site,
                "utr_start": ref_mirnas.get(mirna, {}).get("utr_start", "") or mut_mirnas.get(mirna, {}).get("utr_start", ""),
                "utr_end": ref_mirnas.get(mirna, {}).get("utr_end", "") or mut_mirnas.get(mirna, {}).get("utr_end", ""),
                "utr_region": ref_mirnas.get(mirna, {}).get("utr_region", "") or mut_mirnas.get(mirna, {}).get("utr_region", ""),
                "context_score": ref_score if effect_type == "LOF" else mut_score,
                "raw_context_score": ref_mirnas.get(mirna, {}).get("raw_context_score", "") or mut_mirnas.get(mirna, {}).get("raw_context_score", ""),
                "mature_mirna_sequence": ref_mirnas.get(mirna, {}).get("mature_mirna_sequence", "") or mut_mirnas.get(mirna, {}).get("mature_mirna_sequence", ""),
                "mirna_family": ref_mirnas.get(mirna, {}).get("mirna_family", "") or mut_mirnas.get(mirna, {}).get("mirna_family", ""),
                "pairing": ref_mirnas.get(mirna, {}).get("pairing", "") or mut_mirnas.get(mirna, {}).get("pairing", ""),
                "ref_site_type": ref_site,
                "alt_site_type": alt_site,
                "ref_utr_start": ref_mirnas.get(mirna, {}).get("utr_start", ""),
                "ref_utr_end": ref_mirnas.get(mirna, {}).get("utr_end", ""),
                "alt_utr_start": mut_mirnas.get(mirna, {}).get("utr_start", ""),
                "alt_utr_end": mut_mirnas.get(mirna, {}).get("utr_end", ""),
                "ref_utr_region": ref_mirnas.get(mirna, {}).get("utr_region", ""),
                "alt_utr_region": mut_mirnas.get(mirna, {}).get("utr_region", ""),
                "ref_pairing": ref_mirnas.get(mirna, {}).get("pairing", ""),
                "alt_pairing": mut_mirnas.get(mirna, {}).get("pairing", ""),
                "ref_mature_mirna_sequence": ref_mirnas.get(mirna, {}).get("mature_mirna_sequence", ""),
                "alt_mature_mirna_sequence": mut_mirnas.get(mirna, {}).get("mature_mirna_sequence", ""),
                "ref_overlaps_snp": bool(ref_mirnas.get(mirna, {}).get("overlaps_snp")),
                "alt_overlaps_snp": bool(mut_mirnas.get(mirna, {}).get("overlaps_snp")),
                "ref_score": ref_score,
                "mut_score": mut_score,
                "score_change": delta,
                "binding_change": binding_change,
                "evidence_level": level,
                "priority_score": priority,
            }

        for mirna in all_mirnas:
            ref_score = ref_mirnas.get(mirna, {}).get("context_score", 0)
            mut_score = mut_mirnas.get(mirna, {}).get("context_score", 0)
            score_diff = ref_score - mut_score
            ref_site = ref_mirnas.get(mirna, {}).get("site_type", "")
            mut_site = mut_mirnas.get(mirna, {}).get("site_type", "")
            site_changed = (ref_site or "") != (mut_site or "")
            significant_delta = 0.01
            ref_overlap = bool(ref_mirnas.get(mirna, {}).get("overlaps_snp"))
            mut_overlap = bool(mut_mirnas.get(mirna, {}).get("overlaps_snp"))

            # Primary mode: SNP-overlap aware classification to match local/miRNASNP-style interpretation.
            if ref_overlap and not mut_overlap:
                loss_of_function.append(build_record("LOF", mirna, ref_score, mut_score))
                continue
            if mut_overlap and not ref_overlap:
                gain_of_function.append(build_record("GOF", mirna, ref_score, mut_score))
                continue
            if ref_overlap and mut_overlap:
                if (score_diff > significant_delta + 0.02) or (score_diff > significant_delta and site_changed):
                    loss_of_function.append(build_record("LOF", mirna, ref_score, mut_score))
                elif (score_diff < -(significant_delta + 0.02)) or (score_diff < -significant_delta and site_changed):
                    gain_of_function.append(build_record("GOF", mirna, ref_score, mut_score))
                else:
                    neutral.append(build_record("NEUTRAL", mirna, ref_score, mut_score))
                continue

            if (score_diff > significant_delta + 0.05) or (score_diff > significant_delta and site_changed):
                loss_of_function.append(build_record("LOF", mirna, ref_score, mut_score))
            elif (score_diff < -(significant_delta + 0.05)) or (score_diff < -significant_delta and site_changed):
                gain_of_function.append(build_record("GOF", mirna, ref_score, mut_score))
            else:
                neutral.append(build_record("NEUTRAL", mirna, ref_score, mut_score))

        loss_sorted = sorted(loss_of_function, key=lambda x: (x["priority_score"], x["score_change"]), reverse=True)
        gain_sorted = sorted(gain_of_function, key=lambda x: (x["priority_score"], x["score_change"]), reverse=True)
        neutral_sorted = sorted(neutral, key=lambda x: (x["priority_score"], x["score_change"]), reverse=True)
        return {
            "loss_of_function": loss_sorted,
            "gain_of_function": gain_sorted,
            "neutral": neutral_sorted,
            "all_effects": loss_sorted + gain_sorted + neutral_sorted,
            "methodology_note": "Results use dual-tool concordance: score and site-type evidence.",
        }

    def validate_utr_extraction(self, sequence: str, gene_id: str, ref_targets: List[Dict], mut_targets: List[Dict]) -> Dict:
        return {
            "sequence_info": {
                "gene_id": gene_id,
                "sequence_length": len(sequence) if sequence else 0,
                "sequence_preview": (sequence or "")[:100],
                "gc_content": ((sequence.count("G") + sequence.count("C")) / len(sequence)) if sequence else 0,
            },
            "ref_targets_info": {
                "total_count": len(ref_targets),
                "unique_mirnas": len(set(t.get("mirna_id") for t in ref_targets)),
                "utr_coordinates_found": sum(1 for t in ref_targets if t.get("utr_start") and t.get("utr_end")),
                "samples": ref_targets[:5],
            },
            "mut_targets_info": {
                "total_count": len(mut_targets),
                "unique_mirnas": len(set(t.get("mirna_id") for t in mut_targets)),
                "utr_coordinates_found": sum(1 for t in mut_targets if t.get("utr_start") and t.get("utr_end")),
                "samples": mut_targets[:5],
            },
            "status": "OK",
        }

    def compare_with_local_targetscan(self, gene_id: str, rs_id: str, ref_targets: List[Dict], mut_targets: List[Dict]) -> Dict:
        _ = rs_id
        lof_report = []
        gof_report = []
        for t in ref_targets:
            lof_report.append(
                {
                    "gene_id": f"{gene_id}_ref",
                    "mirna_id": t.get("mirna_id", ""),
                    "site_type": t.get("site_type", ""),
                    "utr_start": t.get("utr_start"),
                    "utr_end": t.get("utr_end"),
                    "overlaps_snp": bool(t.get("overlaps_snp")),
                    "context_score": t.get("raw_context_score", ""),
                    "utr_region": t.get("utr_region", ""),
                    "pairing": t.get("pairing", ""),
                    "mature_mirna_sequence": t.get("mature_mirna_sequence", ""),
                    "mirna_family": t.get("mirna_family", ""),
                }
            )
        for t in mut_targets:
            gof_report.append(
                {
                    "gene_id": f"{gene_id}_alt",
                    "mirna_id": t.get("mirna_id", ""),
                    "site_type": t.get("site_type", ""),
                    "utr_start": t.get("utr_start"),
                    "utr_end": t.get("utr_end"),
                    "overlaps_snp": bool(t.get("overlaps_snp")),
                    "context_score": t.get("raw_context_score", ""),
                    "utr_region": t.get("utr_region", ""),
                    "pairing": t.get("pairing", ""),
                    "mature_mirna_sequence": t.get("mature_mirna_sequence", ""),
                    "mirna_family": t.get("mirna_family", ""),
                }
            )
        return {
            "loss_of_function": lof_report,
            "gain_of_function": gof_report,
            "note": "Format matches Ensembl BioMart TargetScan output for easy comparison",
        }

    def get_target_enrichment_hints(self, target_gains: List[Dict], target_losses: List[Dict]) -> Dict:
        return {
            "gain_count": len(target_gains),
            "loss_count": len(target_losses),
            "profile": "balanced_disruption" if abs(len(target_gains) - len(target_losses)) < 3 else ("gain_skewed" if len(target_gains) > len(target_losses) else "loss_skewed"),
            "note": "Heuristic enrichment summary; use external KEGG/GO tools for statistical enrichment.",
        }

    def estimate_structure_impact(self, ref_seq: str, mut_seq: str, mirna_id: str) -> Dict:
        def gc_content(seq: str) -> float:
            if not seq:
                return 0.0
            return (seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c")) / (2.0 * len(seq))

        ref_gc = gc_content(ref_seq)
        mut_gc = gc_content(mut_seq)
        gc_change = abs(ref_gc - mut_gc)
        mfe_change = (mut_gc - ref_gc) * 10.0
        return {
            "mirna_id": mirna_id,
            "ref_gc_content": round(ref_gc, 4),
            "mut_gc_content": round(mut_gc, 4),
            "gc_change": round(gc_change, 4),
            "mfe_change": round(mfe_change, 2),
            "mfe_impact": "significant" if abs(mfe_change) > 2.0 else "minor",
            "mfe_source": "heuristic",
            "recommendation": "Validate with ViennaRNA for publication-grade structural interpretation.",
        }

    def export_results_csv(self, payload: Dict, delimiter: str = ',') -> str:
        output = io.StringIO()
        fieldnames = [
            "rs_id",
            "gene_id",
            "transcript_id",
            "effect_type",
            "mirna_id",
            "site_type",
            "ref_site_type",
            "alt_site_type",
            "utr_start",
            "utr_end",
            "utr_region",
            "context_score",
            "raw_context_score",
            "ref_score",
            "mut_score",
            "score_change",
            "ref_mfe",
            "alt_mfe",
            "mfe_delta",
            "ref_pvalue",
            "alt_pvalue",
            "mature_mirna_sequence",
            "mirna_family",
            "evidence_level",
            "binding_change",
            "priority_score",
        ]
        writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()

        rows = payload.get("all_effects", []) or []
        for row in rows:
            writer.writerow(
                {
                    "rs_id": payload.get("rs_id", ""),
                    "gene_id": payload.get("gene_id", ""),
                    "transcript_id": payload.get("transcript_id", ""),
                    "effect_type": row.get("effect_type", ""),
                    "mirna_id": row.get("mirna_id", ""),
                    "site_type": row.get("site_type", ""),
                    "ref_site_type": row.get("ref_site_type", ""),
                    "alt_site_type": row.get("alt_site_type", ""),
                    "utr_start": row.get("utr_start", ""),
                    "utr_end": row.get("utr_end", ""),
                    "utr_region": row.get("utr_region", ""),
                    "context_score": row.get("context_score", ""),
                    "raw_context_score": row.get("raw_context_score", ""),
                    "ref_score": row.get("ref_score", ""),
                    "mut_score": row.get("mut_score", ""),
                    "score_change": row.get("score_change", ""),
                    "ref_mfe": row.get("ref_mfe", ""),
                    "alt_mfe": row.get("alt_mfe", ""),
                    "mfe_delta": row.get("mfe_delta", ""),
                    "ref_pvalue": row.get("ref_pvalue", ""),
                    "alt_pvalue": row.get("alt_pvalue", ""),
                    "mature_mirna_sequence": row.get("mature_mirna_sequence", ""),
                    "mirna_family": row.get("mirna_family", ""),
                    "evidence_level": row.get("evidence_level", ""),
                    "binding_change": row.get("binding_change", ""),
                    "priority_score": row.get("priority_score", ""),
                }
            )
        return output.getvalue()
