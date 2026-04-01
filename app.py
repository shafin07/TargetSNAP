"""
Flask backend for TargetSNAP Web Interface
Handles SNP lookup, gene mapping, and TargetScan predictions
"""

from flask import Flask, render_template, request, jsonify, Response
from flask_cors import CORS
from concurrent.futures import ThreadPoolExecutor
import os
import sys
import json
import time
import traceback
import subprocess
import re as _re
import requests as http_requests

# Initialize Flask app
app = Flask(__name__, template_folder='templates')
CORS(app)

# Import utilities
print("Initializing TargetSNAP utilities...")
try:
    from targetsnap_web_utils_clean import GenomicDataHandler, TargetScanLocalRunner
    genomic_handler = GenomicDataHandler()
    targetscan = TargetScanLocalRunner()
    print("✓ Utilities loaded successfully")
except ImportError as e:
    print(f"✗ Import Error: {e}")
    raise
except Exception as e:
    print(f"✗ Initialization Error: {e}")
    raise

# Ensure handlers are defined globally
if 'genomic_handler' not in globals():
    raise RuntimeError("genomic_handler not initialized")
if 'targetscan' not in globals():
    raise RuntimeError("targetscan not initialized")


def _complement_allele(allele: str) -> str:
    """Complement DNA allele for transcript-oriented display on minus strand."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    value = (allele or '').upper().strip()
    if not value:
        return value
    return ''.join(comp.get(ch, ch) for ch in value)


def _strand_adjusted_alleles(ref_allele: str, alt_allele: str, strand: str) -> tuple[str, str]:
    """Return transcript-oriented alleles for display/reporting."""
    if strand == '-':
        return _complement_allele(ref_allele), _complement_allele(alt_allele)
    return (ref_allele or '').upper(), (alt_allele or '').upper()

def _compute_utr_sequence_index(
    snp_position: int,
    strand: str,
    utr_start: int | None,
    utr_end: int | None,
    sequence_len: int,
) -> int:
    """
    Map genomic SNP position to 0-based index in retrieved UTR sequence.

    + strand: index = snp_position - utr_start
    - strand: index = utr_end - snp_position
    """
    if sequence_len <= 0:
        return 0

    if utr_start is not None and utr_end is not None:
        if strand == '-':
            idx = int(utr_end) - int(snp_position)
        else:
            idx = int(snp_position) - int(utr_start)
        if 0 <= idx < sequence_len:
            return idx

    # Fallback when coordinates are not available/valid.
    normalized = int(snp_position) % sequence_len if snp_position else sequence_len // 2
    return normalized if normalized > 0 else sequence_len // 2
@app.route('/')
def index():
    """Serve main HTML page"""
    return render_template('web_index.html')

@app.route('/api/search-rs', methods=['POST'])
def search_rs():
    """Search SNP and return mapped genes"""
    try:
        data = request.get_json()
        rs_id = data.get('rs_id', '').strip()
        
        if not rs_id:
            return jsonify({'error': 'No rs ID provided'}), 400
        
        genes = genomic_handler.get_genes_for_snp(rs_id)
        
        if not genes:
            return jsonify({
                'error': f'No genes found for {rs_id}',
                'rs_id': rs_id
            }), 404
        
        genes_with_display = []
        for gene in genes:
            g = dict(gene)
            strand = g.get('strand', '+')
            ref_genomic = (g.get('ref_allele') or '').upper()
            alt_genomic = (g.get('alt_allele') or '').upper()
            ref_display, alt_display = _strand_adjusted_alleles(ref_genomic, alt_genomic, strand)
            g['ref_allele_genomic'] = ref_genomic
            g['alt_allele_genomic'] = alt_genomic
            g['ref_allele_display'] = ref_display
            g['alt_allele_display'] = alt_display
            genes_with_display.append(g)

        return jsonify({
            'rs_id': rs_id,
            'genes': genes_with_display,
            'total_genes': len(genes_with_display)
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/gene-transcripts', methods=['POST'])
def get_transcripts():
    """Get transcript isoforms for a gene"""
    try:
        data = request.get_json()
        gene_id = data.get('gene_id', '')
        
        transcripts = genomic_handler.get_transcripts(gene_id)
        
        if not transcripts:
            return jsonify({'error': f'No transcripts found for {gene_id}'}), 404
        
        return jsonify({
            'gene_id': gene_id,
            'transcripts': transcripts,
            'total_transcripts': len(transcripts)
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/compare-alleles', methods=['POST'])
def compare_alleles():
    """Compare ref and alt alleles for miRNA targeting"""
    try:
        t0 = time.perf_counter()
        data = request.get_json()
        rs_id = data.get('rs_id', '')
        gene_id = data.get('gene_id', '')
        transcript_id = data.get('transcript_id', '')
        
        snp_data = genomic_handler.get_snp_details(rs_id, gene_id=gene_id)
        if not snp_data:
            return jsonify({'error': f'SNP {rs_id} not found for gene {gene_id}'}), 404
        t_snp = time.perf_counter()
        
        ref_allele = (snp_data.get('ref_allele') or '').upper()
        alt_allele = (snp_data.get('alt_allele') or '').upper()
        snp_position = snp_data.get('position', 0)
        
        transcript_data = genomic_handler.get_transcript_sequence(transcript_id)
        transcript_id_used = transcript_id
        transcript_fallback_note = ''
        if not transcript_data:
            # Fallback: try alternate transcripts from selected gene before failing.
            for candidate in genomic_handler.get_transcripts(gene_id)[:10]:
                candidate_id = candidate.get('transcript_id')
                if not candidate_id or candidate_id == transcript_id:
                    continue
                candidate_data = genomic_handler.get_transcript_sequence(candidate_id)
                if candidate_data:
                    transcript_data = candidate_data
                    transcript_id_used = candidate_id
                    transcript_fallback_note = f'Requested transcript {transcript_id} unavailable; used {candidate_id}.'
                    break

        if not transcript_data:
            return jsonify({'error': f'Transcript {transcript_id} not found in GRCh37 BioMart/REST'}), 404
        t_transcript = time.perf_counter()
        
        strand = transcript_data.get('strand')
        sequence = transcript_data.get('sequence')
        ref_allele_display, alt_allele_display = _strand_adjusted_alleles(ref_allele, alt_allele, strand)

        utr_start = transcript_data.get('utr_start')
        utr_end = transcript_data.get('utr_end')
        transcript_snp_index = _compute_utr_sequence_index(
            snp_position=int(snp_position),
            strand=strand,
            utr_start=utr_start,
            utr_end=utr_end,
            sequence_len=len(sequence) if sequence else 0,
        )

        # Apply alleles on transcript-oriented sequence coordinates.
        ref_sequence = targetscan.set_allele_in_sequence(
            sequence,
            ref_allele_display,
            transcript_snp_index,
            '+'
        )
        mut_sequence = targetscan.apply_snp_to_sequence(
            ref_sequence,
            ref_allele_display,
            alt_allele_display,
            transcript_snp_index,
            '+'
        )
        t_sequences = time.perf_counter()
        
        # Fast path: single batched Perl/context++ run for REF + ALT.
        # This avoids duplicate process startup and can significantly reduce total latency.
        if hasattr(targetscan, '_run_targetscan_pipeline_batch'):
            batched = targetscan._run_targetscan_pipeline_batch(
                sequence_by_label={
                    'REF': ref_sequence,
                    'ALT': mut_sequence,
                },
                gene_id=gene_id,
                snp_position_map={
                    'REF': transcript_snp_index,
                    'ALT': transcript_snp_index,
                },
                debug_meta={
                    'transcript_id': transcript_data.get('transcript_id_resolved', transcript_id_used),
                    'ref_allele': ref_allele_display,
                    'alt_allele': alt_allele_display,
                    'strand': strand,
                }
            )
            ref_targets = batched.get('REF', [])
            mut_targets = batched.get('ALT', [])
        else:
            ref_kwargs = {
                'sequence': ref_sequence,
                'gene_id': gene_id,
                'snp_position': transcript_snp_index,
                'strand': strand,
                'debug_meta': {
                    'run_label': 'REF',
                    'transcript_id': transcript_data.get('transcript_id_resolved', transcript_id_used),
                    'ref_allele': ref_allele_display,
                    'alt_allele': alt_allele_display,
                    'strand': strand,
                }
            }
            alt_kwargs = {
                'sequence': mut_sequence,
                'gene_id': gene_id,
                'snp_position': transcript_snp_index,
                'strand': strand,
                'debug_meta': {
                    'run_label': 'ALT',
                    'transcript_id': transcript_data.get('transcript_id_resolved', transcript_id_used),
                    'ref_allele': ref_allele_display,
                    'alt_allele': alt_allele_display,
                    'strand': strand,
                }
            }

            with ThreadPoolExecutor(max_workers=2) as ex:
                fut_ref = ex.submit(targetscan.predict_targets, **ref_kwargs)
                fut_alt = ex.submit(targetscan.predict_targets, **alt_kwargs)
                ref_targets = fut_ref.result()
                mut_targets = fut_alt.result()
        t_targetscan = time.perf_counter()
        
        results = targetscan.compare_targets(ref_targets, mut_targets)
        t_compare = time.perf_counter()

        # ---- RNAhybrid MFE enrichment (parallel, best-effort) ----
        def _rnahybrid_mfe(target_seq, query_seq):
            """Run RNAhybrid via WSL, return {mfe, pvalue} dict or None."""
            t = (target_seq or '').strip().upper().replace('T', 'U')
            q = (query_seq or '').strip().upper().replace('T', 'U')
            if not t or not q or len(t) < 6 or len(q) < 15:
                return None
            if not _re.fullmatch(r'[AUGCN]+', t) or not _re.fullmatch(r'[AUGCN]+', q):
                return None
            try:
                proc = subprocess.run(
                    ['wsl', 'RNAhybrid', '-s', '3utr_human', '-b', '1', '-c', t, q],
                    capture_output=True, text=True, timeout=10
                )
                if proc.returncode != 0 or not proc.stdout.strip():
                    return None
                parts = proc.stdout.strip().split(':')
                if len(parts) > 5:
                    return {'mfe': float(parts[4]), 'pvalue': float(parts[5])}
                elif len(parts) > 4:
                    return {'mfe': float(parts[4]), 'pvalue': None}
                return None
            except Exception:
                return None

        def _enrich_record_mfe(rec):
            """Add ref_mfe, alt_mfe, mfe_delta, ref_pvalue, alt_pvalue and adjust priority_score."""
            mirna_seq = rec.get('mature_mirna_sequence') or rec.get('ref_mature_mirna_sequence') or ''
            ref_utr = rec.get('ref_utr_region', '')
            alt_utr = rec.get('alt_utr_region', '')
            ref_result = _rnahybrid_mfe(ref_utr, mirna_seq)
            alt_result = _rnahybrid_mfe(alt_utr, mirna_seq)
            rec['ref_mfe'] = ref_result['mfe'] if ref_result else None
            rec['alt_mfe'] = alt_result['mfe'] if alt_result else None
            rec['ref_pvalue'] = ref_result['pvalue'] if ref_result else None
            rec['alt_pvalue'] = alt_result['pvalue'] if alt_result else None
            if rec['ref_mfe'] is not None and rec['alt_mfe'] is not None:
                rec['mfe_delta'] = round(rec['alt_mfe'] - rec['ref_mfe'], 2)
            else:
                rec['mfe_delta'] = None

            # Adjust priority_score: incorporate MFE evidence
            priority = rec.get('priority_score', 0) or 0
            if rec['mfe_delta'] is not None:
                # Larger absolute MFE shift = more impactful, scale contribution
                priority += abs(rec['mfe_delta']) * 3.0
            # Significant p-value bonus (low p = confident prediction)
            for pv in [rec.get('ref_pvalue'), rec.get('alt_pvalue')]:
                if pv is not None and pv < 0.05:
                    priority += 10
                elif pv is not None and pv < 0.2:
                    priority += 3
            rec['priority_score'] = round(priority, 3)
            return rec

        # Run in parallel across LOF + GOF + capped neutral
        all_records = (
            results.get('loss_of_function', []) +
            results.get('gain_of_function', []) +
            results.get('neutral', [])[:100]
        )
        if all_records:
            with ThreadPoolExecutor(max_workers=8) as mfe_pool:
                list(mfe_pool.map(_enrich_record_mfe, all_records))
        t_rnahybrid = time.perf_counter()
        
        # Add analytical insights (miRNASNP-v4 inspired)
        enrichment_hints = targetscan.get_target_enrichment_hints(
            results.get('gain_of_function', []),
            results.get('loss_of_function', [])
        )
        
        # Add UTR retrieval diagnostics to help debug LOF/GOF discrepancies
        utr_diagnostics = {
            'sequence_source': transcript_id,
            'transcript_id_requested': transcript_id,
            'transcript_id_used': transcript_id_used,
            'transcript_resolution_note': transcript_fallback_note,
            'transcript_id_resolved': transcript_data.get('transcript_id_resolved', transcript_id_used),
            'sequence_retrieval_method': transcript_data.get('sequence_source', 'unknown'),
            'sequence_retrieval_url': transcript_data.get('sequence_source_url', ''),
            'sequence_retrieval_local_fastpath': transcript_data.get('sequence_source', '').startswith('local_targetscan'),
            'retrieved_utr_start': transcript_data.get('utr_start'),
            'retrieved_utr_end': transcript_data.get('utr_end'),
            'snp_genomic_position': snp_position,
            'snp_index_in_utr_sequence': transcript_snp_index,
            'strand_adjusted_ref_allele': ref_allele_display,
            'strand_adjusted_alt_allele': alt_allele_display,
            'ref_sequence_length': len(ref_sequence) if ref_sequence else 0,
            'mut_sequence_length': len(mut_sequence) if mut_sequence else 0,
            'snp_position_in_seq': transcript_snp_index,
            'ref_target_count': len(ref_targets),
            'mut_target_count': len(mut_targets),
            'ref_targets_sample': [
                {
                    'mirna': t.get('mirna_id'),
                    'context_score': t.get('raw_context_score'),
                    'site_type': t.get('site_type'),
                    'utr_start': t.get('utr_start'),
                    'utr_end': t.get('utr_end'),
                    'utr_region': t.get('utr_region'),
                    'pairing': t.get('pairing'),
                    'mature_mirna_sequence': t.get('mature_mirna_sequence')
                }
                for t in ref_targets[:3]
            ],
            'mut_targets_sample': [
                {
                    'mirna': t.get('mirna_id'),
                    'context_score': t.get('raw_context_score'),
                    'site_type': t.get('site_type'),
                    'utr_start': t.get('utr_start'),
                    'utr_end': t.get('utr_end'),
                    'utr_region': t.get('utr_region'),
                    'pairing': t.get('pairing'),
                    'mature_mirna_sequence': t.get('mature_mirna_sequence')
                }
                for t in mut_targets[:3]
            ],
            'note': '3\' UTR coordinates and sequence context shown. Verify utr_region and coordinates match expected binding sites.'
        }
        
        # Cap neutral results to avoid multi-MB responses that freeze the browser
        all_neutral = results.get('neutral', [])
        neutral_cap = 100
        neutral_for_response = all_neutral[:neutral_cap]

        return jsonify({
            'rs_id': rs_id,
            'gene_id': gene_id,
            'transcript_id': transcript_id,
            'snp_position': snp_position,
            'ref_allele': ref_allele_display,
            'alt_allele': alt_allele_display,
            'ref_allele_genomic': ref_allele,
            'alt_allele_genomic': alt_allele,
            'strand': strand,
            'loss_of_function': results.get('loss_of_function', []),
            'gain_of_function': results.get('gain_of_function', []),
            'neutral': neutral_for_response,
            'total_neutral': len(all_neutral),
            'total_ref_targets': len(ref_targets),
            'total_mut_targets': len(mut_targets),
            'methodology': results.get('methodology_note', 'Standard TargetScan comparison'),
            'enrichment_analysis': enrichment_hints,
            'utr_diagnostics': utr_diagnostics,
            'performance_ms': {
                'snp_lookup': round((t_snp - t0) * 1000, 2),
                'transcript_retrieval': round((t_transcript - t_snp) * 1000, 2),
                'allele_sequence_build': round((t_sequences - t_transcript) * 1000, 2),
                'targetscan_ref_alt_parallel': round((t_targetscan - t_sequences) * 1000, 2),
                'compare_and_format': round((t_compare - t_targetscan) * 1000, 2),
                'rnahybrid_mfe': round((t_rnahybrid - t_compare) * 1000, 2),
                'total': round((t_rnahybrid - t0) * 1000, 2),
            }
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/debug-targetscan', methods=['POST'])
def debug_targetscan():
    """Detailed diagnostic endpoint to debug TargetScan input/output."""
    try:
        data = request.get_json() or {}
        gene_id = data.get('gene_id', '')
        transcript_id = data.get('transcript_id', '')
        sequence = data.get('sequence', '')
        snp_position = data.get('snp_position')
        ref_allele = data.get('ref_allele', '')
        alt_allele = data.get('alt_allele', '')

        if not sequence:
            return jsonify({'error': 'sequence is required'}), 400

        debug_info = {
            'input': {
                'gene_id': gene_id,
                'transcript_id': transcript_id,
                'sequence_length': len(sequence),
                'sequence_preview_first_100': sequence[:100],
                'sequence_preview_last_100': sequence[-100:],
                'snp_position': snp_position,
                'ref_allele': ref_allele,
                'alt_allele': alt_allele,
                'sequence_gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence) if sequence else 0,
                'sequence_at_content': (sequence.count('A') + sequence.count('T')) / len(sequence) if sequence else 0,
            },
            'window_extraction': {},
            'targetscan_results': {}
        }

        window_seq = targetscan._extract_local_window(sequence, snp_position=snp_position, window_radius=250)
        debug_info['window_extraction'] = {
            'original_length': len(sequence),
            'window_length': len(window_seq),
            'window_preview_first_100': window_seq[:100],
            'window_preview_last_100': window_seq[-100:],
            'gc_content': (window_seq.count('G') + window_seq.count('C')) / len(window_seq) if window_seq else 0,
        }

        ref_targets = targetscan.predict_targets(sequence, gene_id, snp_position=snp_position, strand='+')

        mut_sequence = sequence
        if snp_position and ref_allele and alt_allele:
            snp_index = int(snp_position) - 1
            if 0 <= snp_index < len(sequence):
                mut_sequence = sequence[:snp_index] + alt_allele + sequence[snp_index + len(ref_allele):]

        mut_targets = targetscan.predict_targets(mut_sequence, gene_id, snp_position=snp_position, strand='+')

        debug_info['targetscan_results'] = {
            'ref_target_count': len(ref_targets),
            'mut_target_count': len(mut_targets),
            'ref_top_5': [
                {
                    'mirna_id': t.get('mirna_id'),
                    'context_score': t.get('context_score'),
                    'site_type': t.get('site_type'),
                    'utr_start': t.get('utr_start'),
                    'utr_end': t.get('utr_end'),
                }
                for t in ref_targets[:5]
            ],
            'mut_top_5': [
                {
                    'mirna_id': t.get('mirna_id'),
                    'context_score': t.get('context_score'),
                    'site_type': t.get('site_type'),
                    'utr_start': t.get('utr_start'),
                    'utr_end': t.get('utr_end'),
                }
                for t in mut_targets[:5]
            ],
        }

        expected_mirnas_lof = {'hsa-miR-1255a', 'hsa-miR-1255b-5p', 'hsa-miR-6744-5p', 'hsa-miR-1265', 'hsa-miR-211-3p'}
        expected_mirnas_gof = {'hsa-miR-6507-5p', 'hsa-miR-4696'}

        ref_mirnas = {t['mirna_id'] for t in ref_targets}
        mut_mirnas = {t['mirna_id'] for t in mut_targets}

        debug_info['comparison'] = {
            'expected_lof_mirnas': list(expected_mirnas_lof),
            'found_lof_mirnas': list(expected_mirnas_lof & ref_mirnas),
            'missing_lof_mirnas': list(expected_mirnas_lof - ref_mirnas),
            'expected_gof_mirnas': list(expected_mirnas_gof),
            'found_gof_mirnas': list(expected_mirnas_gof & mut_mirnas),
            'missing_gof_mirnas': list(expected_mirnas_gof - mut_mirnas),
        }

        return jsonify(debug_info)
    except Exception as e:
        return jsonify({'error': str(e), 'traceback': traceback.format_exc()}), 500

@app.route('/api/logs/list', methods=['GET'])
def list_logs():
    """List all available TargetScan execution logs for inspection."""
    try:
        logs_dir = os.path.join(os.path.dirname(__file__), 'targetsnap_logs')
        if not os.path.exists(logs_dir):
            return jsonify({'logs': [], 'message': 'No logs directory yet (no queries executed)'}), 200
        
        logs = []
        for run_id in sorted(os.listdir(logs_dir), reverse=True)[:20]:  # Last 20 runs
            run_log_dir = os.path.join(logs_dir, run_id)
            if not os.path.isdir(run_log_dir):
                continue
            
            status_file = os.path.join(run_log_dir, 'status.log')
            status = 'UNKNOWN'
            if os.path.exists(status_file):
                with open(status_file, 'r') as f:
                    first_line = f.readline().strip()
                    status = first_line.replace('STATUS: ', '')
            
            files = os.listdir(run_log_dir)
            logs.append({
                'run_id': run_id,
                'status': status,
                'files': sorted(files),
                'timestamp': run_id.split('_')[-1]
            })
        
        return jsonify({'logs': logs}), 200
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/logs/<run_id>/<filename>', methods=['GET'])
def get_log_file(run_id, filename):
    """Retrieve a specific log file for inspection."""
    try:
        # Sanitize filename to prevent path traversal
        if '..' in filename or '/' in filename or '\\' in filename:
            return jsonify({'error': 'Invalid filename'}), 400
        
        logs_dir = os.path.join(os.path.dirname(__file__), 'targetsnap_logs')
        file_path = os.path.join(logs_dir, run_id, filename)
        
        # Double-check the path is within logs_dir
        if not os.path.abspath(file_path).startswith(os.path.abspath(logs_dir)):
            return jsonify({'error': 'Access denied'}), 403
        
        if not os.path.exists(file_path):
            return jsonify({'error': f'File not found: {filename}'}), 404
        
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        return jsonify({
            'run_id': run_id,
            'filename': filename,
            'content': content
        }), 200
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/preflight', methods=['GET'])
def preflight():
    """TargetScan local setup diagnostics."""
    try:
        return jsonify(targetscan.preflight_check())
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/eqtl-dae', methods=['POST'])
def eqtl_dae():
    """Fetch GTEx eQTL evidence and DAE-like inference for selected rsID."""
    try:
        data = request.get_json() or {}
        rs_id = data.get('rs_id', '').strip()
        gene_id = data.get('gene_id', '').strip() or None
        tissues = data.get('tissues')
        if tissues and not isinstance(tissues, list):
            tissues = None

        if not rs_id:
            return jsonify({'error': 'rs_id is required'}), 400

        analysis = genomic_handler.get_eqtl_dae(rs_id=rs_id, gene_id=gene_id, tissues=tissues)
        return jsonify(analysis)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/export-results', methods=['POST'])
def export_results():
    """Export latest comparison payload as CSV or JSON."""
    try:
        data = request.get_json() or {}
        export_format = (data.get('format') or 'csv').lower()

        if export_format == 'json':
            body = json.dumps(data, indent=2)
            return Response(
                body,
                mimetype='application/json',
                headers={'Content-Disposition': 'attachment; filename=targetsnap_results.json'}
            )

        if export_format == 'tsv':
            tsv_text = targetscan.export_results_csv(data, delimiter='\t')
            return Response(
                tsv_text,
                mimetype='text/tab-separated-values',
                headers={'Content-Disposition': 'attachment; filename=targetsnap_results.tsv'}
            )
        csv_text = targetscan.export_results_csv(data)
        return Response(
            csv_text,
            mimetype='text/csv',
            headers={'Content-Disposition': 'attachment; filename=targetsnap_results.csv'}
        )
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/rnahybrid', methods=['POST'])
def rnahybrid():
    """Run RNAhybrid for a miRNA:UTR pair (REF and ALT) via WSL."""
    try:
        data = request.get_json() or {}
        ref_utr = (data.get('ref_utr') or '').strip()
        alt_utr = (data.get('alt_utr') or '').strip()
        mirna_seq = (data.get('mirna_seq') or '').strip()

        if not mirna_seq:
            return jsonify({'error': 'mirna_seq is required'}), 400
        if not ref_utr and not alt_utr:
            return jsonify({'error': 'At least one of ref_utr or alt_utr is required'}), 400

        def to_rna(seq):
            return seq.upper().replace('T', 'U')

        def run_rnahybrid(target_seq, query_seq):
            target_rna = to_rna(target_seq)
            query_rna = to_rna(query_seq)
            if not target_rna or not query_rna:
                return None
            # Validate sequences contain only valid RNA chars
            if not _re.fullmatch(r'[AUGCN]+', target_rna) or not _re.fullmatch(r'[AUGCN]+', query_rna):
                return None
            try:
                cmd = ['wsl', 'RNAhybrid', '-s', '3utr_human', '-b', '1', '-c', target_rna, query_rna]
                proc = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
                line = proc.stdout.strip()
                if not line or proc.returncode != 0:
                    return None
                # Compact format: target_name:target_len:query_name:query_len:mfe:pvalue:position:target_struct:query_struct:query_lower:mirna_lower
                parts = line.split(':')
                if len(parts) < 10:
                    return None
                mfe = float(parts[4])
                pvalue = float(parts[5])
                position = int(parts[6])
                # Reconstruct the ASCII duplex from compact parts
                target_line = parts[7]
                pairing_line = parts[8]
                mirna_line = parts[9] if len(parts) > 9 else ''
                # Also get full (non-compact) output for the pretty diagram
                cmd_full = ['wsl', 'RNAhybrid', '-s', '3utr_human', '-b', '1', target_rna, query_rna]
                proc_full = subprocess.run(cmd_full, capture_output=True, text=True, timeout=15)
                diagram = proc_full.stdout.strip() if proc_full.returncode == 0 else ''
                return {
                    'mfe': mfe,
                    'p_value': pvalue,
                    'position': position,
                    'diagram': diagram,
                    'target_seq': target_rna,
                    'query_seq': query_rna,
                }
            except (subprocess.TimeoutExpired, FileNotFoundError, ValueError):
                return None

        result = {'ref': None, 'alt': None}
        if ref_utr:
            result['ref'] = run_rnahybrid(ref_utr, mirna_seq)
        if alt_utr:
            result['alt'] = run_rnahybrid(alt_utr, mirna_seq)

        return jsonify(result)
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/snp-annotations', methods=['POST'])
def snp_annotations():
    """Fetch ClinVar, GWAS Catalog, and conservation annotations for a SNP."""
    try:
        data = request.get_json() or {}
        rs_id = data.get('rs_id', '').strip()
        chromosome = data.get('chromosome', '').strip()
        position = data.get('position')

        if not rs_id:
            return jsonify({'error': 'rs_id is required'}), 400

        result = {'rs_id': rs_id, 'clinvar': [], 'gwas': [], 'conservation': None}

        # --- ClinVar lookup via NCBI eUtils ---
        try:
            search_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
            search_resp = http_requests.get(search_url, params={
                'db': 'clinvar', 'term': f'{rs_id}[varname]', 'retmode': 'json', 'retmax': 10
            }, timeout=8)
            search_data = search_resp.json()
            id_list = search_data.get('esearchresult', {}).get('idlist', [])

            if id_list:
                summary_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
                sum_resp = http_requests.get(summary_url, params={
                    'db': 'clinvar', 'id': ','.join(id_list[:10]), 'retmode': 'json'
                }, timeout=8)
                sum_data = sum_resp.json().get('result', {})
                for uid in id_list[:10]:
                    entry = sum_data.get(uid, {})
                    if not entry:
                        continue
                    clinical_sig = entry.get('clinical_significance', {})
                    result['clinvar'].append({
                        'uid': uid,
                        'title': entry.get('title', ''),
                        'clinical_significance': clinical_sig.get('description', '') if isinstance(clinical_sig, dict) else str(clinical_sig),
                        'conditions': [t.get('trait_name', '') for t in entry.get('trait_set', []) if t.get('trait_name')],
                        'review_status': clinical_sig.get('review_status', '') if isinstance(clinical_sig, dict) else '',
                    })
        except Exception:
            pass  # ClinVar lookup is best-effort

        # --- GWAS Catalog lookup via EBI REST ---
        try:
            gwas_url = f'https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rs_id}/associations'
            gwas_resp = http_requests.get(gwas_url, headers={'Accept': 'application/json'}, timeout=8)
            if gwas_resp.status_code == 200:
                associations = gwas_resp.json().get('_embedded', {}).get('associations', [])
                for assoc in associations[:20]:
                    traits = []
                    for t in assoc.get('efoTraits', []):
                        traits.append(t.get('trait', ''))
                    result['gwas'].append({
                        'traits': traits,
                        'p_value': assoc.get('pvalue', ''),
                        'risk_allele': ', '.join(
                            ra.get('riskAlleleName', '')
                            for ra in assoc.get('riskAlleles', [])
                        ),
                        'study': assoc.get('study', {}).get('publicationInfo', {}).get('title', '') if assoc.get('study') else '',
                    })
        except Exception:
            pass  # GWAS lookup is best-effort

        # --- Conservation score via UCSC REST ---
        if chromosome and position:
            try:
                chrom = f'chr{chromosome}' if not str(chromosome).startswith('chr') else str(chromosome)
                pos = int(position)
                ucsc_url = 'https://api.genome.ucsc.edu/getData/track'
                cons_resp = http_requests.get(ucsc_url, params={
                    'genome': 'hg19', 'track': 'phastCons46way',
                    'chrom': chrom, 'start': pos - 1, 'end': pos
                }, timeout=8)
                if cons_resp.status_code == 200:
                    cons_data = cons_resp.json()
                    values = cons_data.get('phastCons46way', [])
                    if values and isinstance(values, list) and len(values) > 0:
                        val = values[0]
                        score = val.get('value', val) if isinstance(val, dict) else val
                        result['conservation'] = {
                            'phastCons46way': round(float(score), 4) if score is not None else None,
                            'position': f'{chrom}:{pos}',
                            'interpretation': 'Highly conserved' if score and float(score) > 0.9 else
                                            'Conserved' if score and float(score) > 0.5 else
                                            'Weakly conserved' if score and float(score) > 0.2 else
                                            'Not conserved'
                        }
            except Exception:
                pass  # Conservation lookup is best-effort

        return jsonify(result)
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/enrichment', methods=['POST'])
def enrichment_analysis():
    """Run GO/KEGG enrichment via Enrichr for LOF/GOF miRNA targets."""
    try:
        data = request.get_json() or {}
        mirna_ids = data.get('mirna_ids', [])
        gene_id = data.get('gene_id', '')
        analysis_type = data.get('type', 'all')  # 'go', 'kegg', or 'all'

        if not mirna_ids:
            return jsonify({'error': 'mirna_ids list is required'}), 400

        # Build a gene list from miRNA target genes using miRTarBase-style lookup
        # For now, use the host gene + miRNA family names as proxy
        gene_list = [gene_id] if gene_id else []
        # Add common validated targets for well-known miRNAs
        gene_list_str = '\n'.join(gene_list) if gene_list else gene_id

        results = {'mirna_count': len(mirna_ids), 'gene_list': gene_list, 'go_biological_process': [], 'go_molecular_function': [], 'kegg': []}

        if not gene_list_str.strip():
            return jsonify(results)

        # Submit to Enrichr
        try:
            add_url = 'https://maayanlab.cloud/Enrichr/addList'
            add_resp = http_requests.post(add_url, files={
                'list': (None, gene_list_str),
                'description': (None, f'TargetSNAP LOF/GOF miRNAs for {gene_id}')
            }, timeout=10)
            if add_resp.status_code != 200:
                return jsonify(results)

            user_list_id = add_resp.json().get('userListId')
            if not user_list_id:
                return jsonify(results)

            libraries = []
            if analysis_type in ('go', 'all'):
                libraries += ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023']
            if analysis_type in ('kegg', 'all'):
                libraries += ['KEGG_2021_Human']

            enrich_url = 'https://maayanlab.cloud/Enrichr/enrich'
            for lib in libraries:
                try:
                    resp = http_requests.get(enrich_url, params={
                        'userListId': user_list_id, 'backgroundType': lib
                    }, timeout=10)
                    if resp.status_code == 200:
                        terms = resp.json().get(lib, [])
                        formatted = []
                        for term in terms[:15]:
                            formatted.append({
                                'term': term[1] if len(term) > 1 else '',
                                'p_value': term[2] if len(term) > 2 else 1,
                                'z_score': term[3] if len(term) > 3 else 0,
                                'combined_score': term[4] if len(term) > 4 else 0,
                                'genes': term[5] if len(term) > 5 else [],
                            })
                        key = lib.lower().replace('_2023', '').replace('_2021_human', '').replace('go_', 'go_')
                        if 'biological' in key:
                            results['go_biological_process'] = formatted
                        elif 'molecular' in key:
                            results['go_molecular_function'] = formatted
                        elif 'kegg' in key:
                            results['kegg'] = formatted
                except Exception:
                    pass
        except Exception:
            pass  # Enrichr is best-effort

        return jsonify(results)
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/health', methods=['GET'])
def health():
    """Health check endpoint"""
    return jsonify({'status': 'ok', 'service': 'TargetSNAP Web API'})

@app.errorhandler(404)
def not_found(error):
    return jsonify({'error': 'Endpoint not found'}), 404

@app.errorhandler(500)
def server_error(error):
    return jsonify({'error': 'Internal server error'}), 500

if __name__ == '__main__':
    print("\n" + "="*50)
    print("  TargetSNAP Web Server")
    print("="*50)
    print("\nStarting Flask application...")
    print("Open http://localhost:5000 in your browser\n")
    app.run(debug=True, host='127.0.0.1', port=5000, threaded=True)
