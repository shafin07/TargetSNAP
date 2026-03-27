"""TargetSNAP – Flask web application.

Provides a web interface to:
1. Enter a dbSNP rs ID.
2. View mapped genes and their transcript isoforms.
3. Select a gene/transcript to run TargetScan on both ref and alt UTRs.
4. View Gain-of-Function / Loss-of-Function miRNA results.
"""

from flask import Flask, jsonify, render_template, request

from targetscan.analysis import classify
from targetscan.snp_lookup import get_mapped_genes, get_snp_info
from targetscan.strand_handler import apply_snp, effective_alleles, snp_offset_in_utr
from targetscan.targetscan_runner import run_targetscan
from targetscan.utr_retrieval import get_utr_sequence

app = Flask(__name__)


# ---------------------------------------------------------------------------
# Pages
# ---------------------------------------------------------------------------

@app.route("/")
def index():
    """Render the main search page."""
    return render_template("index.html")


# ---------------------------------------------------------------------------
# API endpoints
# ---------------------------------------------------------------------------

@app.route("/api/snp/<rs_id>")
def api_snp(rs_id):
    """Return SNP metadata for *rs_id*."""
    info = get_snp_info(rs_id)
    if info is None:
        return jsonify({"error": f"SNP {rs_id} not found"}), 404
    return jsonify(info)


@app.route("/api/genes/<rs_id>")
def api_genes(rs_id):
    """Return genes with 3'UTR overlapping *rs_id*."""
    genes = get_mapped_genes(rs_id)
    if not genes:
        return jsonify({"error": "No genes with 3'UTR found for this SNP"}), 404
    return jsonify(genes)


@app.route("/api/analyze", methods=["POST"])
def api_analyze():
    """Run TargetScan on ref and alt UTR, return GOF/LOF results.

    Expected JSON body::

        {
            "rs_id": "rs11552978",
            "transcript_id": "ENST00000...",
            "gene_symbol": "DCBLD2",
            "strand": "-",
            "utr_start": 98514785,
            "utr_end": 98518215,
            "snp_position": 98517568,
            "ref_allele": "G",
            "alt_allele": "A",
            "chr": "3"
        }
    """
    data = request.get_json(force=True)

    gene_symbol = data["gene_symbol"]
    transcript_id = data["transcript_id"]
    strand = data["strand"]
    utr_start = int(data["utr_start"])
    utr_end = int(data["utr_end"])
    snp_pos = int(data["snp_position"])
    ref_allele = data["ref_allele"]
    alt_allele = data["alt_allele"]
    chrom = data["chr"]

    # 1. Fetch the reference UTR sequence from UCSC hg19
    try:
        ref_seq = get_utr_sequence(chrom, utr_start, utr_end)
    except Exception as exc:
        return jsonify({"error": f"Failed to fetch UTR sequence: {exc}"}), 500

    # 2. Build the alternate UTR sequence
    offset = snp_offset_in_utr(snp_pos, utr_start, utr_end, strand)
    eff_ref, eff_alt = effective_alleles(ref_allele, alt_allele, strand)

    # Validate that the reference base matches what we fetched
    actual_base = ref_seq[offset] if 0 <= offset < len(ref_seq) else "?"
    if actual_base != eff_ref.upper():
        note = (
            f"Warning: expected ref base '{eff_ref}' at offset {offset} "
            f"but found '{actual_base}' in the fetched sequence."
        )
    else:
        note = None

    alt_seq = apply_snp(ref_seq, offset, alt_allele, strand)

    # 3. Run TargetScan on both sequences
    try:
        ref_results = run_targetscan(gene_symbol, transcript_id, ref_seq)
    except Exception as exc:
        return jsonify({"error": f"TargetScan ref run failed: {exc}"}), 500

    try:
        alt_results = run_targetscan(gene_symbol, transcript_id, alt_seq)
    except Exception as exc:
        return jsonify({"error": f"TargetScan alt run failed: {exc}"}), 500

    # 4. Classify miRNAs as GOF / LOF
    classification = classify(ref_results, alt_results)

    response = {
        "gene_symbol": gene_symbol,
        "transcript_id": transcript_id,
        "strand": strand,
        "snp_position": snp_pos,
        "ref_allele": ref_allele,
        "alt_allele": alt_allele,
        "eff_ref_allele": eff_ref,
        "eff_alt_allele": eff_alt,
        "utr_length": len(ref_seq),
        "snp_offset": offset,
        "gof": classification["gof"],
        "lof": classification["lof"],
        "unchanged": classification["unchanged"],
    }
    if note:
        response["note"] = note

    return jsonify(response)


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    app.run(debug=True, port=5000)
