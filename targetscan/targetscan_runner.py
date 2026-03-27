"""Execute TargetScan 7.0 Perl scripts and parse results.

TargetScan expects tab-separated input files and produces tab-separated
output.  This module creates temporary input files, invokes the Perl
scripts, and parses the results into Python dicts.
"""

import os
import subprocess
import tempfile
import uuid

from config import (
    MIRNA_FAMILY_FILE,
    TARGETSCAN_SCRIPT,
    TEMP_DIR,
)


def _write_utr_file(gene_symbol, transcript_id, species_id, sequence, path):
    """Write a TargetScan UTR input file.

    Format: ``<transcript_id>\\t<species_id>\\t<sequence>``
    """
    with open(path, "w") as fh:
        fh.write(f"{transcript_id}\t{species_id}\t{sequence}\n")


def _write_mirna_file(mirna_family_file, path):
    """Copy or symlink the miRNA family file into the working dir.

    The miRNA family file is expected to be already in TargetScan
    format: ``<miR_family>\\t<seed_region>\\t<species_id>``.
    """
    # If the source and destination are different, copy the file.
    if os.path.abspath(mirna_family_file) != os.path.abspath(path):
        with open(mirna_family_file, "r") as src, open(path, "w") as dst:
            dst.write(src.read())


def run_targetscan(gene_symbol, transcript_id, utr_sequence, species_id="9606"):
    """Run the TargetScan 7.0 prediction script.

    Parameters
    ----------
    gene_symbol : str
        Gene symbol, used only for labelling.
    transcript_id : str
        Ensembl transcript identifier.
    utr_sequence : str
        3'UTR nucleotide sequence (A/C/G/T).
    species_id : str
        NCBI taxonomy ID.  ``"9606"`` for human.

    Returns
    -------
    list[dict]
        Each dict represents one predicted miRNA–target site with keys
        such as ``miRNA_family``, ``seed_match_type``, ``UTR_start``,
        ``UTR_end``, and ``context_score``.
    """
    run_id = uuid.uuid4().hex[:8]
    utr_path = os.path.join(TEMP_DIR, f"utr_{run_id}.txt")
    mirna_path = os.path.join(TEMP_DIR, f"mirna_{run_id}.txt")
    out_path = os.path.join(TEMP_DIR, f"out_{run_id}.txt")

    try:
        _write_utr_file(gene_symbol, transcript_id, species_id, utr_sequence, utr_path)
        _write_mirna_file(MIRNA_FAMILY_FILE, mirna_path)

        cmd = [
            "perl",
            TARGETSCAN_SCRIPT,
            mirna_path,
            utr_path,
            out_path,
        ]
        subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            timeout=120,
        )

        return _parse_output(out_path)
    finally:
        for p in (utr_path, mirna_path, out_path):
            if os.path.exists(p):
                os.remove(p)


def _parse_output(output_path):
    """Parse TargetScan tab-separated output into a list of dicts."""
    results = []
    if not os.path.exists(output_path):
        return results

    with open(output_path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if header is None:
                header = cols
                continue
            row = dict(zip(header, cols))
            results.append(row)
    return results
