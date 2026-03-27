"""Retrieve 3'UTR sequences from the UCSC DAS server (hg19)."""

import xml.etree.ElementTree as ET

import requests

from config import UCSC_DAS_URL

_REQUEST_TIMEOUT = 30


def get_utr_sequence(chrom, start, end):
    """Fetch the genomic DNA sequence from UCSC hg19 DAS.

    Parameters
    ----------
    chrom : str
        Chromosome name, e.g. ``"1"`` or ``"chr1"``.
    start : int
        1-based start coordinate.
    end : int
        1-based end coordinate.

    Returns
    -------
    str
        Upper-cased DNA sequence.
    """
    if not str(chrom).startswith("chr"):
        chrom = f"chr{chrom}"

    params = {"segment": f"{chrom}:{start},{end}"}
    resp = requests.get(UCSC_DAS_URL, params=params, timeout=_REQUEST_TIMEOUT)
    resp.raise_for_status()

    root = ET.fromstring(resp.text)
    dna_elem = root.find(".//DNA")
    if dna_elem is None or dna_elem.text is None:
        raise ValueError(
            f"No sequence returned for {chrom}:{start}-{end}"
        )
    return dna_elem.text.strip().replace("\n", "").upper()
