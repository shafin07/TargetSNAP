"""Look up SNP information from Ensembl GRCh37 REST API.

Given an rs ID, retrieve the SNP's chromosomal position, ref/alt alleles,
and the genes whose 3'UTR regions overlap the SNP.
"""

import requests

from config import ENSEMBL_REST_URL

# Timeout in seconds for all outbound HTTP requests
_REQUEST_TIMEOUT = 30


def get_snp_info(rs_id):
    """Return SNP metadata from Ensembl GRCh37.

    Parameters
    ----------
    rs_id : str
        dbSNP rs identifier, e.g. ``"rs11552978"``.

    Returns
    -------
    dict or None
        ``{"id", "chr", "position", "alleles", "strand",
          "mappings": [{"chr", "start", "end", "allele_string", "strand"}]}``
    """
    url = f"{ENSEMBL_REST_URL}/variation/human/{rs_id}"
    params = {"content-type": "application/json"}
    resp = requests.get(url, params=params, timeout=_REQUEST_TIMEOUT)
    if resp.status_code != 200:
        return None
    data = resp.json()

    mappings = data.get("mappings", [])
    if not mappings:
        return None

    primary = mappings[0]
    return {
        "id": data.get("name", rs_id),
        "chr": primary.get("seq_region_name"),
        "position": primary.get("start"),
        "alleles": primary.get("allele_string", "").split("/"),
        "strand": primary.get("strand"),
        "mappings": mappings,
    }


def get_mapped_genes(rs_id):
    """Return genes whose 3'UTR overlaps *rs_id*.

    Uses the Ensembl ``/overlap/region`` endpoint filtered to
    ``feature=transcript`` so every transcript isoform is included.

    Returns
    -------
    list[dict]
        Each dict: ``{"gene_symbol", "gene_id", "transcript_id",
        "strand", "utr_start", "utr_end", "biotype"}``.
    """
    snp = get_snp_info(rs_id)
    if snp is None:
        return []

    chrom = snp["chr"]
    pos = snp["position"]

    url = f"{ENSEMBL_REST_URL}/overlap/region/human/{chrom}:{pos}-{pos}"
    params = {
        "content-type": "application/json",
        "feature": "transcript",
    }
    resp = requests.get(url, params=params, timeout=_REQUEST_TIMEOUT)
    if resp.status_code != 200:
        return []

    transcripts = resp.json()
    results = []
    seen = set()
    for tx in transcripts:
        tx_id = tx.get("id", "")
        if tx_id in seen:
            continue
        seen.add(tx_id)

        utr_info = _get_3utr_for_transcript(tx_id)
        if utr_info is None:
            continue

        strand_int = tx.get("strand", 1)
        results.append(
            {
                "gene_symbol": tx.get("external_name", ""),
                "gene_id": tx.get("Parent", ""),
                "transcript_id": tx_id,
                "strand": "+" if strand_int == 1 else "-",
                "utr_start": utr_info["start"],
                "utr_end": utr_info["end"],
                "biotype": tx.get("biotype", ""),
            }
        )
    return results


def _get_3utr_for_transcript(transcript_id):
    """Return the 3'UTR coordinates for *transcript_id*, or ``None``."""
    url = f"{ENSEMBL_REST_URL}/lookup/id/{transcript_id}"
    params = {"content-type": "application/json", "expand": "1"}
    resp = requests.get(url, params=params, timeout=_REQUEST_TIMEOUT)
    if resp.status_code != 200:
        return None

    data = resp.json()
    utrs = [
        ex
        for ex in data.get("UTR", data.get("Exon", []))
        if isinstance(ex, dict)
    ]

    three_prime = data.get("three_prime_UTR")
    if three_prime:
        return {"start": three_prime["start"], "end": three_prime["end"]}

    # Fallback: try fetching the sequence feature for the 3'UTR
    url2 = f"{ENSEMBL_REST_URL}/overlap/id/{transcript_id}"
    params2 = {"content-type": "application/json", "feature": "three_prime_UTR"}
    resp2 = requests.get(url2, params=params2, timeout=_REQUEST_TIMEOUT)
    if resp2.status_code == 200:
        utr_list = resp2.json()
        if utr_list:
            starts = [u["start"] for u in utr_list]
            ends = [u["end"] for u in utr_list]
            return {"start": min(starts), "end": max(ends)}

    if utrs:
        return {"start": utrs[-1].get("start"), "end": utrs[-1].get("end")}

    return None
