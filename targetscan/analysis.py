"""
analysis.py – TargetSNAP core analysis.

Workflow
--------
1. Look up the rs ID via the Ensembl Variation REST API to get genomic
   coordinates, reference and alternate alleles, and the most-severe
   functional consequence.
2. Retrieve a genomic sequence window (±WINDOW nt) around the SNP from
   Ensembl Sequence REST API (+ strand only).
3. For every gene overlapping the SNP, derive the mRNA-oriented window and
   apply the alt allele to create the alternate sequence.
4. Search both sequences for all miRNA seed-match sites (8mer, 7mer-m8,
   7mer-A1, 6mer) from a curated set of 50 human miRNA families.
5. Report:
     - LOF  – sites present in the reference but absent in the alternate
               (SNP disrupts an existing miRNA target).
     - GOF  – sites absent in the reference but present in the alternate
               (SNP creates a new miRNA target).
"""

from __future__ import annotations

import requests

from .strand_handler import complement, reverse_complement

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ENSEMBL_REST = "https://rest.ensembl.org"
WINDOW = 100               # nucleotides on each side of the SNP
FILTER_DISTANCE = 50       # max nt between SNP and a seed-match change to report
REQUEST_TIMEOUT = 15       # seconds

# ---------------------------------------------------------------------------
# Curated set of conserved human miRNA mature sequences (miRBase v22).
# Seed positions 2-8 (1-indexed) are used to compute UTR seed-match motifs.
# ---------------------------------------------------------------------------

MIRNA_SEQUENCES: dict[str, str] = {
    # let-7 / miR-98 family  (seed 2-8: GAGGUAG → UTR 7mer-m8: CTACCTC)
    "hsa-let-7a-5p":   "UGAGGUAGUAGGUUGUAUAGUU",
    "hsa-let-7b-5p":   "UGAGGUAGUGAGGUUGUAUAGUU",
    "hsa-let-7c-5p":   "UGAGGUAGUAGGUUGUAUGGUU",
    "hsa-let-7d-5p":   "AGAGGUAGUGAGGUUGUAUAGUU",
    "hsa-let-7e-5p":   "UGAGGUAGGAGGUUGUAUAGUU",
    "hsa-let-7f-5p":   "UGAGGUAGUGAUGUUGUAUAGUU",
    "hsa-let-7g-5p":   "UGAGGUAGUAGUUUGUACAGUU",
    "hsa-let-7i-5p":   "UGAGGUAGUAGUUUGUGCUGUU",
    "hsa-miR-98-5p":   "UGAGGUAGUAAGUUGUAUUGUU",
    # miR-17/20/93/106 family (seed: AAAGUGC → GCACTTT)
    "hsa-miR-17-5p":   "CAAAGUGCUUACAGUGCAGGUAG",
    "hsa-miR-20a-5p":  "UAAAGUGCUUAUAGUGCAGGUAG",
    "hsa-miR-20b-5p":  "CAAAGUGCUCAUAGUGCAGGUAG",
    "hsa-miR-93-5p":   "CAAAGUGCUGUUCGUGCAGGUAG",
    "hsa-miR-106a-5p": "AAAAGUGCUUACAGUGCAGGUAG",
    "hsa-miR-106b-5p": "UAAAGUGCUGACAGUGCAGAU",
    # miR-19 family  (seed: UGCAAAU → ATTTGCA)
    "hsa-miR-19a-3p":  "UGUGCAAAUCUAUGCAAAACUGA",
    "hsa-miR-19b-3p":  "UGUGCAAAUCCAUGCAAAACUGA",
    # miR-21  (seed: AGCUUAU → ATAAGCT)
    "hsa-miR-21-5p":   "UAGCUUAUCAGACUGAUGUUGA",
    # miR-23 family  (seed: UCACAUU → AATGTGA)
    "hsa-miR-23a-3p":  "AUCACAUUGCCAGGGAUUUCC",
    "hsa-miR-23b-3p":  "AUCACAUUGCCAGGGAUUACC",
    # miR-24  (seed: GGCUCAG → CTGAGCC)
    "hsa-miR-24-3p":   "UGGCUCAGUUCAGCAGGAACAG",
    # miR-27 family  (seed: UCACAGU → ACTGTGA)
    "hsa-miR-27a-3p":  "UUCACAGUGGCUAAGUUCCGC",
    "hsa-miR-27b-3p":  "UUCACAGUGGCUAAGUUCUGC",
    # miR-29 family  (seed: AGCACCA → TGGTGCT)
    "hsa-miR-29a-3p":  "UAGCACCAUCUGAAAUCGGUUA",
    "hsa-miR-29b-3p":  "UAGCACCAUUUGAAAUCAGUGUU",
    "hsa-miR-29c-3p":  "UAGCACCAUUUGAAAUCGGUUA",
    # miR-30 family  (seed: GUAAACA → TGTTTAC)
    "hsa-miR-30a-5p":  "UGUAAACAUCCUCGACUGGAAG",
    "hsa-miR-30b-5p":  "UGUAAACAUCCUACACUCAGC",
    "hsa-miR-30c-5p":  "UGUAAACAUCCUACACUCUCAGC",
    "hsa-miR-30d-5p":  "UGUAAACAUCCCCGACUGGAAG",
    "hsa-miR-30e-5p":  "UGUAAACAUCCUUGACUGGAAG",
    # miR-34 family  (seed: GGCAGUG → CACTGCC)
    "hsa-miR-34a-5p":  "UGGCAGUGUGGUUAGCUGGUUGU",
    "hsa-miR-34b-5p":  "AGGCAGUGUCAUUAGCUGAUUGU",
    "hsa-miR-34c-5p":  "AGGCAGUGUAGUUAGCUGAUUGC",
    # miR-92 family  (seed: AUUGCAC → GTGCAAT)
    "hsa-miR-92a-3p":  "UAUUGCACUUGUCCCGGCCUGU",
    "hsa-miR-92b-3p":  "UAUUGCACUCGUCCCGGCCUCC",
    # miR-96/182/183 family
    "hsa-miR-96-5p":   "UUUGGCACUAGCACAUUUUUGCU",
    "hsa-miR-182-5p":  "UUUGGCAAUGGUAGAACUCACAC",
    "hsa-miR-183-5p":  "UAUGGCACUGGUAGAAUUCACU",
    # miR-1/206 family  (seed: GGAAUGU → ACATTCC)
    "hsa-miR-1-3p":    "UGGAAUGUAAAGAAGUAUGUAU",
    "hsa-miR-206":     "UGGAAUGUAAGGAAGUGUGUGG",
    # miR-7  (seed: GGAAGAC → GTCTTCC)
    "hsa-miR-7-5p":    "UGGAAGACUAGUGAUUUUGUUGU",
    # miR-9  (seed: CUUUGGU → ACCAAAG)
    "hsa-miR-9-5p":    "UCUUUGGUUAUCUAGCUGUAUGA",
    # miR-10 family  (seed: ACCCUGU → ACAGGGT)
    "hsa-miR-10a-5p":  "UACCCUGUAGAUCCGAAUUUGUG",
    "hsa-miR-10b-5p":  "UACCCUGUAGAACCGAAUUUGUG",
    # miR-122  (seed: GGAGUGU → ACACTCC)
    "hsa-miR-122-5p":  "UGGAGUGUGACAAUGGUGUUUG",
    # miR-125 family  (seed: UCCCUGA → TCAGGGA)
    "hsa-miR-125a-5p": "UCCCUGAGACCCUUUAACCUGUG",
    "hsa-miR-125b-5p": "UCCCUGAGACCCUAACUUGUGA",
    # miR-126  (seed: CGUACCG → CGGTACG)
    "hsa-miR-126-3p":  "UCGUACCGUGAGUAAUAAUGCG",
    # miR-133 family  (seed: UUGGUCCC → GGGACCAA)
    "hsa-miR-133a-3p": "UUUGGUCCCCUUCAACCAGCUG",
    "hsa-miR-133b":    "UUUGGUCCCCUUCAACCAGCUA",
    # miR-143  (seed: GAGAUGA → TCATCTC)
    "hsa-miR-143-3p":  "UGAGAUGAAGCACUGUAGCUCA",
    # miR-145  (seed: UCCAGUU → AACTGGA)
    "hsa-miR-145-5p":  "GUCCAGUUUUCCCAGGAAUCCCU",
    # miR-146 family  (seed: UGAGAAC → GTTCTCA)
    "hsa-miR-146a-5p": "UGAGAACUGAAUUCCAUGGGUU",
    "hsa-miR-146b-5p": "UGAGAACUGAAUUCCAUAGGCU",
    # miR-150  (seed: UCUCCCAA → TTGGGAGA)
    "hsa-miR-150-5p":  "UCUCCCAACCCUUGUACCAGUG",
    # miR-155  (seed: UAAUGCU → AGCATTA)
    "hsa-miR-155-5p":  "UUAAUGCUAAUUGUGAUAGGGGU",
    # miR-181 family  (seed: ACAUUCA → TGAATGT)
    "hsa-miR-181a-5p": "AACAUUCAACGCUGUCGGUGAGU",
    "hsa-miR-181b-5p": "AACAUUCAUUGCUGUCGGUGGG",
    "hsa-miR-181c-5p": "AACAUUCAACCUGUCGGUGAGU",
    # miR-200/141/429 family  (seed: AAUACUG → CAGTATT)
    "hsa-miR-200a-3p": "UAACACUGUCUGGUAAAGAUGG",
    "hsa-miR-200b-3p": "UAAUACUGCCUGGUAAUGAUGA",
    "hsa-miR-200c-3p": "UAAUACUGCCGGGUAAUGAUGGA",
    "hsa-miR-141-3p":  "UAACACUGUCUGGUAAAGAUGG",
    "hsa-miR-429":     "UAAUACUGUCUGGUAAAACCGU",
    # miR-210  (seed: UGUGCGU → ACGCACA)
    "hsa-miR-210-3p":  "CUGUGCGUGUGACAGCGGCUGA",
    # miR-221/222 family  (seed: GCUACAU → ATGTAGC)
    "hsa-miR-221-3p":  "AGCUACAUUGUCUGCUGGGUUUC",
    "hsa-miR-222-3p":  "AGCUACAUCUGGCUACUGGGU",
    # miR-223  (seed: GUCAGUU → AACTGAC)
    "hsa-miR-223-3p":  "UGUCAGUUUGUCAAAUACCCCA",
    # miR-375  (seed: UUGUUCG → CGAACAA)
    "hsa-miR-375":     "UUUGUUCGUUCGGCUCGCGUGA",
}

# ---------------------------------------------------------------------------
# Pre-compute seed-match motifs (as they appear in the 3′ UTR, DNA)
# ---------------------------------------------------------------------------

def _rc_window(mirna_seq: str, start: int, end: int) -> str:
    """Reverse complement of miRNA positions [start:end] (0-indexed, DNA)."""
    dna = mirna_seq.upper().replace("U", "T")
    return reverse_complement(dna[start:end])


def _build_seed_table() -> dict[str, dict[str, str]]:
    table: dict[str, dict[str, str]] = {}
    for name, seq in MIRNA_SEQUENCES.items():
        # positions 2-8 (1-indexed) → Python slice [1:8]
        table[name] = {
            "8mer":     "A" + _rc_window(seq, 1, 8),
            "7mer-m8":  _rc_window(seq, 1, 8),
            "7mer-A1":  "A" + _rc_window(seq, 1, 7),
            "6mer":     _rc_window(seq, 1, 7),
        }
    return table


MIRNA_SEEDS: dict[str, dict[str, str]] = _build_seed_table()

# ---------------------------------------------------------------------------
# Ensembl REST helpers
# ---------------------------------------------------------------------------

_HEADERS = {"Content-Type": "application/json"}


def _get(url: str, params: dict | None = None):
    r = requests.get(url, headers=_HEADERS, params=params, timeout=REQUEST_TIMEOUT)
    r.raise_for_status()
    return r.json()


def lookup_variant(rsid: str) -> dict:
    """Fetch variant metadata from Ensembl (coordinates, alleles, consequence)."""
    return _get(f"{ENSEMBL_REST}/variation/human/{rsid}")


def get_sequence(chrom: str, start: int, end: int) -> str:
    """Return the + strand genomic sequence for the given 1-based region."""
    url = f"{ENSEMBL_REST}/sequence/region/human/{chrom}:{start}..{end}:1"
    return _get(url)["seq"]


def get_overlapping_genes(chrom: str, pos: int) -> list[dict]:
    """Return gene features overlapping *pos* (Ensembl overlap endpoint)."""
    url = f"{ENSEMBL_REST}/overlap/region/human/{chrom}:{pos}-{pos}"
    try:
        data = _get(url, params={"feature": "gene", "content-type": "application/json"})
        return [g for g in data if isinstance(g, dict)]
    except Exception:
        return []

# ---------------------------------------------------------------------------
# miRNA seed-match searching
# ---------------------------------------------------------------------------

SiteMatch = dict  # keys: mirna, site_type, seed, start, end


def find_seed_matches(sequence: str) -> list[SiteMatch]:
    """Return every miRNA seed-match site found in *sequence* (DNA, 5′→3′)."""
    seq = sequence.upper()
    matches: list[SiteMatch] = []
    seen: set[tuple] = set()
    for mirna, sites in MIRNA_SEEDS.items():
        for site_type, seed in sites.items():
            if not seed:
                continue
            pos = 0
            while True:
                idx = seq.find(seed, pos)
                if idx == -1:
                    break
                key = (mirna, site_type, idx)
                if key not in seen:
                    seen.add(key)
                    matches.append(
                        {
                            "mirna": mirna,
                            "site_type": site_type,
                            "seed": seed,
                            "start": idx,
                            "end": idx + len(seed),
                        }
                    )
                pos = idx + 1
    return matches


def _diff_matches(
    ref_matches: list[SiteMatch], alt_matches: list[SiteMatch]
) -> tuple[list[SiteMatch], list[SiteMatch], list[SiteMatch]]:
    """Return (lof, gof, neutral) site lists."""

    def key(m: SiteMatch) -> tuple:
        return (m["mirna"], m["site_type"], m["start"])

    ref_map = {key(m): m for m in ref_matches}
    alt_map = {key(m): m for m in alt_matches}

    lof = [ref_map[k] for k in ref_map if k not in alt_map]
    gof = [alt_map[k] for k in alt_map if k not in ref_map]
    neutral = [ref_map[k] for k in ref_map if k in alt_map]
    return lof, gof, neutral

# ---------------------------------------------------------------------------
# Main analysis entry point
# ---------------------------------------------------------------------------

AnalysisResult = dict


def analyze_snp(rsid: str) -> AnalysisResult:
    """Look up *rsid*, retrieve sequence context, and classify miRNA site changes.

    Returns a dict with keys:
        rsid, chromosome, position, ref, alt, consequence,
        genes, analyses, error.

    Each element in *analyses* covers one gene strand orientation and
    contains: strand, ref_seq, alt_seq, snp_idx, ref_allele, alt_allele,
    lof, gof, neutral_count, seq_mismatch.
    """
    rsid = rsid.strip()
    if not rsid.lower().startswith("rs"):
        rsid = "rs" + rsid

    result: AnalysisResult = {
        "rsid": rsid,
        "error": None,
        "chromosome": None,
        "position": None,
        "ref": None,
        "alt": [],
        "consequence": None,
        "genes": [],
        "analyses": [],
    }

    # ── 1. Variant metadata ──────────────────────────────────────────────────
    try:
        variant = lookup_variant(rsid)
    except requests.exceptions.HTTPError as exc:
        if exc.response is not None and exc.response.status_code == 404:
            result["error"] = f"rs ID '{rsid}' was not found in Ensembl. Please check the identifier."
        else:
            result["error"] = f"Ensembl Variation API error: {exc}"
        return result
    except requests.exceptions.RequestException as exc:
        result["error"] = f"Network error while retrieving variant data: {exc}"
        return result

    mappings = variant.get("mappings", [])
    if not mappings:
        result["error"] = f"No genomic mapping found for {rsid}."
        return result

    # Prefer GRCh38 assembly
    mapping = next(
        (m for m in mappings if m.get("assembly_name") == "GRCh38"), mappings[0]
    )
    chrom: str = mapping["seq_region_name"]
    pos: int = int(mapping["start"])
    allele_str: str = mapping.get("allele_string", "N/N")
    alleles = allele_str.split("/")
    ref = alleles[0] if alleles else "N"
    alts = alleles[1:] if len(alleles) > 1 else []

    result.update(
        {
            "chromosome": chrom,
            "position": pos,
            "ref": ref,
            "alt": alts,
            "consequence": variant.get("most_severe_consequence", "unknown"),
        }
    )

    if not alts:
        result["error"] = "No alternate allele reported for this variant."
        return result
    if len(ref) != 1 or any(len(a) != 1 for a in alts):
        result["error"] = (
            "TargetSNAP currently supports only single-nucleotide substitution (SNV) variants."
        )
        return result

    # ── 2. Overlapping genes ─────────────────────────────────────────────────
    try:
        raw_genes = get_overlapping_genes(chrom, pos)
        result["genes"] = [
            {
                "id": g.get("id", "?"),
                "name": g.get("external_name") or g.get("id", "?"),
                "strand": "+" if g.get("strand", 1) == 1 else "-",
                "biotype": g.get("biotype", "?"),
            }
            for g in raw_genes
        ]
    except Exception:
        raw_genes = []

    # ── 3. + strand genomic sequence window ─────────────────────────────────
    win_start = max(1, pos - WINDOW)
    win_end = pos + WINDOW

    try:
        ref_seq_plus: str = get_sequence(chrom, win_start, win_end)
    except requests.exceptions.RequestException as exc:
        result["error"] = f"Could not retrieve genomic sequence: {exc}"
        return result

    # 0-based index of SNP in the + strand window
    snp_idx_plus = pos - win_start

    # ── 4. Analyse per unique gene strand ───────────────────────────────────
    unique_strands: set[int] = {
        1 if g.get("strand", "+") in (1, "+") else -1
        for g in raw_genes
    }
    if not unique_strands:
        unique_strands = {1}  # default to + strand if no gene data

    alt_allele = alts[0]  # analyse first alt allele

    for strand_int in sorted(unique_strands):
        strand_char = "+" if strand_int == 1 else "-"

        if strand_int == 1:
            ref_seq = ref_seq_plus
            snp_idx = snp_idx_plus
            ref_allele_local = ref
            alt_allele_local = alt_allele
        else:
            ref_seq = reverse_complement(ref_seq_plus)
            # mirror the SNP index
            snp_idx = len(ref_seq_plus) - snp_idx_plus - 1
            ref_allele_local = complement(ref)
            alt_allele_local = complement(alt_allele)

        # Sanity-check the reference base at the SNP position
        actual_base = ref_seq[snp_idx].upper() if snp_idx < len(ref_seq) else "?"
        seq_mismatch = actual_base != ref_allele_local.upper()

        # Build alt sequence
        alt_seq = ref_seq[:snp_idx] + alt_allele_local + ref_seq[snp_idx + 1:]

        ref_matches = find_seed_matches(ref_seq)
        alt_matches = find_seed_matches(alt_seq)

        lof, gof, neutral = _diff_matches(ref_matches, alt_matches)

        # Keep only changes within ±50 nt of the SNP (the seed must overlap
        # the substituted position or be very close to it)
        lof = [m for m in lof if abs(m["start"] - snp_idx) <= FILTER_DISTANCE]
        gof = [m for m in gof if abs(m["start"] - snp_idx) <= FILTER_DISTANCE]

        result["analyses"].append(
            {
                "strand": strand_char,
                "ref_seq": ref_seq,
                "alt_seq": alt_seq,
                "snp_idx": snp_idx,
                "ref_allele": ref_allele_local,
                "alt_allele": alt_allele_local,
                "lof": sorted(lof, key=lambda m: m["start"]),
                "gof": sorted(gof, key=lambda m: m["start"]),
                "neutral_count": len(neutral),
                "seq_mismatch": seq_mismatch,
            }
        )

    return result
