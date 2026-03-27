"""Gain-of-Function / Loss-of-Function analysis.

Compare the TargetScan Context++ scores between a reference (ref) and
an alternate (alt) UTR sequence to classify each miRNA family as:

* **Loss of Function (LOF)** – the miRNA binding site is present with
  the reference allele but absent or weaker with the alternate allele.
* **Gain of Function (GOF)** – the miRNA binding site is absent or
  weaker with the reference allele but present or stronger with the
  alternate allele.
* **Unchanged** – no meaningful change in binding affinity.
"""


def _score(row):
    """Extract a numeric context/score value from a result row."""
    for key in ("context_score", "context++_score", "weighted_context++_score"):
        val = row.get(key)
        if val is not None:
            try:
                return float(val)
            except (ValueError, TypeError):
                continue
    return 0.0


def _index_by_mirna(results):
    """Index a list of TargetScan result rows by miRNA family.

    Returns a dict mapping ``miRNA_family`` → ``list[row]``.
    """
    index = {}
    for row in results:
        family = row.get("miRNA_family_ID") or row.get("miR_Family", "")
        index.setdefault(family, []).append(row)
    return index


def classify(ref_results, alt_results, threshold=0.0):
    """Classify miRNA families as GOF, LOF, or unchanged.

    Parameters
    ----------
    ref_results, alt_results : list[dict]
        Output from :func:`targetscan_runner.run_targetscan` for the
        reference and alternate UTR sequences respectively.
    threshold : float
        Minimum absolute score difference to be considered meaningful.

    Returns
    -------
    dict
        ``{"gof": [...], "lof": [...], "unchanged": [...]}``.
        Each entry is a dict with ``miRNA_family``, ``ref_score``,
        ``alt_score``, ``delta``, ``classification``, ``ref_details``,
        and ``alt_details``.
    """
    ref_idx = _index_by_mirna(ref_results)
    alt_idx = _index_by_mirna(alt_results)

    all_families = sorted(set(ref_idx) | set(alt_idx))

    gof, lof, unchanged = [], [], []

    for family in all_families:
        ref_rows = ref_idx.get(family, [])
        alt_rows = alt_idx.get(family, [])

        ref_best = min((_score(r) for r in ref_rows), default=0.0)
        alt_best = min((_score(r) for r in alt_rows), default=0.0)
        delta = alt_best - ref_best

        entry = {
            "miRNA_family": family,
            "ref_score": ref_best,
            "alt_score": alt_best,
            "delta": round(delta, 4),
            "ref_details": ref_rows,
            "alt_details": alt_rows,
        }

        if ref_rows and not alt_rows:
            entry["classification"] = "LOF"
            lof.append(entry)
        elif alt_rows and not ref_rows:
            entry["classification"] = "GOF"
            gof.append(entry)
        elif delta < -threshold:
            entry["classification"] = "GOF"
            gof.append(entry)
        elif delta > threshold:
            entry["classification"] = "LOF"
            lof.append(entry)
        else:
            entry["classification"] = "Unchanged"
            unchanged.append(entry)

    return {"gof": gof, "lof": lof, "unchanged": unchanged}
