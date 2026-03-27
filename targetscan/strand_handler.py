"""Strand-aware nucleotide handling.

For genes on the **minus (−) strand**:

* The reported ref/alt alleles must be complemented
  (e.g. ``G>A`` becomes ``C>T``).
* The SNP position within the UTR is counted from the **end** of the UTR
  (``end − snp_pos``), and the change is applied at that offset from the
  **beginning** of the stored sequence.

For genes on the **plus (+) strand** the alleles are used as-is, and the
position is simply ``snp_pos − start``.
"""

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(base):
    """Return the Watson-Crick complement of *base*."""
    return base.translate(_COMPLEMENT)


def snp_offset_in_utr(snp_pos, utr_start, utr_end, strand):
    """Return the 0-based offset of the SNP inside the UTR sequence.

    Parameters
    ----------
    snp_pos : int
        Chromosomal position of the SNP.
    utr_start, utr_end : int
        Chromosomal start/end of the 3'UTR.
    strand : str
        ``"+"`` or ``"-"``.

    Returns
    -------
    int
        0-based index into the UTR sequence string.
    """
    if strand == "-":
        return utr_end - snp_pos
    return snp_pos - utr_start


def apply_snp(sequence, offset, alt_base, strand):
    """Return *sequence* with the nucleotide at *offset* replaced.

    If the gene is on the minus strand the *alt_base* is first
    complemented before substitution.
    """
    if strand == "-":
        alt_base = complement(alt_base)
    seq_list = list(sequence)
    seq_list[offset] = alt_base.upper()
    return "".join(seq_list)


def effective_alleles(ref, alt, strand):
    """Return ``(eff_ref, eff_alt)`` adjusted for strand.

    On the minus strand both alleles are complemented.
    """
    if strand == "-":
        return complement(ref), complement(alt)
    return ref, alt
