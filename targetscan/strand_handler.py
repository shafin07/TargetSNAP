"""
strand_handler.py – DNA strand utilities for TargetSNAP.

For minus-strand genes the mRNA sequence is the reverse complement of the
genomic + strand.  Helper functions here let the rest of the code stay
strand-agnostic.
"""

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def complement(seq: str) -> str:
    """Return the complement of a DNA sequence (same 5′→3′ direction)."""
    return seq.translate(_COMP)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (5′→3′ of opposite strand)."""
    return complement(seq)[::-1]


def snp_offset_in_utr(snp_pos: int, utr_start: int, utr_end: int, strand: str) -> int:
    """Return the 0-based offset of a SNP inside a UTR segment.

    Parameters
    ----------
    snp_pos   : 1-based genomic coordinate of the SNP.
    utr_start : 1-based start of the UTR (always < utr_end, + strand convention).
    utr_end   : 1-based end of the UTR.
    strand    : '+' or '-'.

    Returns
    -------
    0-based index into the UTR sequence *in mRNA orientation*.
    """
    if strand == "+":
        return snp_pos - utr_start
    else:
        # On the minus strand the mRNA reads from utr_end down to utr_start;
        # position utr_end is index 0 in the RC sequence.
        return utr_end - snp_pos
