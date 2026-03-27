"""Unit tests for TargetSNAP core modules."""

import json
import unittest

from targetscan.analysis import classify
from targetscan.strand_handler import (
    apply_snp,
    complement,
    effective_alleles,
    snp_offset_in_utr,
)


class TestComplement(unittest.TestCase):
    def test_single_bases(self):
        self.assertEqual(complement("A"), "T")
        self.assertEqual(complement("T"), "A")
        self.assertEqual(complement("C"), "G")
        self.assertEqual(complement("G"), "C")

    def test_lowercase(self):
        self.assertEqual(complement("a"), "t")
        self.assertEqual(complement("g"), "c")


class TestSnpOffset(unittest.TestCase):
    def test_plus_strand(self):
        # SNP at pos 100, UTR 90-110 → offset = 100 - 90 = 10
        self.assertEqual(snp_offset_in_utr(100, 90, 110, "+"), 10)

    def test_minus_strand(self):
        # SNP at 98517568, UTR 98514785-98518215 → offset = 98518215 - 98517568 = 647
        self.assertEqual(
            snp_offset_in_utr(98517568, 98514785, 98518215, "-"), 647
        )


class TestApplySnp(unittest.TestCase):
    def test_plus_strand_no_complement(self):
        seq = "AACCGGTT"
        result = apply_snp(seq, 2, "T", "+")
        self.assertEqual(result, "AATCGGTT")

    def test_minus_strand_complement(self):
        seq = "AACCGGTT"
        # alt_allele='A' → complement='T' for minus strand
        result = apply_snp(seq, 2, "A", "-")
        self.assertEqual(result, "AATCGGTT")


class TestEffectiveAlleles(unittest.TestCase):
    def test_plus_strand(self):
        self.assertEqual(effective_alleles("G", "A", "+"), ("G", "A"))

    def test_minus_strand(self):
        # G>A on minus strand → C>T
        self.assertEqual(effective_alleles("G", "A", "-"), ("C", "T"))


class TestClassify(unittest.TestCase):
    def _row(self, family, score):
        return {"miRNA_family_ID": family, "context_score": str(score)}

    def test_lof_site_lost(self):
        ref = [self._row("miR-21", -0.3)]
        alt = []
        result = classify(ref, alt)
        self.assertEqual(len(result["lof"]), 1)
        self.assertEqual(result["lof"][0]["classification"], "LOF")

    def test_gof_site_gained(self):
        ref = []
        alt = [self._row("miR-21", -0.3)]
        result = classify(ref, alt)
        self.assertEqual(len(result["gof"]), 1)
        self.assertEqual(result["gof"][0]["classification"], "GOF")

    def test_unchanged(self):
        ref = [self._row("miR-21", -0.3)]
        alt = [self._row("miR-21", -0.3)]
        result = classify(ref, alt)
        self.assertEqual(len(result["unchanged"]), 1)

    def test_gof_by_score_decrease(self):
        # alt score more negative → stronger binding → GOF
        ref = [self._row("miR-21", -0.1)]
        alt = [self._row("miR-21", -0.5)]
        result = classify(ref, alt)
        self.assertEqual(len(result["gof"]), 1)

    def test_lof_by_score_increase(self):
        # alt score less negative → weaker binding → LOF
        ref = [self._row("miR-21", -0.5)]
        alt = [self._row("miR-21", -0.1)]
        result = classify(ref, alt)
        self.assertEqual(len(result["lof"]), 1)


class TestFlaskApp(unittest.TestCase):
    """Smoke-test the Flask endpoints (without external API calls)."""

    def setUp(self):
        from app import app
        app.config["TESTING"] = True
        self.client = app.test_client()

    def test_index_page(self):
        resp = self.client.get("/")
        self.assertEqual(resp.status_code, 200)
        self.assertIn(b"TargetSNAP", resp.data)


if __name__ == "__main__":
    unittest.main()
