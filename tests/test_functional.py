import os
import logging
import sys
from unittest import TestCase

from seqseek import BUILD37, BUILD38

from stranding import stranding

# these are the sequences in the seqseek fixtures for chromosme 1 on build 37 and 38
# they are reverse complements of each other
hg19_seq = ('CGGCTGTCCAAGGAGCTGCAGGCGGCGCAGGCCCGGCTGGGCGCGGACATGGAGGACGTGT'
            'GCGGCCGCCTGGTGCAGTACCGCGGCGAGGTGCAGGCCATGCTCGGCCAGAGCACCGAGG')
hg38_seq = ('CCTCGGTGCTCTGGCCGAGCATGGCCTGCACCTCGCCGCGGTACTGCACCAGGCGGCCGCA'
            'CACGTCCTCCATGTCCGCGCCCAGCCGGGCCTGCGCCGCCTGCAGCTCCTTGGACAGCCG')


class SeqSeekTestCase(TestCase):
    def setUp(self):
        os.environ["SEQSEEK_DATA_DIR"] = os.path.join('tests', 'seqseek_fixtures')
        self.strander = stranding.GenomeStranding()
        self.strict_strander = stranding.GenomeStranding(tolerance=1.0)

        console = logging.StreamHandler(stream=sys.stderr)
        root = logging.getLogger('')
        root.setLevel('DEBUG')
        if not root.handlers:
            root.addHandler(console)


class TestBasicStranding(SeqSeekTestCase):
    PERFECT_5P = 'GCGCGGACATGGAGGACGTG'
    PERFECT_3P = 'GCGGCCGCCTGGTGCAGTAC'

    def test_exact_forward_alignment_5p(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert 1 == self.strict_strander.strand_flanks(_5p, '', BUILD37, 1, 60)

    def test_exact_forward_alignment_3p(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert 1 == self.strict_strander.strand_flanks('', _3p, BUILD37, 1, 60)

    def test_exact_reverse_alignment_5p(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert -1 == self.strict_strander.strand_flanks(_5p, '', BUILD38, 1, 60)

    def test_exact_reverse_alignment_3p(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert -1 == self.strict_strander.strand_flanks('', _3p, BUILD38, 1, 60)

    def test_perfect_score_forward_alignment_5p(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert 1 == self.strict_strander.strand_flanks(_5p, '', BUILD37, 1, 61, window=1)

    def test_perfect_score_forward_alignment_3p(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert 1 == self.strict_strander.strand_flanks('', _3p, BUILD37, 1, 61, window=1)

    def test_perfect_score_reverse_alignment_5p(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert -1 == self.strict_strander.strand_flanks(_5p, '', BUILD38, 1, 61, window=1)

    def test_offset_forward_alignment(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert 1 == self.strander.strand_flanks(_5p, _3p, BUILD37, 1, 50, window=10)

    def test_offset_reverse_alignment(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert -1 == self.strander.strand_flanks(_5p, _3p, BUILD38, 1, 50, window=10)

    def test_alignment_beyond_offset_boundary(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        with self.assertRaises(stranding.Unstrandable):
            self.strander.strand_flanks(_5p, _3p, BUILD38, 1, 50, window=5)

    def test_alignment_too_close_to_contig_boundary(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        with self.assertRaises(stranding.MissingReferenceFlank):
            self.strander.strand_flanks(_5p, _3p, BUILD38, 1, 10)

    def test_strict_stranding(self):
        _5p = self.PERFECT_5P + 'A'
        _3p = self.PERFECT_3P + 'A'
        strander = stranding.GenomeStranding(tolerance=1.0)
        with self.assertRaises(stranding.Unstrandable):
            self.assertEqual(1, strander.strand_flanks(_5p, _3p, BUILD37, 1, 60))

    def test_chromosome_0(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        try:
            self.strander.strand_flanks(_5p, _3p, BUILD37, 0, 60, window=10)
        except stranding.Unstrandable as e:
            assert str(e) == 'Chromosome 0 is unmapped'

    def test_position_0(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        try:
            self.strander.strand_flanks(_5p, _3p, BUILD37, 1, 0, window=10)
        except stranding.Unstrandable as e:
            assert str(e) == 'Position 0 is unmapped'

    def test_short_flanks(self):
        with self.assertRaises(stranding.FlanksTooShort):
            self.strander.strand_flanks('atcg', 'atcg', BUILD37, 1, 1)


class TestFuzzyStranding(SeqSeekTestCase):

    def test_forward_fuzzy_stranding(self):
        _5p = 'GCGCAGGCCCGGCTGCGCGCGGTCATGGAGGACGTGT'
        _3p = 'GCGGCCGCCTGGTCGGCAGTAC'
        assert 1 == self.strander.strand_flanks(_5p, _3p, BUILD37, 1, 60, window=10)
        assert -1 == self.strander.strand_flanks(_5p, _3p, BUILD38, 1, 60, window=10)

    def test_inconsistent_stranding(self):
        _5p = "GCGGCCGCT"  # this is present on both the forward and reverse strand in the same region
        strander = stranding.GenomeStranding(min_flank_length=7)
        with self.assertRaises(stranding.InconsistentAlignment):
            strander.strand_flanks(_5p, '', BUILD37, 1, 65, window=20)

