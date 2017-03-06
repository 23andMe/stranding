import os
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


class TestBasicStranding(SeqSeekTestCase):
    PERFECT_5P = 'GCGCGGACATGGAGGACGTG'
    PERFECT_3P = 'GCGGCCGCCTGGTGCAGTAC'

    def test_perfect_forward_alignment(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert 1 == self.strander.strand_flanks(_5p, _3p, BUILD37, 1, 60)
        assert 1 == self.strander.strand_flanks('', _3p, BUILD37, 1, 60)

    def test_perfect_reverse_alignment(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert -1 == self.strander.strand_flanks('', _3p, BUILD38, 1, 60, window=10)

    def test_perfect_reverse_complement_3p(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert -1 == self.strander.strand_flanks(_5p, '', BUILD38, 1, 52, window=10)

    def test_offset_forward_alignment(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert 1 == self.strander.strand_flanks(_5p, _3p, BUILD37, 1, 50)

    def test_offset_reverse_alignment(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        assert -1 == self.strander.strand_flanks(_5p, _3p, BUILD38, 1, 50)

    def test_alignment_beyond_offset_boundary(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        with self.assertRaises(stranding.Unstrandable):
            self.strander.strand_flanks(_5p, _3p, BUILD38, 1, 40, window=10)

    def test_alignment_too_close_to_contig_boundary(self):
        _5p, _3p = self.PERFECT_5P, self.PERFECT_3P
        with self.assertRaises(stranding.MissingReferenceFlank):
            self.strander.strand_flanks(_5p, _3p, BUILD38, 1, 10)

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


class TestFuzzyStranding(SeqSeekTestCase):

    def test_fuzzy_stranding(self):
        _5p = 'GCGCAGGCCCGGCTGCGCGCGGTCATGGAGGACGTGT'
        _3p = 'tttttttttttttttttttttt'
        assert 1 == self.strander.strand_flanks(_5p, _3p, BUILD37, 1, 60, window=10)
        assert -1 == self.strander.strand_flanks(_5p, _3p, BUILD38, 1, 60, window=10)

    def test_inconsistent_stranding(self):
        _5p = "CGCCTGA"  # this is present on both the forward and reverse strand in the same region
        with self.assertRaises(stranding.InconsistentAlignment):
            self.strander.strand_flanks(_5p, '', BUILD37, 1, 80)

