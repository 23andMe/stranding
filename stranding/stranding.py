import logging

from Bio.pairwise2 import align, format_alignment
from Bio.Seq import Seq

from seqseek import Chromosome

from .exceptions import (MissingReferenceFlank,
                        InconsistentAlignment,
                        Unstrandable,
                        FlanksTooShort)


logger = logging.getLogger("stranding")


DEFAULT_MIN_FLANK_LENGTH = 10
DEFAULT_WINDOW_EXTENSION = 20
DEFAULT_SCORE_MULTIPLIER = 1.5
DEFAULT_IDENTITY_CUTOFF_RATIO = 0.90
DEFAULT_MATCH_SCORE = 2
DEFAULT_MISMATCH_PENALTY = -1


class GenomeStranding(object):

    def __init__(self,
                 min_flank_length=DEFAULT_MIN_FLANK_LENGTH,
                 score_multiplier=DEFAULT_SCORE_MULTIPLIER,
                 identity_cutoff_ratio=DEFAULT_IDENTITY_CUTOFF_RATIO,
                 match_score=DEFAULT_MATCH_SCORE,
                 mismatch_penalty=DEFAULT_MISMATCH_PENALTY):

            self.min_flank_length = min_flank_length
            self.score_multiplier = score_multiplier
            self.identity_cutoff_ratio = identity_cutoff_ratio
            self.match_score = match_score
            self.mismatch_penalty = mismatch_penalty
            #ssw_scoring = swalign.NucleotideScoringMatrix(self.match_score, self.mismatch_penalty)
            #self.aligner = swalign.LocalAlignment(ssw_scoring)

    def is_high_scoring(self, alignment, query):
        a1, a2, score, begin, end = alignment
        return score > len(query) * self.score_multiplier

    def align(self, a, b):
        return align.localms(a, b, self.match_score, self.mismatch_penalty, self.mismatch_penalty, self.mismatch_penalty)

    def strand_flanks(self, _5p, _3p, build, chr_name, pos, window=DEFAULT_WINDOW_EXTENSION):
        """
        Given a 5' flank, a 3' flank, the assembly, chromosome, and position
        determine whether a strand flip is required to get to forward genome stranding.

        Consider BLAT if mapping information is unknown.
        """

        # sanity checks
        if pos == 0:
            raise Unstrandable('Position 0 is unmapped')
        elif chr_name in ('0', 0):
            raise Unstrandable('Chromosome 0 is unmapped')
        elif max(_5p, _3p) < self.min_flank_length:
            raise FlanksTooShort('Min flank length is specified as %d' % self.min_flank_length)

        # chromosome-specific conventions
        loop = chr_name == 'MT'
        chr_name = chr_name if chr_name != 'XY' else 'X'

        # reference sequences
        try:
            chromosome = Chromosome(chr_name, build, loop=loop)
            ref_5p = chromosome.sequence(pos - window - len(_5p), pos + window)
            ref_3p = chromosome.sequence(pos - window, pos + len(_3p) + window)
        except ValueError:
            raise MissingReferenceFlank()

        ref_5p_RC = str(Seq(ref_5p).reverse_complement())
        ref_3p_RC = str(Seq(ref_3p).reverse_complement())

        # alignments named w.r.t the reference sequence
        fwd_5p_alignments = self.align(ref_5p, _5p)
        fwd_3p_alignments = self.align(ref_3p, _3p)
        rev_5p_alignments = self.align(ref_5p_RC, _3p)
        rev_3p_alignments = self.align(ref_3p_RC, _5p)

        good_fwd_5p_alignments = [a for a in fwd_5p_alignments if self.is_high_scoring(a, _5p)]
        good_fwd_3p_alignments = [a for a in fwd_3p_alignments if self.is_high_scoring(a, _3p)]
        good_rev_5p_alignments = [a for a in rev_5p_alignments if self.is_high_scoring(a, _3p)]
        good_rev_3p_alignments = [a for a in rev_3p_alignments if self.is_high_scoring(a, _5p)]

        good_fwd_alignments = good_fwd_5p_alignments + good_fwd_3p_alignments
        good_rev_alignments = good_rev_5p_alignments + good_rev_3p_alignments

        if good_fwd_alignments and good_rev_alignments:
            raise InconsistentAlignment()
        elif good_fwd_alignments:
            return 1
        elif good_rev_alignments:
            return -1
        raise Unstrandable()

