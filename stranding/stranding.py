import logging
import warnings

from Bio.pairwise2 import align, format_alignment
from Bio.Seq import Seq

from seqseek import Chromosome

from .exceptions import (MissingReferenceFlank,
                        InconsistentAlignment,
                        Unstrandable,
                        FlanksTooShort)


LOGGER = logging.getLogger("stranding")
DEFAULT_MIN_FLANK_LENGTH = 15
DEFAULT_WINDOW_EXTENSION = 0
DEFAULT_IDENTITY_CUTOFF_RATIO = 0.77
DEFAULT_MATCH_SCORE = 2
DEFAULT_MISMATCH_PENALTY = -1
DEFAULT_GAP_OPEN_PENALTY = -5


class GenomeStranding(object):

    def __init__(self,
                 min_flank_length=DEFAULT_MIN_FLANK_LENGTH,
                 identity_cutoff_ratio=DEFAULT_IDENTITY_CUTOFF_RATIO,
                 match_score=DEFAULT_MATCH_SCORE,
                 mismatch_penalty=DEFAULT_MISMATCH_PENALTY,
                 gap_open_penalty=DEFAULT_GAP_OPEN_PENALTY):

            self.min_flank_length = min_flank_length
            self.identity_cutoff_ratio = identity_cutoff_ratio
            self.match_score = match_score
            self.mismatch_penalty = mismatch_penalty
            self.gap_open_penalty = gap_open_penalty
            if self.min_flank_length < DEFAULT_MIN_FLANK_LENGTH:
                warnings.warn('Short flank lengths may lead to inaccurate alignments')



    def is_high_scoring(self, score, query):
        if len(query) < self.min_flank_length:
            return False
        return score > len(query) * self.match_score * self.identity_cutoff_ratio

    def is_perfect_score(self, score, query):
        return score == len(query) * self.match_score and len(query) > self.min_flank_length

    def align(self, ref, query, score_only=True):
        alignment = align.localms(ref, query, self.match_score, self.mismatch_penalty,
                                  self.gap_open_penalty, self.mismatch_penalty,
                                  score_only=score_only)
        if score_only and not alignment:
            return 0
        return alignment

    def align_and_log(self, ref, query):
        alignments = self.align(ref, query, False)
        for a in alignments:
            if self.is_high_scoring(a, query):
                LOGGER.error(format_alignment(*a))

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
        elif max(len(_5p), len(_3p)) < self.min_flank_length:
            raise FlanksTooShort('At least one flank must be longer than the specified'
                                 'minimum flank length of s %d' % self.min_flank_length)

        # chromosome-specific conventions
        loop = chr_name == 'MT'
        chr_name = chr_name if chr_name != 'XY' else 'X'
        max_length = max(len(_5p), len(_3p))

        # reference sequences
        try:
            chromosome = Chromosome(chr_name, build, loop=loop)
            ref_5p = chromosome.sequence(pos - window - max_length, pos + window)
            ref_3p = chromosome.sequence(pos - window, pos + max_length + window + 1)
        except ValueError:
            raise MissingReferenceFlank(
                'Could not find flanks for %s %d %d' % (chr_name, pos, window))

	# cheap
        if window == 0 and (_3p == ref_3p or _5p == ref_5p):
            return 1

        ref_5p_RC = str(Seq(ref_5p).reverse_complement())
        ref_3p_RC = str(Seq(ref_3p).reverse_complement())

        if window == 0 and (_3p == ref_5p_RC or _5p == ref_3p_RC):
            return -1

        if window == 0 and self.identity_cutoff_ratio == 1.0:
            raise Unstrandable('Strict stranding failed')

	# expensive
        fwd_5p_score = self.align(ref_5p, _5p)
        if self.is_perfect_score(fwd_5p_score, _5p):
            return 1

        fwd_3p_score = self.align(ref_3p, _3p)
        if self.is_perfect_score(fwd_3p_score, _3p):
            return 1

        rev_5p_score = self.align(ref_5p_RC, _3p)
        if self.is_perfect_score(rev_5p_score, _3p):
            return -1

        rev_3p_score = self.align(ref_3p_RC, _5p)
        if self.is_perfect_score(rev_3p_score, _5p):
            return -1

        strands = []
        if self.is_high_scoring(fwd_5p_score, _5p):
            strands.append(1)
        if self.is_high_scoring(fwd_3p_score, _3p):
            strands.append(1)
        if self.is_high_scoring(rev_5p_score, _3p):
            strands.append(-1)
        if self.is_high_scoring(rev_3p_score, _5p):
            strands.append(-1)

        if len(set(strands)) > 1:
            LOGGER.error('Forward alignments')
            self.align_and_log(ref_5p, _5p)
            self.align_and_log(ref_3p, _3p)
            LOGGER.error('Reverse alignments')
            self.align_and_log(ref_5p_RC, _3p)
            self.align_and_log(ref_3p_RC, _5p)
            raise InconsistentAlignment('Inconsistent alignments')
        elif strands:
            return strands[0]
        raise Unstrandable('No matching alignments')
