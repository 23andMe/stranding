import logging

import swalign

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
            ssw_scoring = swalign.NucleotideScoringMatrix(self.match_score, self.mismatch_penalty)
            self.aligner = swalign.LocalAlignment(ssw_scoring)

    def is_high_scoring(self, alignment):
        print self.identity_cutoff_ratio
        print self.score_multiplier
        return (alignment.identity > self.identity_cutoff_ratio and
                alignment.score > len(alignment.query) * self.score_multiplier)


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

        # alignments named w.r.t the reference sequence
        fwd_5p_alignment = self.aligner.align(ref_5p, _5p)
        fwd_3p_alignment = self.aligner.align(ref_3p, _3p)
        rev_5p_alignment = self.aligner.align(swalign.revcomp(ref_5p), _3p, rc=True)
        rev_3p_alignment = self.aligner.align(swalign.revcomp(ref_3p), _5p, rc=True)

        alignments = [fwd_5p_alignment, fwd_3p_alignment, rev_5p_alignment, rev_3p_alignment]
        high_scoring_alignments = [a for a in alignments if self.is_high_scoring(a)]
        strands = [-1 if a.rc else 1 for a in high_scoring_alignments]

        if not high_scoring_alignments:
            for alignment in alignments:
                alignment.dump()
            raise Unstrandable()

        if not len(set(strands)) == 1:
            raise InconsistentAlignment()

        return strands[0]
