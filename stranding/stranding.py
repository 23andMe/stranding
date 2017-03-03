import logging

import swalign

from seqseek import Chromosome

from .exceptions import (MissingReferenceFlank,
                        InconsistentAlignment,
                        Unstrandable,
                        FlanksTooShort)


logger = logging.getLogger("stranding")


MIN_FLANK_LENGTH = 10
DEFAULT_MAX_OFFSET = 20
SCORE_MULTIPLIER = 1.5
IDENTITY_CUTOFF_RATIO = 0.90
MATCH_SCORE = 2
MISMATCH_PENALTY = -1
SSW_SCORING = swalign.NucleotideScoringMatrix(MATCH_SCORE, MISMATCH_PENALTY)
ALIGNER = swalign.LocalAlignment(SSW_SCORING)


def is_high_scoring(alignment):
    return (alignment.identity > IDENTITY_CUTOFF_RATIO and
            alignment.score > len(alignment.query) * SCORE_MULTIPLIER)


def genome_stranding(_5p, _3p, build, chr_name, pos, max_offset=DEFAULT_MAX_OFFSET):
    """
    Given a 5' flank, a 3' flank, the assembly, chromosome, and position
    determine whether a strand flip is required to get to forward genome stranding.

    A better way of doing this would be to integrate something like BLAT which can
    quickly map arbitrary sequences in a genome assembly but has setup and licencing issues.
    """
    # sanity checks
    if pos == 0:
        raise Unstrandable('Position 0 is unmapped')
    elif chr_name in ('0', 0):
        raise Unstrandable('Chromosome 0 is unmapped')
    elif max(_5p, _3p) < MIN_FLANK_LENGTH:
        raise FlanksTooShort('Min flank length is specified as 10')

    # chromosome-specific conventions
    loop = chr_name == 'MT'
    chr_name = chr_name if chr_name != 'XY' else 'X'

    # reference sequences
    try:
        chromosome = Chromosome(chr_name, build, loop=loop)
        ref_5p = chromosome.sequence(pos - max_offset - len(_5p), pos + max_offset)
        ref_3p = chromosome.sequence(pos - max_offset, pos + len(_3p) + max_offset)
    except ValueError:
        raise MissingReferenceFlank()

    # alignments named w.r.t the reference sequence
    fwd_5p_alignment = ALIGNER.align(ref_5p, _5p)
    fwd_3p_alignment = ALIGNER.align(ref_3p, _3p)
    rev_5p_alignment = ALIGNER.align(swalign.revcomp(ref_5p), _3p, rc=True)
    rev_3p_alignment = ALIGNER.align(swalign.revcomp(ref_3p), _5p, rc=True)

    alignments = [fwd_5p_alignment, fwd_3p_alignment, rev_5p_alignment, rev_3p_alignment]
    high_scoring_alignments = [a for a in alignments if is_high_scoring(a)]
    strands = [-1 if a.rc else 1 for a in high_scoring_alignments]

    for alignment in high_scoring_alignments:
        logger.info(alignment.dump())

    if not high_scoring_alignments:
        for alignment in alignments:
            logger.debug(alignment.query)
            logger.debug(alignment.ref)
            logger.debug(alignment.dump())
        raise Unstrandable()

    if not len(set(strands)) == 1:
        raise InconsistentAlignment()

    return strands[0]
