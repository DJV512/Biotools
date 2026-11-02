def _align_protein(
        seq1: str,
        seq2: str
    ) -> tuple[str, str]:
    """
    Perform pairwise protein sequence alignment using the Gotoh algorithm.

    Uses BLOSUM62 substitution scores and affine gap penalties.

    Implemented in Rust. For use only with validated Protein objects.

    Args:
        seq1: First protein sequence.
        seq2: Second protein sequence.

    Returns:
        A tuple containing the two aligned protein sequences,
        including gap characters.
    """
    ...

def _merge(
        left_names:list[str], 
        left_sequences: list[str],
        right_names:list[str],
        right_sequences: list[str],
        gap_open:float,
        gap_extend:float,
    ) -> tuple[list[str], list[str]]:
    '''
    Globally align two profiles, made up of already aligned sequences, using affine gap penalties.

    Inputs are two Profiles to merge (separated out into their constituent lists of names and sequences to pass to rust),
    as well as the gap open and gap extend penalties to use. This internal function does not perform input validation.

    Returns a tuple of a list of sequence names, and a list of sequence alignments.
    '''
    ...