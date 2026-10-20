def _align_protein(seq1: str, seq2: str) -> tuple[str, str]:
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