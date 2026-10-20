lazy from collections import defaultdict
lazy from itertools import combinations
lazy from pathlib import Path
lazy from . import constants as btc
lazy from . import models as btm
lazy from . import pairwise_align as btpa
lazy from . import parsers as btp


def motifs(protein1:btm.Protein, protein2:btm.Protein) -> dict:
    '''
    Given two protein sequences, this method looks for similar stretches fo amino acids, ranging from length 4 to length 40.

    Matches count as +1, while conservative substitutions count for +0.5. Any two peptides where total_score / length > threshold is considered a match.

    Threshold is more stringent for shorter peptides, and more lenient for longer peptides. Defined in constant THRESHOLD.
    '''

    # Determine shortest sequence and set that as seq1
    length1 = len(protein1.sequence)
    length2 = len(protein2.sequence)

    if length1 > length2:
        seq1 = protein2.sequence
        seq2 = protein1.sequence
    else:
        seq1 = protein1.sequence
        seq2 = protein2.sequence

    # Initialize the dictionary that will hold matches
    similars = defaultdict(list)

    # Loop that looks at peptides of length 4 to 40 in the two sequences, compares them, and looks for regions of similarity.
    for size in range(4, min(40, length1, length2)+1):

        for min_size, max_size, cutoff in btc.THRESHOLD:
            if min_size <= size <= max_size:
                threshold = size * cutoff
                break
        else:
            raise ValueError(f"No threshold defined for motif size {size}.")

        for n in range(0, len(seq1) - size + 1):
            peptide1 = seq1[n:n+size]
            for x in range(0, len(seq2) - size + 1):
                count = 0
                peptide2 = seq2[x:x+size]
                for z in range(size):
                    if peptide1[z] == peptide2[z]:
                        count += 1
                        continue
                    if peptide2[z] in btc.CONSERVATIVES[peptide1[z]]:
                        count += 0.5
                if count >= threshold:
                    similars[(size, threshold)].append(((peptide1, n+1, n + size), (peptide2, x+1, x + size), count))

    return similars


def identity(
        file:str
    ) -> dict[tuple[str,str], float]:

    if Path(file).exists():
        sequences = btp.parse_fasta_file(file)
    else:
        sequences = btp.parse_fasta_string(file)

    identities = {}

    for protein1, protein2 in combinations(sequences, 2):
        align1, align2 = btpa.pairwise_protein(protein1, protein2)
        length = len(align1)
        identical = sum(1 for aa1, aa2 in zip(align1, align2, strict=True) if aa1==aa2)
        percent_identity = identical/length
        identities[(protein1.name, protein2.name)] = percent_identity

    return identities