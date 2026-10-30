from . import constants as btc
from . import decorators as btd
from . import models as btm
from math import inf
from operator import itemgetter

def _sigma_dna(str1:str, str2:str) -> int:
    '''
    Private method for determining a sigma score between two DNA bases.
    '''
    if str1 == str2:
        return btc.DNA_MATCH
    else:
        return btc.DNA_MISMATCH


def _global_dna_backtrack(sequence1:str, sequence2:str, score:dict, i:int, j:int, align1:list[str]=None, align2:list[str]=None) -> tuple[list[str], list[str]]:
    '''
    Private recursive method for taking a pairwise score table and backtracking to assemble the optimal alignment.
    '''

    if align1 == None and align2 == None:
        align1 = []
        align2 = []


    #Base Cases
    if i == 0 or j == 0:

        if i > 0:
            for x in range(1, i+1):
                align1.append(sequence1[x])
                align2.append("-")

        if j > 0:
            for x in range(1, j+1):
                align1.append("-")
                align2.append(sequence2[x])

        return align1, align2
        
    # Recursive calls based on score matching
    if score[(i,j)] == score[(i-1, j-1)] + _sigma_dna(sequence1[i], sequence2[j]):
        _global_dna_backtrack(sequence1, sequence2, score, i-1, j-1, align1, align2)
        align1.append(sequence1[i])
        align2.append(sequence2[j])
        return align1, align2


    elif score[(i,j)] == score[(i-1, j)] + btc.DNA_GAP:
        _global_dna_backtrack(sequence1, sequence2, score, i-1, j, align1, align2)
        align1.append(sequence1[i])
        align2.append("-")
        return align1, align2


    else:
        _global_dna_backtrack(sequence1, sequence2, score, i, j-1, align1, align2)
        align1.append("-")
        align2.append(sequence2[j])
        return align1, align2
 

@btd.validate_dna()
def pairwise_dna(dna1:str, dna2:str) -> tuple[list[str], list[str]]:
    '''
    Performs a pairwise alignment of two DNA sequences.\n
    Maximum sequence length for alignment ~ 1000 bases.\n
    Required arguments: dna1, dna2
    '''

    #Add a space at the beginning to represent the 0th position
    dna1 = " " + dna1.upper()
    dna2 = " " + dna2.upper()
    length1 = len(dna1)
    length2 = len(dna2)

    #Pre-populate the score dicitonary with the horizontal 0 row
    score = {(0,j):btc.DNA_GAP*j for j in range(length2)}

    #Iterate over the entire length of sequence 1
    for i in range(1,length1):

        #Pre-populate the score matrix with the vertical 0 row
        score[i,0] = btc.DNA_GAP * i

        #Iterate over the entire length of sequence 2
        for j in range(1, length2):

            #Choose best score and add it to the score matrix
            score[i,j] = max(
                score[i - 1 ,j] + btc.DNA_GAP,
                score[i, j - 1] + btc.DNA_GAP,
                score[i - 1, j - 1] + _sigma_dna(dna1[i], dna2[j])
            ) 

    return _global_dna_backtrack(dna1, dna2, score, length1-1, length2-1)


def _global_protein_backtrack(
        sequence1:str, 
        sequence2:str, 
        score_diag:dict[tuple[int, int], tuple[float, str]], 
        score_hori:dict[tuple[int, int], tuple[float, str]], 
        score_vert:dict[tuple[int, int], tuple[float, str]], 
        i:int, 
        j:int
    ) -> tuple[list[str], list[str]]:
    '''
    Private method for taking pairwise score matrices and backtracking through them to assemble the optimal alignment.
    '''

    align1 = []
    align2 = []

    _, predecessor = max(
        (score_diag[(i, j)][0], "diag"),
        (score_hori[(i, j)][0], "hori"),
        (score_vert[(i, j)][0], "vert"),
        key = itemgetter(0)
    )

    while i > 0 or j > 0:

        if predecessor == "diag":
            predecessor = score_diag[(i, j)][1]
            align1.append(sequence1[i])
            align2.append(sequence2[j])
            i -= 1
            j -= 1

        elif predecessor == "hori":
            predecessor = score_hori[(i, j)][1]
            align1.append("-")
            align2.append(sequence2[j])
            j -= 1

        else:
            predecessor = score_vert[(i, j)][1]
            align1.append(sequence1[i])
            align2.append("-")
            i -= 1

    align1.reverse()
    align2.reverse()
    
    return align1, align2
    
        
 

def pairwise_protein(protein1:btm.Protein, protein2:btm.Protein) -> tuple[list[str], list[str]]:
    '''
    Performs a pairwise alignment of two amino acid sequences.

    Residue-residue scoring is determined by BLOSUM62.

    Affine gap scoring: GAP_OPEN costs 10, while GAP_EXTEND only costs 1.
    '''

    #Add a space at the beginning to represent the 0th position
    seq1 = " " + protein1.sequence
    seq2 = " " + protein2.sequence
    length1 = len(seq1)
    length2 = len(seq2)


    #Initialize the affine gap scoring matrices
    score_diag = {
        (0, 0): (0, "start")
    }
    score_vert = {
        (0, 0): (-inf, "impossible") 
    }
    score_hori = {
        (0, 0): (-inf, "impossible")
    }
    for j in range(1, length2):
        score_diag[(0, j)] = (-inf, "impossible")
        score_vert[(0, j)] = (-inf, "impossible")
        score_hori[(0, j)] = (btc.PROTEIN_GAP_OPEN + (j-1) * btc.PROTEIN_GAP_EXTEND, "diag" if j == 1 else "hori")

    for i in range(1, length1):
        score_diag[(i, 0)] = (-inf, "impossible")
        score_vert[(i, 0)] = (btc.PROTEIN_GAP_OPEN + (i-1) * btc.PROTEIN_GAP_EXTEND, "diag" if i == 1 else "vert")
        score_hori[(i, 0)] = (-inf, "impossible")


    #Iterate over the entire length of sequence 1
    for i in range(1,length1):

        #Iterate over the entire length of sequence 2
        for j in range(1, length2):

            #Choose best score for each of the three scoring matrices and add them
            #Current pair is a residue-residue alignment
            residue_blosum_score = btc.BIOTOOLS_COMPARE[seq1[i]][seq2[j]]
            score_diag[(i, j)] = max(
                (score_diag[(i - 1, j - 1)][0] + residue_blosum_score, "diag"),
                (score_hori[(i - 1, j - 1)][0] + residue_blosum_score, "hori"),
                (score_vert[(i - 1, j - 1)][0] + residue_blosum_score, "vert"),
                key = itemgetter(0)
            )
            
            #Current pair is a gap in sequence 2
            score_vert[(i, j)] = max(
                (score_diag[(i - 1, j)][0] + btc.PROTEIN_GAP_OPEN, "diag"),
                (score_vert[(i - 1, j)][0] + btc.PROTEIN_GAP_EXTEND, "vert"),
                key = itemgetter(0)
            )

            #Current pair is a gap in sequence 1
            score_hori[(i, j)] = max(
                (score_diag[(i, j - 1)][0] + btc.PROTEIN_GAP_OPEN, "diag"),
                (score_hori[(i, j - 1)][0] + btc.PROTEIN_GAP_EXTEND, "hori"),
                key = itemgetter(0)
            )

    return _global_protein_backtrack(seq1, seq2, score_diag, score_hori, score_vert, length1-1, length2-1)