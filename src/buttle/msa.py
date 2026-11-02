from itertools import combinations
from collections import Counter
from pathlib import Path
from operator import itemgetter
from math import inf
from copy import copy
from . import models as btm
from . import constants as btc
from . import parsers as btp
from . import pairwise_align as btpa
from . import buttle_rust as btrust


def _build_pairwise_alignments(
        sequences: list[btm.Protein],
        backend: str,
    ) -> list[btm.AlignedProteinPair]:
    '''
    Private function for performing pairwise alignments of all pairs of Proteins in the sequences list.
    '''
    pairwise_alignments = []
    for protein1, protein2 in combinations(sequences, 2):

        if backend == "rust":
            align1, align2 = btrust._align_protein(protein1.sequence, protein2.sequence)
            align1 = list(align1)
            align2 = list(align2)
        else:
            align1, align2 = btpa.pairwise_protein(protein1, protein2)

        distance = _distance_score(align1, align2)
        new_alignment = btm.AlignedProteinPair(
            protein1 = protein1,
            protein2 = protein2,
            align1 = align1,
            align2 = align2,
            distance = distance
        )
        pairwise_alignments.append(new_alignment)
    return pairwise_alignments


def _distance_score(
        seq1:list[str], 
        seq2:list[str]
    ) -> float:
    '''
    Private function for determining distance between two protein sequences.
    '''
    length = len(seq1)
    score = 0
    for residue1, residue2 in zip(seq1, seq2):
        if residue1 != "-" and residue2 != "-":
            if residue1 == residue2:
                score += 1
    similarity = score/length
    distance = 1 - similarity
    return distance


def _build_original_distance_matrix(
        sequences:list[btm.Protein], 
        pairwise_alignments:list[btm.AlignedProteinPair]
    ) -> dict[str,dict[str, float]]:
    '''
    Private function for constructing the first distance matrix.
    '''
    distance_matrix = {protein1.name:{protein2.name:0.0 for protein2 in sequences} for protein1 in sequences}
    for alignment in pairwise_alignments:
        distance_matrix[alignment.protein1.name][alignment.protein2.name] = alignment.distance
        distance_matrix[alignment.protein2.name][alignment.protein1.name] = alignment.distance
    return distance_matrix


def _get_members(
        item: btm.Protein | btm.Cluster
    ) -> list[btm.Protein]:
    '''
    Private function for returning a list of all members associated with item. The item could be a Protein object or a Cluster object.
    '''
    if isinstance(item, btm.Protein):
        return [item]
    return item.members


def _get_distance(
        item1:btm.Protein | btm.Cluster,
        item2:btm.Protein | btm.Cluster,
        original_distance_matrix:dict[str, dict[str, float]]
    ) -> float:
    '''
    Private function for determining the distance between two items, each of which may be a Protein Object, or a Cluster Object.
    '''
    members1 = _get_members(item1)
    members2 = _get_members(item2)
    distances = []
    for member1 in members1:
        for member2 in members2:
            distances.append(original_distance_matrix[member1.name][member2.name])
    return sum(distances)/len(distances)


def _get_id(
        item: btm.Protein | btm.Cluster
    ) -> str:

    if isinstance(item, btm.Protein):
        return item.name
    return item.id


def _build_cluster_distance_matrix(
        items:list[btm.Protein | btm.Cluster],
        original_distance_matrix:dict[str, dict[str, float]]
    ) -> dict[str, dict[str, float]]:
    '''
    Private function for rebuilding the distance matrix after every Cluster merge.
    '''
    distance_matrix = {_get_id(item1):{_get_id(item2):0.0 for item2 in items} for item1 in items}
    for item1, item2 in combinations(items, 2):
        distance = _get_distance(item1, item2, original_distance_matrix)
        id1 = _get_id(item1)
        id2 = _get_id(item2)
        distance_matrix[id1][id2] = distance
        distance_matrix[id2][id1] = distance
    return distance_matrix


def _find_closest_pair(
        matrix:dict[str, dict[str, float]]
    ) -> tuple[str, str]:
    '''
    Private function for finding the closest pair in the distance matrix.
    '''
    closest = float(inf)
    for key1 in matrix:
        for key2 in matrix[key1]:
            if key1 != key2 and matrix[key1][key2] < closest:
                best_pair = (key1, key2)
                closest = matrix[key1][key2]
    return best_pair


def _build_guide_tree(
        sequences:list[btm.Protein],
        original_distance_matrix:dict[str, dict[str, float]]
    ) -> btm.Cluster:
    '''
    Private function for building the guide tree before assemling the alignment.
    '''

    cluster_number = 1
    distance_matrix = original_distance_matrix.copy()
    items = sequences.copy()
    while len(items) > 1:

        # Find the two most similar active nodes, and remove those two nodes from the list of items.
        best_pair = _find_closest_pair(distance_matrix)
        left, right = None, None
        for item in items.copy():
            if _get_id(item) == best_pair[0]:
                left = item
                items.remove(item)
                if right is not None:
                    break
            if _get_id(item) == best_pair[1]:
                right = item
                items.remove(item)
                if left is not None:
                    break

        # Merge them into a Cluster
        members = _get_members(left) + _get_members(right)
        tree = btm.Cluster(
            id=f"__cluster_{cluster_number}",
            name=f"{left.name}-{right.name}",
            left=left,
            right=right,
            members=members
        )
        cluster_number += 1

        # Construct a new distance matrix unless we are left with just one Cluster, which signals the end of the tree construction.
        items.append(tree)
        if len(items) > 1:
            distance_matrix = _build_cluster_distance_matrix(items, original_distance_matrix)

    return tree


def _msa_sigma(
        left_residues: list[str],
        right_residues: list[str]
    ) -> float:
    '''
    Private function for determining the average pairwise score between two profile columns.
    '''
    comparisons = len(left_residues) * len(right_residues)
    score = 0
    for left_residue in left_residues:
        for right_residue in right_residues:
            if left_residue == "-" and right_residue == "-":
                continue
            elif left_residue == "-" or right_residue == "-":
                score += btc.PROTEIN_GAP
            else:
                score += btc.BUTTLE_COMPARE[left_residue][right_residue]
    return score/comparisons


def _global_msa_backtrack(
        left:btm.Profile, 
        right:btm.Profile, 
        score_diag:dict[tuple[int, int], tuple[float, str]],
        score_vert:dict[tuple[int, int], tuple[float, str]],
        score_hori:dict[tuple[int, int], tuple[float, str]],
        i:int, 
        j:int,
    ) -> btm.Profile:
    '''
    Private method for taking affine profile-profile score tables and backtracking to assemble the optimal alignment.
    '''

    #Initialize the final_profile
    total_sequences = len(left.alignment)+len(right.alignment)
    final_profile = btm.Profile(name=left.name + right.name, alignment=[[] for _ in range(total_sequences)])

    _, predecessor = max(
        (score_diag[i,j][0], "diag"),
        (score_vert[i,j][0], "vert"),
        (score_hori[i,j][0], "hori"),
        key = itemgetter(0)
    )

    #Base case
    while i > 0 or j > 0:

        if predecessor == "diag":
            predecessor = score_diag[i, j][1]
            for z, alignment in enumerate(left.alignment):
                final_profile.alignment[z].append(alignment[i])
            for z, alignment in enumerate(right.alignment):
                final_profile.alignment[len(left.alignment) + z].append(alignment[j])
            i -= 1
            j -= 1

        elif predecessor == "hori":
            predecessor = score_hori[i, j][1]
            for z, alignment in enumerate(left.alignment):
                final_profile.alignment[z].append("-")
            for z, alignment in enumerate(right.alignment):
                final_profile.alignment[len(left.alignment) + z].append(alignment[j])
            j -= 1

        else:
            predecessor = score_vert[i, j][1]
            for z, alignment in enumerate(left.alignment):
                final_profile.alignment[z].append(alignment[i])
            for z, alignment in enumerate(right.alignment):
                final_profile.alignment[len(left.alignment) + z].append("-")
            i -= 1

    for alignment in final_profile.alignment:
        alignment.reverse()

    return final_profile


def _merge(
        left: btm.Profile,
        right: btm.Profile,
        gap_open: int,
        gap_extend: int,
    ) -> btm.Profile:
    '''
    Private function for merging Proteins and Profiles during multiple alignment.
    '''

    left = btm.Profile(name=left.name, alignment=[[" ", *row] for row in left.alignment])
    right = btm.Profile(name=right.name, alignment=[[" ", *row] for row in right.alignment])
    left_profiles = left.alignment
    right_profiles = right.alignment
    left_length = len(left_profiles[0])
    right_length = len(right_profiles[0])

    score_diag = {}
    score_hori = {}
    score_vert = {}

    # Initialize the affine gap scoring matrices
    score_diag[0, 0] = (0, "start")
    score_hori[0, 0] = (-inf, "impossible")
    score_vert[0, 0] = (-inf, "impossible")

    # Top row scores
    for j in range(1, right_length):
        score_diag[0, j] = (-inf, "impossible")
        score_vert[0, j] = (-inf, "impossible")
        score_hori[0, j] = (gap_open + (j-1) * gap_extend, "diag" if j == 1 else "hori")

    # First column scores
    for i in range(1, left_length):
        score_diag[i, 0] = (-inf, "impossible")
        score_hori[i, 0] = (-inf, "impossible")
        score_vert[i, 0] = (gap_open + (i-1) * gap_extend, "diag" if i == 1 else "vert")

    # Interior scores
    # Iterate over the length of profile 1
    for i in range(1, left_length):
        left_column = [profile[i] for profile in left_profiles]

        # Iterate over the length of profile 2
        for j in range(1, right_length):
            right_column = [profile[j] for profile in right_profiles]

            # Current pair is a residue column vs residue column
            column_v_columm_blosum_score = _msa_sigma(left_column, right_column)
            score_diag[i, j] = max(
                (score_diag[i - 1, j - 1][0] + column_v_columm_blosum_score, "diag"),
                (score_hori[i - 1, j - 1][0] + column_v_columm_blosum_score, "hori"),
                (score_vert[i - 1, j - 1][0] + column_v_columm_blosum_score, "vert"),
                key = itemgetter(0)
            )

            # Current pair is a gap column in profile 2
            score_vert[i, j] = max(
                (score_diag[i - 1, j][0] + gap_open, "diag"),
                (score_vert[i - 1, j][0] + gap_extend, "vert"),
                key = itemgetter(0)
            )

            # Current pair is a gap column in profile 1
            score_hori[i, j] = max(
                (score_diag[i, j-1][0] + gap_open, "diag"),
                (score_hori[i, j-1][0] + gap_extend, "hori"),
                key = itemgetter(0)
            )
   
    return _global_msa_backtrack(left, right, score_diag, score_vert, score_hori, left_length-1, right_length-1)


def _walk_tree(
        node:btm.Cluster | btm.Protein,
        backend: str,
        gap_open: int,
        gap_extend: int,
    ) -> btm.Profile:
    '''
    Recursive function to walk through the Cluster tree.
    '''
    if isinstance(node, btm.Protein):
        return btm.Profile(name=[node.name], alignment=[list(node.sequence)])

    left_profile = _walk_tree(node.left, backend, gap_open, gap_extend)
    right_profile = _walk_tree(node.right, backend, gap_open, gap_extend)

    if backend == "rust":
        names, sequences = btrust._merge(
            left_profile.name,
            ["".join(seq) for seq in left_profile.alignment],
            right_profile.name,
            ["".join(seq) for seq in right_profile.alignment],
            gap_open,
            gap_extend,
        )

        return btm.Profile(name=names, alignment=[list(seq) for seq in sequences],
)
    else:
        return _merge(left_profile, right_profile, gap_open, gap_extend)


def _move_reference_sequence_to_top(
        profile:btm.Cluster,
        reference_sequence:btm.Protein
    ) -> btm.Cluster:
    '''
    Moves the sequence that was first in the user input to the top of the alignment.

    There it will be used as the reference sequence for which --diff is checked against.
    '''

    if profile.left == reference_sequence:
        return profile

    if isinstance(profile.left, btm.Protein) or reference_sequence not in profile.left.members:
        profile.swap_nodes()
        _move_reference_sequence_to_top(profile.left, reference_sequence)

    elif isinstance(profile.left, btm.Cluster):
        _move_reference_sequence_to_top(profile.left, reference_sequence)

    return profile


def msa(
        file:str,
        backend:str,
        gap_open:int,
        gap_extend:int,
    ) -> btm.Profile:
    r'''
    Function that takes in sequences in FASTA format, and returns a multiple sequence alignemnt.

    Input to MSA must be a string that represents a FASTA file location, or a string formatted like a FASTA file.

    Command-line usage with a FASTA-formatted string:

    "buttle msa $'>name1\nsequence1\n>name2\nsequence2\n>name3\nsequence3'"
    '''

    if Path(file).exists():
        sequences = btp.parse_fasta_file(file, allow_x = True)
    else:
        sequences = btp.parse_fasta_string(file, allow_x = True)

    reference_sequence = sequences[0]

    pairwise_alignments = _build_pairwise_alignments(sequences, backend)

    original_distance_matrix = _build_original_distance_matrix(sequences, pairwise_alignments)

    tree = _build_guide_tree(sequences, original_distance_matrix)

    arranged_tree = _move_reference_sequence_to_top(tree, reference_sequence)

    return _walk_tree(arranged_tree, backend, gap_open, gap_extend)



def add_consensus(
        profile:btm.Profile
) -> btm.Profile:
    '''
    If the --consensus flag is provided, this function is called.

    It calculates the consensus residue at each position and adds the complete consensus sequence to the bottom of the profile.
    '''
    
    consensus = []

    aligned_length = len(profile.alignment[0])

    for i in range(aligned_length):
        residues = []
        for alignment in profile.alignment:
            residues.append(alignment[i])
        count = Counter(residues)
        consensus.append(count.most_common(1)[0][0])

    profile.name.append("Consensus")
    profile.alignment.append(consensus)

    return profile

