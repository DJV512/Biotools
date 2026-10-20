from itertools import combinations
from collections import Counter
from pathlib import Path
from math import inf
from copy import copy
from . import models as btm
from . import pairwise_align as btpa
from . import constants as btc
from . import parsers as btp
import biotools_rust as btrust



def _build_pairwise_alignments(
        sequences: list[btm.Protein]
    ) -> list[btm.AlignedProteinPair]:
    '''
    Private function for performing pairwise alignments of all pairs of Proteins in the sequences list.
    '''
    pairwise_alignments = []
    for protein1, protein2 in combinations(sequences, 2):
        align1, align2 = btrust._align_protein(protein1.sequence, protein2.sequence)
        align1 = list(align1)
        align2 = list(align2)
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
                score += btc.BIOTOOLS_COMPARE[left_residue][right_residue]
    return score/comparisons


def _global_msa_backtrack(
        left:btm.Profile, 
        right:btm.Profile, 
        traceback:dict[tuple[int, int], str], 
        i:int, 
        j:int, 
        final_profile:btm.Profile|None = None
    ) -> btm.Profile:
    '''
    Private recursive method for taking a profile-profile score table and backtracking to assemble the optimal alignment.
    '''

    #Initialize the final_profile at the first call to this function
    if final_profile is None:
        total_sequences = len(left.alignment)+len(right.alignment)
        final_profile = btm.Profile(name=left.name + right.name, alignment=[[] for _ in range(total_sequences)])

    #Base case
    if i == 0 or j == 0:

        if i > 0:
            for x in range (1, i+1):
                for z, alignment in enumerate(left.alignment):
                    final_profile.alignment[z].append(alignment[x])
                for z, alignment in enumerate(right.alignment):
                    final_profile.alignment[len(left.alignment) + z].append("-")

        if j > 0:
            for x in range (1, j+1):
                for z, alignment in enumerate(left.alignment):
                    final_profile.alignment[z].append("-")
                for z, alignment in enumerate(right.alignment):
                    final_profile.alignment[len(left.alignment) + z].append(alignment[x])

        return final_profile

    #Recursive calls based on the traceback matrix
    if traceback[i,j] == "diag":
        _global_msa_backtrack(left, right, traceback, i-1, j-1, final_profile)
        for z, alignment in enumerate(left.alignment):
            final_profile.alignment[z].append(alignment[i])
        for z, alignment in enumerate(right.alignment):
            final_profile.alignment[len(left.alignment) + z].append(alignment[j])
        return final_profile

    elif traceback[i,j] == "vert":
        _global_msa_backtrack(left, right, traceback, i-1, j, final_profile)
        for z, alignment in enumerate(left.alignment):
            final_profile.alignment[z].append(alignment[i])
        for z, alignment in enumerate(right.alignment):
            final_profile.alignment[len(left.alignment) + z].append("-")
        return final_profile

    else:
        _global_msa_backtrack(left, right, traceback, i, j-1, final_profile)
        for z, alignment in enumerate(left.alignment):
            final_profile.alignment[z].append("-")
        for z, alignment in enumerate(right.alignment):
            final_profile.alignment[len(left.alignment) + z].append(alignment[j])
        return final_profile




def _merge(
        left: btm.Profile,
        right: btm.Profile
    ) -> tuple[btm.Profile, dict[tuple[int, int], str]]:
    '''
    Private function for merging Proteins and Profiles during multiple alignment.
    '''

    left_profiles = left.alignment
    right_profiles = right.alignment

    for profile in left_profiles:
        profile.insert(0, " ")
    for profile in right_profiles:
        profile.insert(0, " ")

    left_length = len(left_profiles[0])
    right_length = len(right_profiles[0])

    traceback = {}

    #Pre-populate the score matrix with the horizontal first row and vertical first column
    score = {(0, 0): 0}

    left_gap_column = ["-"] * len(left_profiles)
    right_gap_column = ["-"] * len(right_profiles)

    # Top row scores
    for j in range(1, right_length):
        right_column = [profile[j] for profile in right_profiles]

        score[0, j] = (
            score[0, j - 1]
            + _msa_sigma(left_gap_column, right_column)
        )

    # First column scores
    for i in range(1, left_length):
        left_column = [profile[i] for profile in left_profiles]

        score[i, 0] = (
            score[i - 1, 0]
            + _msa_sigma(left_column, right_gap_column)
        )

    # Interior scores
    for i in range(1, left_length):
        left_column = [profile[i] for profile in left_profiles]

        for j in range(1, right_length):
            right_column = [profile[j] for profile in right_profiles]

            vertical_score = score[i - 1, j] + _msa_sigma(left_column, right_gap_column)
            horizontal_score = score[i, j - 1] + _msa_sigma(left_gap_column, right_column)
            diagonal_score = score[i - 1, j - 1] + _msa_sigma(left_column, right_column)
            best_score = max(vertical_score, horizontal_score, diagonal_score)

            score[i,j] = best_score

            if diagonal_score == best_score:
                traceback[i,j] = "diag"
            elif vertical_score == best_score:
                traceback[i,j] = "vert"
            else:
                traceback[i,j] = "horz"
                
    return _global_msa_backtrack(left, right, traceback, left_length-1, right_length-1)


def _walk_tree(
        node:btm.Cluster | btm.Protein
    ) -> btm.Profile:
    '''
    Recursive function to walk through the Cluster tree.
    '''
    if isinstance(node, btm.Protein):
        return btm.Profile(name=[node.name], alignment=[list(node.sequence)])

    left_profile = _walk_tree(node.left)
    right_profile = _walk_tree(node.right)

    return _merge(left_profile, right_profile)


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
        file:str
    ) -> btm.Profile:
    r'''
    Function that takes in sequences in FASTA format, and returns a multiple sequence alignemnt.

    Input to MSA must be a string that represents a FASTA file location, or a string formatted like a FASTA file.

    Command-line usage with a FASTA-formatted string:

    "biotools msa $'>name1\nsequence1\n>name2\nsequence2\n>name3\nsequence3'"
    '''

    if Path(file).exists():
        sequences = btp.parse_fasta_file(file, allow_x = True)
    else:
        sequences = btp.parse_fasta_string(file, allow_x = True)

    reference_sequence = sequences[0]
    
    pairwise_alignments = _build_pairwise_alignments(sequences)

    original_distance_matrix = _build_original_distance_matrix(sequences, pairwise_alignments)

    tree = _build_guide_tree(sequences, original_distance_matrix)

    arranged_tree = _move_reference_sequence_to_top(tree, reference_sequence)

    return _walk_tree(arranged_tree)



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

