import biotools.pairwise_align as btpa
import biotools.models as btm
import biotools.parsers as btparse
import biotools_rust as btrust
import pytest
from hypothesis import given, Verbosity, settings, strategies as st
from itertools import combinations
from pathlib import Path


def test_sigma_dna_match():
    base1 = "A"
    base2 = "A"
    score = btpa._sigma_dna(base1, base2)
    assert score == 8


def test_sigma_dna_mismatch():
    base1 = "A"
    base2 = "G"
    score = btpa._sigma_dna(base1, base2)
    assert score == -5


def test_pairwise_dna_case_insensitive():
    dna1 = "acagtgca"
    dna2 = "actgtata"
    align1, align2 = btpa.pairwise_dna(dna1, dna2)
    assert (align1, align2) == (["A", "C", "A", "G", "T", "G", "C", "A"],["A", "C", "T", "G", "T", "A", "T", "A"])


def test_pairwise_dna_leading_gap_first():
    dna1 = "AGTGTC"
    dna2 = "GTGATC"
    align1, align2 = btpa.pairwise_dna(dna1, dna2)
    assert (align1, align2) == (["A", "G", "T", "G", "-", "T", "C"],["-", "G", "T", "G", "A", "T", "C"])


def test_pairwise_dna_leading_gap_second():
    dna2 = "AGTGTC"
    dna1 = "GTGATC"
    align1, align2 = btpa.pairwise_dna(dna1, dna2)
    assert (align1, align2) == (["-", "G", "T", "G", "A", "T", "C"], ["A", "G", "T", "G", "-", "T", "C"])


def test_pairwise_dna_trailing_gap():
    dna1 = "AGTGTCA"
    dna2 = "AGTGATC"
    align1, align2 = btpa.pairwise_dna(dna1, dna2)
    assert (align1, align2) == (["A", "G", "T", "G", "-", "T", "C", "A"],["A", "G", "T", "G", "A", "T", "C", "-"])


def test_pairwise_dna_identical():
    dna1 = "AGTGTC"
    dna2 = "AGTGTC"
    align1, align2 = btpa.pairwise_dna(dna1, dna2)
    assert (align1, align2) == (["A", "G", "T", "G", "T", "C"],["A", "G", "T", "G", "T", "C"])


def test_pairwise_dna_single_mismatch():
    dna1 = "AGTGTC"
    dna2 = "AGAGTC"
    align1, align2 = btpa.pairwise_dna(dna1, dna2)
    assert (align1, align2) == (["A", "G", "T", "G", "T", "C"],["A", "G", "A", "G", "T", "C"])


def test_pairwise_dna_bad_dna1():
    dna1 = "AGTGTR"
    dna2 = "AGAGTC"
    with pytest.raises(ValueError):
        btpa.pairwise_dna(dna1, dna2)


def test_pairwise_dna_bad_dna2():
    dna1 = "AGTGTC"
    dna2 = "AGAGTR"
    with pytest.raises(ValueError):
        btpa.pairwise_dna(dna1, dna2)


def test_pairwise_dna_keyword_arguments():
    dna1 = "AGTGTC"
    dna2 = "AGTGTC"
    align1, align2 = btpa.pairwise_dna(dna1=dna1, dna2=dna2)
    assert (align1, align2) == (["A", "G", "T", "G", "T", "C"],["A", "G", "T", "G", "T", "C"])


def test_pairwise_dna_keyword_arguments_dna2_bad():
    dna1 = "AGTGTC"
    dna2 = "AGMGTC"
    with pytest.raises(ValueError):
        btpa.pairwise_dna(dna1=dna1, dna2=dna2)


def test_pairwise_dna_doesnt_accept_n():
    dna1 = "AGTGTC"
    dna2 = "AGNGTC"
    with pytest.raises(ValueError):
        btpa.pairwise_dna(dna1=dna1, dna2=dna2)


def test_pairwise_protein_case_insesitive():
    protein1 = btm.Protein("","avmyerf")
    protein2 = btm.Protein("","llmfdkf")
    align1, align2 = btpa.pairwise_protein(protein1, protein2)
    assert (align1, align2) == (["A", "V", "M", "Y", "E", "R", "F"],["L", "L", "M", "F", "D", "K", "F"])


def test_pairwise_protein_leading_gap_first():
    protein1 = btm.Protein("","AVMYERF")
    protein2 = btm.Protein("","QAVMFDKF")
    align1, align2 = btpa.pairwise_protein(protein1, protein2)
    assert (align1, align2) == (["-","A", "V", "M", "Y", "E", "R", "F"],["Q", "A", "V", "M", "F", "D", "K", "F"])


def test_pairwise_protein_leading_gap_second():
    protein1 = btm.Protein("","QAVMYERF")
    protein2 = btm.Protein("","AVMFDKF")
    align1, align2 = btpa.pairwise_protein(protein1, protein2)
    assert (align1, align2) == (["Q", "A", "V", "M", "Y", "E", "R", "F"],["-", "A", "V", "M", "F", "D", "K", "F"])


def test_pairwise_protein_single_point_mutation():
    protein1 = btm.Protein("","AVMYERF")
    protein2 = btm.Protein("","AVMYGRF")
    align1, align2 = btpa.pairwise_protein(protein1, protein2)
    assert (align1, align2) == (["A", "V", "M", "Y", "E", "R", "F"],["A", "V", "M", "Y", "G", "R", "F"])


def test_pairwise_protein_both_gap_types():
    protein1 = btm.Protein("","AQVMYERF")
    protein2 = btm.Protein("","AVMYERPF")
    align1, align2 = btpa.pairwise_protein(protein1, protein2)
    assert (align1, align2) == (["A", "Q", "V", "M", "Y", "E", "R", "-", "F"],["A", "-", "V", "M", "Y", "E", "R", "P", "F"])

protein_strategy = st.text(alphabet="ACDEFGHIKLMNPQRSTVWY",
               min_size = 15,
               max_size = 100)

@given(protein_strategy, protein_strategy)
# @settings(max_examples=500)
def test_pairwise_align_hypothesis(protein1, protein2):
    test_protein1 = btm.Protein("", protein1)
    test_protein2 = btm.Protein("", protein2)
    align1, align2 = btpa.pairwise_protein(test_protein1, test_protein2)
    assert len(align1) == len(align2)
    assert "".join(align1).replace("-", "") == test_protein1.sequence
    assert "".join(align2).replace("-", "") == test_protein2.sequence


def test_pairwise_affine_gap_implementation():
    protein1 = btm.Protein("", "AAQPSFCGKLPRMIHLQWTEVFYFLMES")
    protein2 = btm.Protein("", "AAQAFPSFCGRMIHLQVFYLMLIQKMES")

    align1, align2 = btpa.pairwise_protein(protein1, protein2)

    assert "".join(align1).replace("-", "") == protein1.sequence
    assert "".join(align2).replace("-", "") == protein2.sequence
    assert (align1, align2) == (
        ["A", "A", "Q", "-", "-", "P", "S", "F", "C", "G", "K", "L", "P", "R", "M", "I", "H", "L", "Q", "W", "T", "E", "V", "F", "Y", "F", "L", "-", "-", "-", "-", "M", "E", "S"], 
        ["A", "A", "Q", "A", "F", "P", "S", "F", "C", "G", "-", "-", "-", "R", "M", "I", "H", "L", "Q", "-", "-", "-", "V", "F", "Y", "L", "M", "L", "I", "Q", "K", "M", "E", "S"]
        )


def test_pairwise_rust_implementation():
    protein1 = "AAQPSFCGKLPRMIHLQWTEVFYFLMES"
    protein2 = "AAQAFPSFCGRMIHLQVFYLMLIQKMES"

    align1, align2 = btrust._align_protein(protein1, protein2)
    assert align1.replace("-", "") == protein1
    assert align2.replace("-", "") == protein2
    assert (align1, align2) == (
        "AAQ--PSFCGKLPRMIHLQWTEVFYFL----MES",
        "AAQAFPSFCG---RMIHLQ---VFYLMLIQKMES",
    )

def test_rust_matches_python_affine_gap():
    protein1 = btm.Protein("Protein 1", "AAQPSFCGKLPRMIHLQWTEVFYFLMES")
    protein2 = btm.Protein("Protein 2", "AAQAFPSFCGRMIHLQVFYLMLIQKMES")

    python1, python2 = btpa.pairwise_protein(protein1, protein2)
    rust1, rust2 = btrust._align_protein(protein1.sequence, protein2.sequence)

    assert (rust1, rust2) == ("".join(python1), "".join(python2))

proteins = btparse.parse_fasta_file(Path(__file__).parent / "data" / "synthetic_proteins.fasta", allow_x = True)


@pytest.mark.parametrize("protein1, protein2", list(combinations(proteins, 2)))
def test_rust_matches_python(protein1, protein2):
    python1, python2 = btpa.pairwise_protein(protein1, protein2)

    rust1, rust2 = btrust._align_protein(
        protein1.sequence,
        protein2.sequence,
    )

    assert rust1 == "".join(python1)
    assert rust2 == "".join(python2)



AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
protein_sequences = st.text(
    alphabet=AMINO_ACIDS,
    min_size=1,
    max_size=100,
)
@given(seq1=protein_sequences, seq2=protein_sequences)
@settings(max_examples=100, deadline=None)
def test_rust_matches_python_hypothesis(seq1, seq2):
    protein1 = btm.Protein("protein1", seq1)
    protein2 = btm.Protein("protein2", seq2)

    python1, python2 = btpa.pairwise_protein(protein1, protein2)
    rust1, rust2 = btrust._align_protein(seq1, seq2)

    assert rust1 == "".join(python1)
    assert rust2 == "".join(python2)
    assert len(rust1) == len(rust2)
    assert rust1.replace("-", "") == seq1
    assert rust2.replace("-", "") == seq2