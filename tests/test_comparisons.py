import buttle.comparisons as btcomp
import buttle.models as btm
import pytest

def test_motifs():
    protein1 = btm.Protein("","AFRTcKKYDgEKLREpLDwD")
    protein2 = btm.Protein("", "AFRTgEKLREDwD")
    similars = btcomp.motifs(protein1, protein2)
    expected = {
        (4, 4.0): [
            (('AFRT', 1, 4), ('AFRT', 1, 4), 4), 
            (('GEKL', 5, 8), ('GEKL', 10, 13), 4), 
            (('EKLR', 6, 9), ('EKLR', 11, 14), 4), 
            (('KLRE', 7, 10), ('KLRE', 12, 15), 4)
        ], 
        (5, 5.0): [
            (('GEKLR', 5, 9), ('GEKLR', 10, 14), 5), 
            (('EKLRE', 6, 10), ('EKLRE', 11, 15), 5)
        ], 
        (6, 5.1): [
            (('GEKLRE', 5, 10), ('GEKLRE', 10, 15), 6)
        ], 
        (7, 5.95): [
            (('TGEKLRE', 4, 10), ('DGEKLRE', 9, 15), 6), 
            (('GEKLRED', 5, 11), ('GEKLREP', 10, 16), 6)
        ]
    }
    assert similars == expected


def test_motifs_none_found():
    protein1 = btm.Protein("","APMEVTQ")
    protein2 = btm.Protein("","RSQKTAV")
    similars = btcomp.motifs(protein1, protein2)
    assert similars == {}


def test_motifs_small_example():
    protein1 = btm.Protein("","APMEVTQ")
    protein2 = btm.Protein("","APMEVSV")
    similars = btcomp.motifs(protein1, protein2)
    expected = {
        (4, 4.0): [
            (('APME', 1, 4), ('APME', 1, 4), 4),
            (('PMEV', 2, 5), ('PMEV', 2, 5), 4)
        ],
        (5, 5.0): [
            (('APMEV', 1, 5), ('APMEV', 1, 5), 5)
        ],
        (6, 5.1): [
            (('APMEVT', 1, 6), ('APMEVS', 1, 6), 5.5)
        ]
    }
    assert similars == expected


def test_identity_fasta_file_input(tmp_path):

    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seq1\n"
        "APLPERDFGNMTREEWQ\n"
        ">seq2\n"
        "APIDKDFGQMPTREDFQ\n"
        ">seq3\n"
        "GPVDREYFQMSKDEYN"
    )

    expected = {
        ("seq1", "seq2"): 10/18,
        ("seq1", "seq3"): 4/17,
        ("seq2", "seq3"): 4/17
    }

    output = btcomp.identity(fasta)
    assert output == expected


def test_identity_string_input():

    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN"

    expected = {
        ("seq1", "seq2"): 10/18,
        ("seq1", "seq3"): 4/17,
        ("seq2", "seq3"): 4/17
    }

    output = btcomp.identity(fasta)
    assert output == expected