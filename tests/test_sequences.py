import buttle.sequences as bts
import pytest
from collections import Counter

def test_transcribe():
    dna = "actgctcctagctatactca"
    transcription = bts.transcribe(dna)
    assert transcription == "ACUGCUCCUAGCUAUACUCA"


def test_transcribe_rejects_invalid_dna():
    dna = "acctagctagucgtagtcagtc"
    with pytest.raises(ValueError):
        bts.transcribe(dna)


def test_translate():
    dna = "gtatcatgtctacgtgcgcggactgtcatag"
    t1, t2, t3 = bts.translate(dna)
    assert t1 == "VSCLRARTVI"
    assert t2 == "YHVYVRGLS*"
    assert t3 == "IMSTCADCH"


def test_translate_rejects_invalid_dna():
    dna = "acctctagucgt"
    with pytest.raises(ValueError):
        bts.translate(dna)


def test_translate_accepts_n():
    dna = "acctagctagncgtagntcagtnc"
    t1, t2, t3 = bts.translate(dna)
    assert t1 == "T*LXVXQX"
    assert t2 == "PS*X*XS"
    assert t3 == "LAXRXSX"


def test_reverse_comp():
    dna = "gtatcatgtctacgtgcgcggactgtcatag"
    rcomp = bts.reverse_comp(dna)
    assert rcomp == "CTATGACAGTCCGCGCACGTAGACATGATAC"


def test_reverse_comp_rejects_invalid_dna():
    dna = "acctagctagucgtagtcagtc"
    with pytest.raises(ValueError):
        bts.reverse_comp(dna)


def test_gc():
    dna = "gtatcatgtctacgtgcgcggactgtcatag"
    length, contents, at_content, gc_content = bts.gc(dna)
    assert length == 31
    assert contents == Counter({'G': 9, 'T': 9, 'C': 7, 'A': 6})
    assert at_content == pytest.approx(48.3871)
    assert gc_content == pytest.approx(51.6129)


def test_length():
    sequence = "adkfhaksdhfkasf"
    length = bts.length(sequence)
    assert length == 15


def test_length_strip():
    sequence = "  adkfhaksdhfkasf "
    length = bts.length(sequence)
    assert length == 15


def test_find_pos():
    sequence = "asdfasdnaspjnfasdpknfasksjnfs"
    motif = "FAS"
    hits = bts.find_pos(sequence, motif)
    assert hits == [4, 14, 21]


def test_find_pos_overlapping():
    sequence = "asdfasdfanaspjnfasdpknfasksjnfs"
    motif = "AsDfA"
    hits = bts.find_pos(sequence, motif)
    assert hits == [1, 5]


def test_find_pos_motif_missing():
    sequence = "slkfjasngoigwrvnosanfalsd"
    motif = "dsafsdf"
    hits = bts.find_pos(sequence, motif)
    assert hits == []