import buttle.models as btm
import pytest
from hypothesis import given, strategies as st

def test_charge():
    protein = btm.Protein("", "ArFTyGPl")
    charge = protein.charge()
    assert charge == 1


def test_charge_no_charged_residues():
    protein = btm.Protein("","AFTYGPL")
    charge = protein.charge()
    assert charge == 0


def test_charge_many_charged_residues():
    protein = btm.Protein("","AFRTKKYDGEKLREPLDWD")
    charge = protein.charge()
    assert charge == 0


def test_protein_rejects_invalid_aas():
    with pytest.raises(ValueError):
        protein = btm.Protein("","ASWPMCLRUFCQW")


@given(st.text(alphabet="ACDEFGHIKLMNPQRSTVWY",
               min_size = 1,
               max_size = 100))
def test_charge_hypothesis(protein):
    test_protein = btm.Protein("", protein)
    expected = test_protein.sequence.count("K") + test_protein.sequence.count("R") - test_protein.sequence.count("D") - test_protein.sequence.count("E")
    assert test_protein.charge() == expected


def test_mw():
    protein = btm.Protein("","AFRTKKYDGEKLREPLDWD")
    mw = protein.mw()
    assert mw == pytest.approx(2367.6321)


def test_ec():
    protein = btm.Protein("","AFRTcKKYDgEKLREpLDwD")
    ec = protein.ec()
    assert ec == 7090


def test_glycosylation():
    protein = btm.Protein("","YQSFRPMNLSQ")
    sites = protein.glycosylation()
    assert sites == [7]


def test_glycosylation_no_sites():
    protein = btm.Protein("","YQNSFRPMNL")
    sites = protein.glycosylation()
    assert sites == []


def test_glycosylation_cant_have_pro_as_x():
    protein = btm.Protein("","YQSFRPMNPS")
    sites = protein.glycosylation()
    assert sites == []


def test_glycosylation_overlapping():
    protein = btm.Protein("","YQSFRPMNNST")
    sites = protein.glycosylation()
    assert sites == [7, 8]
