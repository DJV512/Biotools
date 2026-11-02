from hypothesis import given, Verbosity, settings, strategies as st
import buttle.msa as btmsa
import buttle.models as btm
import pytest

protein_strategy = st.text(alphabet="ACDEFGHIKLMNPQRSTVWY",
               min_size = 80,
               max_size = 120)

proteins_strategy = st.lists(protein_strategy, min_size = 10, max_size = 10)

@given(proteins_strategy)
@settings(verbosity=Verbosity.verbose, max_examples=10, deadline=None)
def test_msa_hypothesis(proteins):

    originals = {}
    sequences = ""
    for n, protein in enumerate(proteins):
        sequences += f">protein{n}\n{protein}\n"
        originals[f"protein{n}"] = protein

    profile = btmsa.msa(sequences, "python", -10, -3)

    lengths = set()
    for name, align in zip(profile.name, profile.alignment, strict=True):
        lengths.add(len(align))
        assert originals[name] == "".join(align).replace("-","")
    assert len(lengths) == 1



def test_msa_correctly_runs_with_string_fasta_input():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"
    profile = btmsa.msa(fasta, "python", -10, -3)

    originals = {
        "seq1": "APLPERDFGNMTREEWQ",
        "seq2": "APIDKDFGQMPTREDFQ",
        "seq3": "GPVDREYFQMSKDEYN",
        "seq4": "WRTDPYMNCQWSLKCF",
        "seq5": "WRSDPFMNCQWTIRCY"
    }

    lengths = set()
    for name, align in zip(profile.name, profile.alignment, strict=True):
        lengths.add(len(align))
        assert originals[name] == "".join(align).replace("-","")
    assert len(lengths) == 1


def test_msa_correctly_runs_with_file_fasta_input(tmp_path):
    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seq1\n"
        "APLPERDFGN\n"
        "MTREEWQ\n"
        "\n"
        ">seq2\n"
        "APIDK\n"
        "DFGQMPTREDFQ\n"
        "\n"
        "\n"
        ">seq3\n"
        "GPVDREYFQMSKD\n"
        "EYN\n"
        "\n"
        ">seq4\n"
        "WRTDPY\n"
        "MNCQWSLKCF\n"
        "\n"
        "\n"
        "\n"
        ">seq5\n"
        "WRS\n"
        "DPFMNCQWTIRCY\n"
        "\n"
        "\n"
    )

    profile = btmsa.msa(fasta, "python", -10, -3)

    originals = {
        "seq1": "APLPERDFGNMTREEWQ",
        "seq2": "APIDKDFGQMPTREDFQ",
        "seq3": "GPVDREYFQMSKDEYN",
        "seq4": "WRTDPYMNCQWSLKCF",
        "seq5": "WRSDPFMNCQWTIRCY"
    }

    lengths = set()
    for name, align in zip(profile.name, profile.alignment, strict=True):
        lengths.add(len(align))
        assert originals[name] == "".join(align).replace("-","")
    assert len(lengths) == 1


def test_add_consensus():
    profile = btm.Profile(
        name=['seq1', 'seq2', 'seq3'], 
        alignment=[
            ['A', 'P', 'L', 'P', 'E', 'R', 'D', 'F', 'G', 'N', 'M', '-', 'T', 'R', 'E', 'E', 'W', 'Q'], 
            ['A', 'P', 'I', '-', 'D', 'K', 'D', 'F', 'G', 'Q', 'M', 'P', 'T', 'R', 'E', 'D', 'F', 'Q'], 
            ['G', 'P', 'V', '-', 'D', 'R', 'E', 'Y', 'F', 'Q', 'M', '-', 'S', 'K', 'D', 'E', 'Y', 'N']
        ]
    )

    expected = btm.Profile(
        name=['seq1', 'seq2', 'seq3', 'Consensus'], 
        alignment=[
            ['A', 'P', 'L', 'P', 'E', 'R', 'D', 'F', 'G', 'N', 'M', '-', 'T', 'R', 'E', 'E', 'W', 'Q'], 
            ['A', 'P', 'I', '-', 'D', 'K', 'D', 'F', 'G', 'Q', 'M', 'P', 'T', 'R', 'E', 'D', 'F', 'Q'], 
            ['G', 'P', 'V', '-', 'D', 'R', 'E', 'Y', 'F', 'Q', 'M', '-', 'S', 'K', 'D', 'E', 'Y', 'N'],
            ['A', 'P', 'L', '-', 'D', 'R', 'D', 'F', 'G', 'Q', 'M', '-', 'T', 'R', 'E', 'E', 'W', 'Q'], 
        ]
    )

    output = btmsa.add_consensus(profile)
    assert output == expected


def test_msa_rust_vs_python_static_example():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"
    python_profile = btmsa.msa(fasta, "python", -10, -3)
    rust_profile = btmsa.msa(fasta, "rust", -10, -3)

    originals = {
        "seq1": "APLPERDFGNMTREEWQ",
        "seq2": "APIDKDFGQMPTREDFQ",
        "seq3": "GPVDREYFQMSKDEYN",
        "seq4": "WRTDPYMNCQWSLKCF",
        "seq5": "WRSDPFMNCQWTIRCY"
    }

    python_lengths = set()
    for name, align in zip(python_profile.name, python_profile.alignment, strict=True):
        python_lengths.add(len(align))
        assert originals[name] == "".join(align).replace("-","")
    assert len(python_lengths) == 1

    rust_lengths = set()
    for name, align in zip(rust_profile.name, rust_profile.alignment, strict=True):
        rust_lengths.add(len(align))
        assert originals[name] == "".join(align).replace("-","")
    assert len(rust_lengths) == 1

    assert python_profile == rust_profile

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
protein_sequences = st.text(
    alphabet=AMINO_ACIDS,
    min_size=1,
    max_size=100,
)
@given(seq1=protein_sequences, seq2=protein_sequences, seq3=protein_sequences, seq4=protein_sequences, seq5=protein_sequences, seq6=protein_sequences)
@settings(max_examples=500, deadline=None)
def test_msa_rust_vs_python_hypothesis(seq1, seq2, seq3, seq4, seq5, seq6):
    fasta = f">seq1\n{seq1}\n>seq2\n{seq2}\n>seq3\n{seq3}\n>seq4\n{seq4}\n>seq5\n{seq5}\n>seq6\n{seq6}"

    python_profile = btmsa.msa(fasta, "python", -10, -3)
    rust_profile = btmsa.msa(fasta, "rust", -10, -3)

    assert python_profile == rust_profile

