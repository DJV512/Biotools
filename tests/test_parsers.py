import biotools.parsers as btp
import biotools.models as btm
import pytest

def test_fasta_reader_minimal_fasta(tmp_path):
    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seq1\n"
        "APLPERDFGNMTREEWQ\n"
        ">seq2\n"
        "APIDKDFGQMPTREDFQ\n"
        ">seq3\n"
        "GPVDREYFQMSKDEYN\n"
        ">seq4\n"
        "WRTDPYMNCQWSLKCF\n"
        ">seq5\n"
        "WRSDPFMNCQWTIRCY\n"
    )

    expected = [
        btm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]

    output = btp.parse_fasta_file(fasta)
    assert output == expected


def test_fasta_reader_blank_lines(tmp_path):
    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seq1\n"
        "APLPERDFGNMTREEWQ\n"
        "\n"
        ">seq2\n"
        "APIDKDFGQMPTREDFQ\n"
        "\n"
        "\n"
        ">seq3\n"
        "GPVDREYFQMSKDEYN\n"
        "\n"
        ">seq4\n"
        "WRTDPYMNCQWSLKCF\n"
        "\n"
        ">seq5\n"
        "WRSDPFMNCQWTIRCY\n"
    )

    expected = [
        btm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]
    
    output = btp.parse_fasta_file(fasta)
    assert output == expected


def test_fasta_reader_multiline(tmp_path):
    fasta = tmp_path / "test.fasta"

    fasta.write_text(
        ">seq1\n"
        "APLPERDFGN\n"
        "MTREEWQ\n"
        "\n"
        ">seq2\n"
        "APIDK\n"
        "DFGQMP\n"
        "TREDFQ\n"
        "\n"
        ">seq3\n"
        "G\n"
        "PVDREYFQMSKD\n"
        "EYN\n"
        "\n"
        ">seq4\n"
        "WRTDPYMNCQWSLKCF\n"
        "\n"
        ">seq5\n"
        "WRS\n"
        "DPFMNCQWTIRCY\n"
    )

    expected = [
        btm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]

    output = btp.parse_fasta_file(fasta)
    assert output == expected


def test_fasta_reader_stress_test(tmp_path):
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

    expected = [
        btm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]
    
    output = btp.parse_fasta_file(fasta)
    assert output == expected

def test_fasta_reader_one_sequence_rejected(tmp_path):
    fasta = tmp_path / "test.fasta"
    
    fasta.write_text(
        ">seq1\n"
        "APLPERDFGN\n"
        "MTREEWQ\n"
        "\n"
    )

    with pytest.raises(ValueError):
        output = btp.parse_fasta_file(fasta)


def test_fasta_reader_duplicate_names_rejected(tmp_path):
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
        ">seq3\n"
        "WRS\n"
        "DPFMNCQWTIRCY\n"
        "\n"
        "\n"
    )

    with pytest.raises(ValueError):
        output = btp.parse_fasta_file(fasta)


def test_parse_fasta_string_works_with_proper_input():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"

    expected = [
        btm.Protein(name='seq1', sequence='APLPERDFGNMTREEWQ'), 
        btm.Protein(name='seq2', sequence='APIDKDFGQMPTREDFQ'), 
        btm.Protein(name='seq3', sequence='GPVDREYFQMSKDEYN'), 
        btm.Protein(name='seq4', sequence='WRTDPYMNCQWSLKCF'), 
        btm.Protein(name='seq5', sequence='WRSDPFMNCQWTIRCY')
    ]

    output = btp.parse_fasta_string(fasta)
    assert output == expected


def test_parse_sequences_rejects_regular_string():
    fasta = "WRSDPFMNCQWTIRCY"

    with pytest.raises(ValueError):
        output = btp.parse_fasta_string(fasta)


def test_parse_sequences_rejects_missing_carrot():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\nseq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"

    with pytest.raises(ValueError):
        output = btp.parse_fasta_string(fasta)


def test_parse_sequences_rejects_missing_sequence():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq5\nWRSDPFMNCQWTIRCY\n"

    with pytest.raises(ValueError):
        output = btp.parse_fasta_string(fasta)


def test_parse_sequences_rejects_single_sequence():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ"

    with pytest.raises(ValueError):
        output = btp.parse_fasta_string(fasta)


def test_parse_sequences_rejects_duplicate_names():
    fasta = ">seq1\nAPLPERDFGNMTREEWQ\n>seq2\nAPIDKDFGQMPTREDFQ\n>seq3\nGPVDREYFQMSKDEYN\n>seq4\nWRTDPYMNCQWSLKCF\n>seq2\nWRSDPFMNCQWTIRCY\n"

    with pytest.raises(ValueError):
        output = btp.parse_fasta_string(fasta)

def test_parse_fasta_string_empty_input():
    with pytest.raises(ValueError, match="No sequences found"):
        btp.parse_fasta_string("")