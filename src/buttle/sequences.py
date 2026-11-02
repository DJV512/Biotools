from . import constants as btc
from . import decorators as btd
lazy from collections import Counter

@btd.validate_dna()
def transcribe(dna:str) -> str:
    '''
    Given a DNA sequence, returns the RNA transcription.\n
    Required argument: dna
    '''

    #Transcribe the DNA
    rna = ""
    for base in dna:
        base = base.upper()
        if base == "A":
            rna += "A"
        elif base == "C":
            rna += "C"
        elif base == "G":
            rna += "G"
        else:
            rna += "U"
    return rna


@btd.validate_dna(allow_n=True)
def translate(dna:str) -> tuple[str, str, str]:
    '''
    Given a DNA sequence, returns the protein translations in all three forward frames.\n
    Required argument: dna
    '''

    # Initialize the three translations with empty strings.
    proteinf1 = ""
    proteinf2 = ""
    proteinf3 = ""

    # Frame 1
    codons = [dna[i:i+3] for i in range(0, len(dna), 3)]
    for codon in codons:
        if len(codon) == 3:
            try:
                proteinf1 += btc.CODON_TABLE[codon.upper()]
            except KeyError:
                proteinf1 += "X"

    # Frame 2
    codons = [dna[i:i+3] for i in range(1, len(dna), 3)]
    for codon in codons:
        if len(codon) == 3:
            try:
                proteinf2 += btc.CODON_TABLE[codon.upper()]
            except KeyError:
                proteinf2 += "X"

    # Frame 3
    codons = [dna[i:i+3] for i in range(2, len(dna), 3)]
    for codon in codons:
        if len(codon) == 3:
            try:
                proteinf3 += btc.CODON_TABLE[codon.upper()]
            except KeyError:
                proteinf3 += "X"

    return proteinf1, proteinf2, proteinf3


@btd.validate_dna()
def reverse_comp(dna:str) -> str:
    '''
    Given a DNA sequence, returns the reverse complement DNA sequence.\n
    Required argument: dna
    '''

    # Reverse transcribe DNA
    rev_dna = ""
    dna = dna[::-1].upper()
    for base in dna:
        if base == "A":
            rev_dna += "T"
        elif base == "C":
            rev_dna += "G"
        elif base == "G":
            rev_dna += "C"
        else:
            rev_dna += "A"
    return rev_dna


@btd.validate_dna()
def gc(dna:str) -> tuple[int, Counter, float, float]:
    '''
    Given a DNA sequence, returns information about the base contents.
    '''
    dna = dna.upper()
    length = len(dna)
    contents = Counter(dna)
    at_content = 100 * (contents["A"] + contents["T"]) / length
    gc_content = 100 * (contents["G"] + contents["C"]) / length
    return length, contents, at_content, gc_content


def length(sequence:str) -> int:
    '''
    CLI wrapper for built in len function.\n
    Required arguments: sequence
    '''

    return len(sequence.strip())


def find_pos(sequence:str, search_string:str) -> list:
    '''
    Given a sequence and a search string, this method will find and return the position of the search string in the sequence, if it exists.\n
    Required arguments: sequence, search_string
    '''

    search_length = len(search_string)
    sequence_length = len(sequence)
    hits = []
    for i in range(0, sequence_length-search_length+1):
        sequence_part = sequence[i:i+search_length]
        if sequence_part.lower() == search_string.lower():
            hits.append(i+1)
    return hits


