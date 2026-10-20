from copy import deepcopy

# Styles for terminal sequence output
RESET           = "\033[0m"                     # Resets all styles
BLACK_TERM      = "\033[1;38;2;0;0;0m"          # Black text
RED_TERM        = "\033[1;38;2;255;0;0m"        # Red text
BLUE_TERM       = "\033[1;38;2;0;0;255m"        # Blue text
LIGHT_BLUE_TERM = "\033[1;38;2;0;0;200m"        # Light Blue text
GREEN_TERM      = "\033[1;38;2;0;200;0m"        # Green text
BROWN_TERM      = "\033[1;38;2;160;130;0m"      # Brown text
PURPLE_TERM     = "\033[1;38;2;128;0;160m"      # Purple text
WHITE_TERM      = "\033[1;38;2;255;255;255m"    # Brown text
BG_WHITE        = "\033[48;2;255;255;255m"      # White backround
BOLD            = "\033[1m"                     # Bold text

# Styles for terminal help output
BLUE_HELP       = "\033[1;38;2;59;142;234m"
GREEN_HELP      = "\033[1;38;2;35;209;139m"
CYAN_HELP       = "\033[1;38;2;40;184;219m"
MAGENTA_HELP    = "\033[1;38;2;214;112;214m"
YELLOW_HELP     = "\033[1;38;2;245;245;67m"

# DNA Alignment Terminal Colors
TERMINAL_COLORS_DNA = {
        "match": GREEN_TERM,
        "gap": BROWN_TERM,
        "mismatch": RED_TERM,
    }

# Protein Alignment Terminal Colors
TERMINAL_COLORS_PROTEIN = {
    "match": GREEN_TERM,
    "cons": PURPLE_TERM,
    "non-cons": RED_TERM,
    "gap": BLACK_TERM,
    "white": WHITE_TERM
}

# Styles for image output
BLACK      = (0,0,0)
RED        = (255,0,0)
GREEN      = (0,200,0)
BROWN      = (160,130,0)
WHITE      = (255,255,255)
PURPLE     = (128, 0, 160)
TEAL       = (0, 130, 140)

PADDING = 50

# DNA Alignment Image Colors
DNA_IMAGE_COLORS = {
        "match": GREEN,
        "gap": BROWN,
        "mismatch": RED,
    }

# Protein Alignment Image Colors
PROTEIN_IMAGE_COLORS = {
    "match": GREEN,
    "cons": PURPLE,
    "non-cons": RED,
    "gap": BLACK,
    "white": WHITE
}

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Dictionary that contains key:value pairs of codon:amino_acid for protien translation.
CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

# Dictionary that stores the molecular weight of each amino acid.
MW = {
    "A":71.0788,
    "C":103.1448,
    "D":115.0886,
    "E":129.1155,
    "F":147.1766,
    "G":57.0520,
    "H":137.1412,
    "I":113.1595,
    "K":128.1742,
    "L":113.1595,
    "M":131.1986,
    "N":114.1039,
    "P":97.1167,
    "Q":128.1308,
    "R":156.1876,
    "S":87.0782,
    "T":101.1051,
    "V":99.1326,
    "W":186.2133,
    "Y":163.1760,
}

# Dictionary where each key is an amino acid, and its value is a list of amino acids that would be considered conservative substitutions.
# For motifs to use when finding regions of similarity
CONSERVATIVES = {
        "A":[],
        "C":[],
        "D":["E"],
        "E":["D"],
        "F":["Y", "W"],
        "G":[],
        "H":[],
        "I":["L", "M", "V"],
        "K":["R"],
        "L":["I", "M", "V"],
        "M":["I", "L", "V"],
        "N":["Q"],
        "P":[],
        "Q":["N"],
        "R":["K"],
        "S":["T"],
        "T":["S"],
        "V":["I", "L", "M"],
        "W":["Y", "F"],
        "X":[],
        "Y":["F", "W"],
    }

# Tuple of tuples, where each internal tuple has a lower bound, upper bound, and threshold for a match. Used for motif matching.
THRESHOLD = (
    (4, 5, 1.00),
    (6, 8, 0.85),
    (9, 12, 0.80),
    (13, 20, 0.75),
    (21, 40, 0.70),
)

BLOSUM62 = {
    "A": {"A": 4,  "R": -1, "N": -2, "D": -2, "C": 0,  "Q": -1, "E": -1, "G": 0,  "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1,  "T": 0,  "W": -3, "Y": -2, "V": 0},
    "R": {"A": -1, "R": 5,  "N": 0,  "D": -2, "C": -3, "Q": 1,  "E": 0,  "G": -2, "H": 0,  "I": -3, "L": -2, "K": 2,  "M": -1, "F": -3, "P": -2, "S": -1, "T": -1, "W": -3, "Y": -2, "V": -3},
    "N": {"A": -2, "R": 0,  "N": 6,  "D": 1,  "C": -3, "Q": 0,  "E": 0,  "G": 0,  "H": 1,  "I": -3, "L": -3, "K": 0,  "M": -2, "F": -3, "P": -2, "S": 1,  "T": 0,  "W": -4, "Y": -2, "V": -3},
    "D": {"A": -2, "R": -2, "N": 1,  "D": 6,  "C": -3, "Q": 0,  "E": 2,  "G": -1, "H": -1, "I": -3, "L": -4, "K": -1, "M": -3, "F": -3, "P": -1, "S": 0,  "T": -1, "W": -4, "Y": -3, "V": -3},
    "C": {"A": 0,  "R": -3, "N": -3, "D": -3, "C": 9,  "Q": -3, "E": -4, "G": -3, "H": -3, "I": -1, "L": -1, "K": -3, "M": -1, "F": -2, "P": -3, "S": -1, "T": -1, "W": -2, "Y": -2, "V": -1},
    "Q": {"A": -1, "R": 1,  "N": 0,  "D": 0,  "C": -3, "Q": 5,  "E": 2,  "G": -2, "H": 0,  "I": -3, "L": -2, "K": 1,  "M": 0,  "F": -3, "P": -1, "S": 0,  "T": -1, "W": -2, "Y": -1, "V": -2},
    "E": {"A": -1, "R": 0,  "N": 0,  "D": 2,  "C": -4, "Q": 2,  "E": 5,  "G": -2, "H": 0,  "I": -3, "L": -3, "K": 1,  "M": -2, "F": -3, "P": -1, "S": 0,  "T": -1, "W": -3, "Y": -2, "V": -2},
    "G": {"A": 0,  "R": -2, "N": 0,  "D": -1, "C": -3, "Q": -2, "E": -2, "G": 6,  "H": -2, "I": -4, "L": -4, "K": -2, "M": -3, "F": -3, "P": -2, "S": 0,  "T": -2, "W": -2, "Y": -3, "V": -3},
    "H": {"A": -2, "R": 0,  "N": 1,  "D": -1, "C": -3, "Q": 0,  "E": 0,  "G": -2, "H": 8,  "I": -3, "L": -3, "K": -1, "M": -2, "F": -1, "P": -2, "S": -1, "T": -2, "W": -2, "Y": 2,  "V": -3},
    "I": {"A": -1, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -3, "E": -3, "G": -4, "H": -3, "I": 4,  "L": 2,  "K": -3, "M": 1,  "F": 0,  "P": -3, "S": -2, "T": -1, "W": -3, "Y": -1, "V": 3},
    "L": {"A": -1, "R": -2, "N": -3, "D": -4, "C": -1, "Q": -2, "E": -3, "G": -4, "H": -3, "I": 2,  "L": 4,  "K": -2, "M": 2,  "F": 0,  "P": -3, "S": -2, "T": -1, "W": -2, "Y": -1, "V": 1},
    "K": {"A": -1, "R": 2,  "N": 0,  "D": -1, "C": -3, "Q": 1,  "E": 1,  "G": -2, "H": -1, "I": -3, "L": -2, "K": 5,  "M": -1, "F": -3, "P": -1, "S": 0,  "T": -1, "W": -3, "Y": -2, "V": -2},
    "M": {"A": -1, "R": -1, "N": -2, "D": -3, "C": -1, "Q": 0,  "E": -2, "G": -3, "H": -2, "I": 1,  "L": 2,  "K": -1, "M": 5,  "F": 0,  "P": -2, "S": -1, "T": -1, "W": -1, "Y": -1, "V": 1},
    "F": {"A": -2, "R": -3, "N": -3, "D": -3, "C": -2, "Q": -3, "E": -3, "G": -3, "H": -1, "I": 0,  "L": 0,  "K": -3, "M": 0,  "F": 6,  "P": -4, "S": -2, "T": -2, "W": 1,  "Y": 3,  "V": -1},
    "P": {"A": -1, "R": -2, "N": -2, "D": -1, "C": -3, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -3, "L": -3, "K": -1, "M": -2, "F": -4, "P": 7,  "S": -1, "T": -1, "W": -4, "Y": -3, "V": -2},
    "S": {"A": 1,  "R": -1, "N": 1,  "D": 0,  "C": -1, "Q": 0,  "E": 0,  "G": 0,  "H": -1, "I": -2, "L": -2, "K": 0,  "M": -1, "F": -2, "P": -1, "S": 4,  "T": 1,  "W": -3, "Y": -2, "V": -2},
    "T": {"A": 0,  "R": -1, "N": 0,  "D": -1, "C": -1, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1,  "T": 5,  "W": -2, "Y": -2, "V": 0},
    "W": {"A": -3, "R": -3, "N": -4, "D": -4, "C": -2, "Q": -2, "E": -3, "G": -2, "H": -2, "I": -3, "L": -2, "K": -3, "M": -1, "F": 1,  "P": -4, "S": -3, "T": -2, "W": 11, "Y": 2,  "V": -3},
    "Y": {"A": -2, "R": -2, "N": -2, "D": -3, "C": -2, "Q": -1, "E": -2, "G": -3, "H": 2,  "I": -1, "L": -1, "K": -2, "M": -1, "F": 3,  "P": -3, "S": -2, "T": -2, "W": 2,  "Y": 7,  "V": -1},
    "V": {"A": 0,  "R": -3, "N": -3, "D": -3, "C": -1, "Q": -2, "E": -2, "G": -3, "H": -3, "I": 3,  "L": 1,  "K": -2, "M": 1,  "F": -1, "P": -2, "S": -2, "T": 0,  "W": -3, "Y": -1, "V": 4},
}

#Biotools augmented BLOSUM Comparison Table
BIOTOOLS_COMPARE = deepcopy(BLOSUM62)
BIOTOOLS_COMPARE["X"] = {}
BIOTOOLS_COMPARE["X"]["X"] = 25
for aa in AMINO_ACIDS:
    BIOTOOLS_COMPARE["X"][aa] = -25
    BIOTOOLS_COMPARE[aa]["X"] = -25



# Affine gap scoring for protein alignment
PROTEIN_GAP_OPEN = -10
PROTEIN_GAP_EXTEND = -1

# MSA gap scoring
PROTEIN_GAP = -8

# DNA Pariwise Matching Constants
DNA_MATCH = 8
DNA_MISMATCH = -5
DNA_GAP = -3




MY_HELP=f'''{BLUE_HELP}usage: {MAGENTA_HELP}biotools <command> [options]{RESET}

Command line utilities for DNA and protein sequence analysis.

{BLUE_HELP}options:{RESET}
    {GREEN_HELP}-h{RESET}, {CYAN_HELP}--help{RESET}     Show this help message and exit.

{BLUE_HELP}commands:{RESET}

    {BLUE_HELP}General{RESET}
    {GREEN_HELP}find_pos{RESET}       Find the position of a search string in a protein or DNA sequence, if present.
    {GREEN_HELP}length{RESET}         Return the length of a protein or DNA sequence.
  
    {BLUE_HELP}DNA{RESET}
    {GREEN_HELP}transcribe{RESET}     Transcribe a DNA sequence into RNA.
    {GREEN_HELP}translate{RESET}      Translate a DNA sequence into protein.
    {GREEN_HELP}rcomp{RESET}          Return the reverse complement of a DNA sequence.
    {GREEN_HELP}gc{RESET}             Return information about the base content of a DNA sequence.
    {GREEN_HELP}dna_align{RESET}      Align two given DNA sequences.

    {BLUE_HELP}Protein{RESET}
    {GREEN_HELP}charge{RESET}         Return the net charge of a protein sequence.
    {GREEN_HELP}mw{RESET}             Return the molecular weight of a protein sequence.
    {GREEN_HELP}ec{RESET}             Return the extinction coefficient of a protein sequence.
    {GREEN_HELP}props{RESET}          Return the molecular weight, extinction coefficient, and charge of a protein sequence.
    {GREEN_HELP}glycosylation{RESET}  Returns the protein sequence with N-linked glycosylation sites highlighted, if they are present.
    {GREEN_HELP}motifs{RESET}         Find regions of amino acid similarity between two protein sequences.
    {GREEN_HELP}protein_align{RESET}  Align two given protein sequences.
    {GREEN_HELP}msa{RESET}            Align multiple protein sequences, provided in a FASTA file or a FASTA formmated string.
    '''
