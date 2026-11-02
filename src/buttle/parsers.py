from . import models as btm

def parse_fasta_file(
        file:str,
        allow_x:bool = False
    ) -> list[btm.Protein]:
    '''
    Private function for reading in the FASTA file and returning a list of Protein objects.
    '''
    with open(file, "r") as f:
            data = f.readlines()

    length = len(data)
    sequences = []
    name = ""
    sequence = ""
    for i in range(0, length):

        if data[i].startswith(">"):
            if sequence != "":
                new_seq = btm.Protein(
                    name = name,
                    sequence = sequence,
                    allow_x = allow_x
                )
                sequences.append(new_seq)
                sequence = ""
            name = data[i].replace(">", "").strip()
        else:
            sequence += data[i].strip().upper()

    new_seq = btm.Protein(
        name = name,
        sequence=sequence,
        allow_x = allow_x
    )
    sequences.append(new_seq)
        
    if len(sequences) < 2:
        raise ValueError("MSA alignment requires at least 2 sequences.")
    names = [protein.name for protein in sequences]
    if len(names) != len(set(names)):
        raise ValueError("FASTA sequence names must be unique.")

    return sequences


def parse_fasta_string(
        sequences: str,
        allow_x:bool=False
    ) -> list[btm.Protein]:
    '''
    Private function for accepting sequence data in a single-line-string formatted properly for FASTA.
    '''

    seq_list = sequences.splitlines()

    if not seq_list:
        raise ValueError("No sequences found.")

    if len(seq_list) % 2 != 0:
        raise ValueError("Each FASTA header must be followed by a sequence line.")

    parsed_sequences = []
    for i in range(0, len(seq_list), 2):
        header = seq_list[i]
        sequence = seq_list[i+1].strip()
        if not header.startswith(">") or sequence.startswith(">"):
            raise ValueError("Invalid input. Please format as $'>name1\\nsequence1\\n>name2\\nsequence2\\n>name3\\nsequence3'")

        if not sequence:
            raise ValueError(f"Missing sequence for header {header!r}")
        
        new_seq = btm.Protein(
            name = seq_list[i].replace(">", "").strip(),
            sequence = seq_list[i+1].strip().upper(),
            allow_x = allow_x
        )
        parsed_sequences.append(new_seq)

    if len(parsed_sequences) < 2:
        raise ValueError("MSA alignment requires at least 2 sequences.")
    names = [protein.name for protein in parsed_sequences]
    if len(names) != len(set(names)):
        raise ValueError("FASTA sequence names must be unique.")
    return parsed_sequences