import argparse, sys
lazy from pathlib import Path
lazy from . import msa as btmsa
lazy from . import constants as btc
lazy from . import renderers as btr
lazy from . import sequences as bts
lazy from . import comparisons as btcomp
lazy from . import pairwise_align as bta
lazy from . import models as btm
lazy from . import pairwise_align as btpa
lazy import biotools_rust as btrust


def main():

    # Argument Parser Construction

    parser = argparse.ArgumentParser(
        prog="biotools",
        usage="biotools <command> [options]",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-h","--help",
        action="store_true"
        )
    
    subparsers = parser.add_subparsers(
        title="programs",        
        metavar=""
        )

    parser_transcribe = subparsers.add_parser(
        "transcribe",
        help="Transcribe a DNA sequence into RNA."
    )
    parser_transcribe.add_argument(
        "dna",
        help="Provide a DNA sequence to transcribe into RNA."
    )
    parser_transcribe.set_defaults(func=cli_transcribe)


    parser_translate = subparsers.add_parser(
        "translate",
        help="Translate a DNA sequence into protein."
    )
    parser_translate.add_argument(
        "dna",
        help="Provide a DNA sequence to translate into protein."
    )
    parser_translate.set_defaults(func=cli_translate)


    parser_rcomp = subparsers.add_parser(
        "rcomp",
        help="Return the reverse complement of a DNA sequence."
    )
    parser_rcomp.add_argument(
        "dna",
        help="Provide a DNA sequence to receive back the reverse complement sequence."
    )
    parser_rcomp.set_defaults(func=cli_reverse_comp)


    parser_gc = subparsers.add_parser(
        "gc",
        help = "Return the GC content of a DNA sequence."
    )
    parser_gc.add_argument(
        "dna",
        help="Provide a DNA sequence to return the GC content."
    )
    parser_gc.set_defaults(func=cli_gc)


    parser_charge = subparsers.add_parser(
        "charge",
        help="Return the net charge of a protein sequence."
    )
    parser_charge.add_argument(
        "protein",
        help="Provide a protein sequence and receive back the net charge of that protein."
    )
    parser_charge.set_defaults(func=cli_charge)


    parser_mw = subparsers.add_parser(
        "mw",
        help="Return the molecular weight of a protein sequence."
    )
    parser_mw.add_argument(
        "protein",
        help="Provide a protein sequence and receive back the molecular weight of that protein."
    )
    parser_mw.set_defaults(func=cli_mw)


    parser_ec = subparsers.add_parser(
        "ec",
        help="Return the extinction coefficient of a protein sequence."
    )
    parser_ec.add_argument(
        "protein",
        help="Provide a protein sequence and receive back the extinction coefficient of that protein."
    )
    parser_ec.set_defaults(func=cli_ec)


    parser_props = subparsers.add_parser(
        "props",
        help="Return the molecular weight, extinction coefficient, and charge of a protein sequence."
    )
    parser_props.add_argument(
        "protein",
        help="Provide a protein sequence and receive back the extinction coefficient, molecular weight, and charge of that protein."
    )
    parser_props.set_defaults(func=cli_props)


    parser_glycosylation = subparsers.add_parser(
        "glycosylation",
        help = "Returns the protein sequence with N-linked glycosylation sites highlighted, if they are present."
    )
    parser_glycosylation.add_argument(
        "protein",
        help="Provide a protein sequence to locate potential N-linked glycosylation sites."
    )
    parser_glycosylation.set_defaults(func=cli_glycosylation)


    parser_motifs = subparsers.add_parser(
        "motifs",
        help = "Given 2 protein sequences, looks for regions of amino acid similarity."
    )
    parser_motifs.add_argument(
        "s1",
        help="First sequence for comparison."
    )
    parser_motifs.add_argument(
        "s2",
        help="Second sequence for comparison."
    )
    parser_motifs.set_defaults(func=cli_motifs)


    parser_length = subparsers.add_parser(
        "length",
        help = "Given a protein or DNA sequence, returns the length of the sequence."
    )
    parser_length.add_argument(
        "sequence",
        help = "Sequence to determine the length of."
    )
    parser_length.set_defaults(func=cli_length)


    parser_find_pos = subparsers.add_parser(
        "find_pos",
        help = "Given a protein or DNA sequence, and a search string, returns the position of the search string in the sequence, if it exists."
    )
    parser_find_pos.add_argument(
        "sequence",
        help = "Sequence to search in."
    )
    parser_find_pos.add_argument(
        "search_string",
        help = "Sequence to search for."
    )
    parser_find_pos.set_defaults(func=cli_find_pos)


    parser_pairwise_dna = subparsers.add_parser(
        "dna_align",
        help = "Align two given DNA sequences."
    )
    parser_pairwise_dna.add_argument(
        "dna1",
        help = "Sequence 1 for alignment."
    )
    parser_pairwise_dna.add_argument(
        "dna2",
        help = "Sequence 2 for alignment."
    )
    parser_pairwise_dna.add_argument(
        "--save",
        metavar = "FILE",
        help = "Save a PNG image of the alignment with the provided filename."
    )
    parser_pairwise_dna.add_argument(
        "--diff",
        action="store_true",
        help = "This option will output dots for bases that are identical to the first sequence."
    )
    parser_pairwise_dna.add_argument(
        "--s1",
        metavar = "Seq 1 Start",
        type=int,
        default=1,
        help = "The base number that sequence 1 starts at. Defaults to 1."
    )
    parser_pairwise_dna.add_argument(
        "--s2",
        metavar = "Seq 2 Start",
        type=int,
        default=1,
        help = "The base number that sequence 2 starts at. Defaults to 1."
    )
    parser_pairwise_dna.set_defaults(func=cli_dna_align)


    parser_pairwise_protein = subparsers.add_parser(
            "protein_align",
            help = "Align two given protein sequences."
        )
    parser_pairwise_protein.add_argument(
        "seq1",
        help = "Sequence 1 for alignment."
    )
    parser_pairwise_protein.add_argument(
        "seq2",
        help = "Sequence 2 for alignment."
    )
    parser_pairwise_protein.add_argument(
        "--diff",
        action="store_true",
        help = "This option will output dots for bases that are identical to the first sequence."
    )
    parser_pairwise_protein.add_argument(
        "--name1",
        metavar = "NAME",
        help = "The name of sequence 1. Defaults to '1'."
    )
    parser_pairwise_protein.add_argument(
        "--name2",
        metavar = "NAME",
        help = "The name of sequence 2. Defaults to '2'."
    )
    parser_pairwise_protein.add_argument(
        "--wrap",
        metavar = "N",
        type=_positive_int,
        default=300,
        help = "Specfiy how many residues of the alignment are rendered per line (image only). Defaults to 300."
    )
    parser_pairwise_protein.add_argument(
        "--real-ruler",
        nargs="?",
        const = 1,
        type = _positive_int,
        default = None,
        metavar = "START",
        help = "This option will display the ruler based on the amino acid numbers of the reference sequence, starting at START (defaults to 1)."
    )
    parser_pairwise_protein.add_argument(
        "--save",
        nargs="?",
        const=True,
        metavar = "FILE",
        help = "Save alignment to FILE, or choose a filename interactively."
    )
    parser_pairwise_protein.add_argument(
        "--consensus",
        action = "store_true",
        help = "This option will add a consensus sequence to the bottom of the alignment."
    )
    parser_pairwise_protein.add_argument(
        "--backend",
        choices=("python", "rust"),
        default="rust",
        help = "Choose whether to use the rust implementation (compiled, default) or python (much slower for large numbers of sequences.)"
    )
    parser_pairwise_protein.set_defaults(func=cli_protein_align)


    parser_msa = subparsers.add_parser(
        "msa",
        help = "Aligned multiple protein sequences, provided in a fasta file or a fasta formmated string."
    )
    parser_msa.add_argument(
        "file",
        help = "Input fasta file, or a fasta-formatted string."
    )
    parser_msa.add_argument(
        "--diff",
        action = "store_true",
        help = "This option will output dots for amino acids that are identical to the first sequence."
    )
    parser_msa.add_argument(
        "--consensus",
        action = "store_true",
        help = "This option will add a consensus sequence to the bottom of the alignment."
    )
    parser_msa.add_argument(
        "--wrap",
        metavar = "N",
        type = _positive_int,
        default = 300,
        help = "Specify how many residues of the alignment are rendered per line (image only). Defaults to 300."
    )
    parser_msa.add_argument(
        "--real-ruler",
        nargs="?",
        const = 1,
        type = _positive_int,
        default = None,
        metavar = "START",
        help = "This option will display the ruler based on the amino acid numbers of the reference sequence, starting at START (defaults to 1)."
    )
    parser_msa.add_argument(
        "--save",
        nargs="?",
        const=True,
        metavar = "OUTPUT_FILE",
        help = "Save alignment to FILE, or choose a filename interactively."
    )
    parser_msa.add_argument(
        "--backend",
        choices=("python", "rust"),
        default="rust",
        help = "Choose whether to use the rust implementation (compiled, default) or python (much slower for large numbers of sequences.)"
    )
    parser_msa.add_argument(
        "--gap-extend",
        metavar = "N",
        type = _non_positive_int,
        default = -3,
        help = "Specify the profile merge level affine gap extend penalty. Defaults to -3."
    )
    parser_msa.add_argument(
        "--gap-open",
        metavar = "N",
        type = _non_positive_int,
        default = -10,
        help = "Specify the profile merge level affine gap open penalty. Defaults to -10."
        )
    parser_msa.set_defaults(func=cli_msa)


    parser_identity = subparsers.add_parser(
        "identity",
        help = "Returns all pairwise sequence identities from an input FASTA file or FASTA-formatted string."
    )
    parser_identity.add_argument(
        "file",
        help = "Input fasta file, or a fasta-formatted string."
    )
    parser_identity.set_defaults(func=cli_identity)


    args = parser.parse_args()

    if args.help:
        print(btc.MY_HELP)
        return

    if hasattr(args, "func"):
        try:
            args.func(args)
        except ValueError as error: 
            print(error, file=sys.stderr)
            sys.exit(1)
    else:
        parser.print_help()


# CLI Function Dispatchers

def cli_transcribe(args):
    print(f"{btc.RED_TERM}{bts.transcribe(args.dna)}{btc.RESET}")


def cli_translate(args):
    f1, f2, f3 = bts.translate(args.dna)
    print(f"{btc.RED_TERM}Frame 1: {f1}")
    print(f"Frame 2: {f2}")
    print(f"Frame 3: {f3}{btc.RESET}")


def cli_reverse_comp(args):
    print(f"{btc.RED_TERM}{bts.reverse_comp(args.dna)}{btc.RESET}")


def cli_gc(args):
    length, contents, at_content, gc_content = bts.gc(args.dna)
    print(f"{btc.RED_TERM}Total bases:{btc.RESET} {length}")
    print(f"{btc.RED_TERM}A:{btc.RESET} {contents["A"]}")
    print(f"{btc.RED_TERM}C:{btc.RESET} {contents["C"]}")
    print(f"{btc.RED_TERM}G:{btc.RESET} {contents["G"]}")
    print(f"{btc.RED_TERM}T:{btc.RESET} {contents["T"]}")
    print(f"{btc.RED_TERM}AT Content:{btc.RESET} {at_content:.2f}")
    print(f"{btc.RED_TERM}GC Content:{btc.RESET} {gc_content:.2f}")


def cli_charge(args):
    protein = btm.Protein(name="", sequence=args.protein)
    print(f"{btc.RED_TERM}Charge = {protein.charge()}{btc.RESET}")


def cli_mw(args):
    protein = btm.Protein(name="", sequence=args.protein)
    print(f"{btc.RED_TERM}Molecular Weight = {protein.mw():.1f}{btc.RESET}")


def cli_ec(args):
    protein = btm.Protein(name="", sequence=args.protein)
    print(f"{btc.RED_TERM}Extinction Coefficient = {protein.ec()}{btc.RESET}")


def cli_props(args):
    protein = btm.Protein(name="", sequence=args.protein)
    print(f"{btc.RED_TERM}Molecular Weight = {protein.mw():.1f}{btc.RESET}")
    print(f"{btc.RED_TERM}Charge = {protein.charge()}{btc.RESET}")
    print(f"{btc.RED_TERM}Extinction Coefficient = {protein.ec()}{btc.RESET}")


def cli_glycosylation(args):
    protein = btm.Protein(name="", sequence=args.protein)
    sites = protein.glycosylation()
    if len(sites) == 0:
        print("No N-linked glycosylation sites found.")
    else:    
        for i, aa in enumerate(args.protein):
            aa = aa.upper()
            if i not in sites and i-1 not in sites and i-2 not in sites:
                print(aa, end="")
            else:
                print(f"{btc.RED_TERM}{aa}{btc.RESET}", end="")
        print()


def cli_motifs(args):
    protein1 = btm.Protein(name="", sequence=args.s1)
    protein2 = btm.Protein(name="", sequence=args.s2)
    similars = btcomp.motifs(protein1, protein2)
    if not similars:
        print("No significant similarity found.")
    for key in similars:
        print(f"{btc.RED_TERM}Peptide size = {key[0]}. Threshold of similarity = {key[1]}.{btc.RESET}")
        for ((peptide1, start1, end1), (peptide2, start2, end2), score) in similars[key]:
            print(f"{start1}-{peptide1}-{end1}")
            print(f"{start2}-{peptide2}-{end2}")
            print(f"Score: {score}/{key[0]} ({score/key[0]:.1%})")
            print()


def cli_length(args):
    print(f"Length: {btc.RED_TERM}{bts.length(args.sequence)}{btc.RESET}")


def cli_find_pos(args):
    hits = bts.find_pos(args.sequence, args.search_string)
    if not hits:
        print("Not found")
    else:
        print("Found!")
        for hit in hits:
            print(f"At position: {btc.RED_TERM}{hit}{btc.RED_TERM}")


def cli_dna_align(args):
    align1, align2 = bta.pairwise_dna(args.dna1, args.dna2)
    btr.print_dna_pairwise_alignment([align1, align2], args.diff)

    if args.save:
        btr.save_dna_pairwise_alignment([align1, align2], filename=args.save, start1=args.s1, start2=args.s2, diff=args.diff)


def cli_protein_align(args):
    protein1 = btm.Protein("", args.seq1, allow_x=True)
    protein2 = btm.Protein("", args.seq2, allow_x=True)

    if args.backend == "rust":
        align1, align2 = btrust._align_protein(protein1.sequence, protein2.sequence)
        align1 = list(align1)
        align2 = list(align2)
    else:
        align1, align2 = btpa.pairwise_protein(protein1, protein2)

    profile = btm.Profile(
        name = [args.name1 or "1", args.name2 or "2"],
        alignment = [align1, align2]
    )

    if args.consensus:
        profile = btmsa.add_consensus(profile)
        consensus = True
    else:
        consensus = False

    _output_protein_alignment(profile, args, consensus)


def cli_msa(args):
    profile = btmsa.msa(args.file, args.backend, int(args.gap_open), int(args.gap_extend))
    if args.consensus:
        profile = btmsa.add_consensus(profile)
        consensus = True
    else:
        consensus = False
    _output_protein_alignment(profile, args, consensus)


def cli_identity(args):
    identities = btcomp.identity(args.file)
    for key in identities:
        print(key, f"{identities[key]:.3f}")


def _output_protein_alignment(profile:btm.Profile, args, consensus:bool):
    if args.save:
        filename = args.save

        # If no filename provided, use dialog window to choose location and filename.
        if filename is True:
            from . import dialogs as btdi
            filename = btdi.get_save_filename()
        if not filename:
            print("Save cancelled.")
            return

        # Determine the requested image format
        suffix = Path(filename).suffix.lower()
        if suffix == ".png":
            scale = 1
        elif suffix in {".tif", ".tiff"}:
            scale = 4
        else:
            raise ValueError(f"Unsupported image format {suffix}")

        # Render the alignment properly, then save it as an image
        image = btr.render_protein_alignment_image(profile, scale, wrap=args.wrap, start=args.real_ruler, consensus=consensus, diff=args.diff)
        btr.save_protein_alignment_image(image, suffix, filename=filename)

    else:

        # Terminal output
        btr.print_protein_alignment(profile, consensus, args.diff)


def _positive_int(value:str):
    try:
        value = int(value)
        if value < 1:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Ruler start must be a positive integer.")
    return value


def _non_positive_int(value:str):
    try:
        value = int(value)
        if value > 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("Gap penalty must be a non positive integer.")
    return value


    
if __name__ == "__main__":
    main()