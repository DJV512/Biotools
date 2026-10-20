from . import constants as btc
from . import models as btm
lazy from PIL import Image, ImageDraw, ImageFont
lazy from collections import Counter
lazy from math import ceil

####################################
##### Public Facing Functions ######
####################################

def print_dna_pairwise_alignment(alignments:list[list[str]], diff:bool = False) -> None:
    '''
    Given a list of two aligned sequences, outputs a colorized "image" to the terminal.
    '''

    # Analyze the DNA alignments for matches, mismatches, and gaps
    analyzed = _analyze_dna_alignments(alignments, diff)


    # Render the lines into colorized strings
    lines = []
    for line in analyzed:
        rendered_line = "".join(f"{btc.TERMINAL_COLORS_DNA[status]}{char}" for status, char in line)
        lines.append(rendered_line)


    # Output to the terminal
    for line in lines:
        print(f"{btc.BG_WHITE}{line}{btc.RESET}")


def save_dna_pairwise_alignment(alignments:list[list[str]], filename:str, *, start1:int=1, start2:int=1, diff:bool=False) -> None:
    '''
    Handed a list of two aligned sequences, constructs and saves a PNG image of the alignment.
    '''

    # Analyze the DNA alignments for matches, mismatches, and gaps
    analyzed = _analyze_dna_alignments(alignments, diff)


    # Set Image parameters
    num_lines = len(alignments)
    char_count = len(alignments[0])
    font = ImageFont.truetype("/System/Library/Fonts/Menlo.ttc", 20)
    char_width = font.getlength("A")
    ascent, descent = font.getmetrics()
    line_height = ascent + descent
    image_width = int(char_count * char_width + 2 * btc.PADDING)
    image_height = int(num_lines * line_height + 2 * btc.PADDING + 2 * line_height)
    line_num = 0

    # Create Image
    image = Image.new(mode = "RGB", size = (image_width, image_height), color = (255,255,255))
    draw = ImageDraw.Draw(image)

    # Add top sequence numbers to the Image
    _draw_number_line(draw, font, line_num, start1, char_count, line_height, char_width)

    # Add sequences to the Image
    for line in analyzed:
        line_num += 1
        _draw_dna_sequence_line(draw, font, line, line_num, line_height, char_width)

    # Bottom sequence numbers to the Image
    line_num += 1  
    _draw_number_line(draw, font, line_num, start2, char_count, line_height, char_width)

    # Save Image
    image.save(filename)


def print_protein_alignment(
        final_alignment: btm.Profile,
        consensus: bool,
        diff:bool = False,
    ) -> None:
    '''
    Given a multiple sequence alignment Profile, outputs a colorized "image" to the terminal.
    '''

    # Analyze the protein alignments for matches, mismatches, conservative substitutions, and gaps
    parsed_profile = _analyze_msa(final_alignment, diff, consensus)

    max_length = max(len(x) for x in parsed_profile.name)

    # # Render the lines into colorized strings
    lines = []
    for i, (name, align) in enumerate(zip(parsed_profile.name, parsed_profile.alignment)):
        rendered_line = f"{btc.BLACK_TERM}{name:<{max_length+3}}" + "".join(f"{btc.TERMINAL_COLORS_PROTEIN[status]}{char}" for status, char in align)
        lines.append(rendered_line)


    # # Output to the terminal
    for line in lines:
        print(f"{btc.BG_WHITE}{line}{btc.RESET}")


def render_protein_alignment_image(
        final_alignment: btm.Profile,
        scale: int,
        wrap: int,
        start: int,
        consensus: bool,
        diff: bool = False,
    ) -> Image:
    '''
    Given a multiple sequence alignment Profile, renders and returns a PIL Image object.
    '''

    padding = btc.PADDING * scale

    # Analyze the protein alignments for matches, mismatches, conservative substitutions, and gaps
    parsed_profile = _analyze_msa(final_alignment, diff, consensus)
    max_name_length = max(len(x) for x in final_alignment.name)    

    # Set Image parameters
    name_width = max_name_length + 3
    num_lines = len(final_alignment.name)
    char_count = len(final_alignment.alignment[0])
    font = ImageFont.truetype("/System/Library/Fonts/Menlo.ttc", 30 * scale)
    char_width = font.getlength("A")
    ascent, descent = font.getmetrics()
    line_height = ascent + descent
    necessary_image_height = int(num_lines * line_height + padding + 2 * line_height)
    necessary_alignment_rows = ceil(char_count / wrap)
    image_height = necessary_image_height * necessary_alignment_rows
    max_image_width = int((name_width + wrap) * char_width + 2 * padding)
    necessary_image_width = int((name_width + char_count) * char_width + 2 * padding)
    image_width = min(max_image_width, necessary_image_width)
    line_num = 0
    sequence_number = 0
    position_start = 0
    if start is not None:
        residue_start = start
        alignment_ruler = True
    else:
        residue_start = 1
        alignment_ruler = False

    
    # Create Image
    image = Image.new(mode = "RGB", size = (image_width, image_height), color = (255,255,255))
    draw = ImageDraw.Draw(image)

    for _ in range(necessary_alignment_rows):

        # Add top sequence numbers to the Image
        current_residue = _draw_msa_number_line(draw, font, line_num, position_start, residue_start, wrap, line_height, char_width, name_width, char_count, padding, parsed_profile.alignment[0], alignment_ruler)
        residue_start = current_residue
        position_start += wrap
        line_num += 1

        # Add sequences to the Image
        for name, alignment in zip(parsed_profile.name, parsed_profile.alignment, strict=True):
            line_num += 1
            _draw_msa_line(draw, font, name, alignment, sequence_number, line_num, line_height, char_width, name_width, wrap, padding)

        sequence_number += wrap
        line_num += 2

    return image


def save_protein_alignment_image(
        image: Image,
        suffix: str,
        filename: str
    ) -> None:
    '''
    Accepts a rendered PIL image object and saves it in the proper format.
    '''
    if suffix == ".png":
        image.save(filename)
    elif suffix in [".tif", ".tiff"]:
        image.save(filename, format= "TIFF", compression="tiff_lzw", dpi=(600,600))




####################################
##### Private Helper Functions #####
####################################

def _analyze_dna_alignments(alignments:list[list[str]], diff:bool) -> list[list[tuple[str, str]]]:
    '''
    Using provided alignments, analzye them to determine whether each position results in a match, mismatch, or a gap.
    '''

    line1 = []
    line2 = []

    for base1, base2 in zip(*alignments, strict=True):
        if base1 == base2:
            status = "match"
        elif base1 == "-" or base2 == "-":
            status = "gap"
        else:
            status = "mismatch"

        line1.append((status, base1))
        if diff and base1 == base2:
            line2.append((status, "."))
        else:
            line2.append((status, base2))

    return [line1, line2]


def _draw_number_line(draw:ImageDraw.ImageDraw, font:ImageFont.FreeTypeFont, line_num:int, start:int, char_count:int, line_height:int, char_width:int) -> None:
    '''
    Adds a number line to the image. 
    '''
    for i, x in enumerate(range(start, start + char_count)):
        if x % 10 == 0:
            i = btc.PADDING + i * char_width
            draw.text((i, btc.PADDING + line_num * line_height), str(x), fill = btc.BLACK, font=font)


def _draw_dna_sequence_line(draw:ImageDraw.ImageDraw, font:ImageFont.FreeTypeFont, line:list[tuple[str, str]], line_num:int, line_height:int, char_width:int) -> None:
    '''
    Adds each sequence in the alignment to the image.
    '''
    for j, (status, char) in enumerate(line):
        x = btc.PADDING + j * char_width
        y = btc.PADDING + line_num * line_height
        draw.text((x, y), char, fill = btc.DNA_IMAGE_COLORS[status], font=font)


def _analyze_msa(
        alignment: btm.Profile,
        diff: bool,
        consensus: bool
    ) -> btm.ParsedProfile:
    '''
    Analyzes the multiple alignment and replaces matches in all sequences except the first with ".".
    '''

    number_of_sequences = len(alignment.alignment)
    parsed_profile = btm.ParsedProfile(name=alignment.name, alignment = [[] for _ in range(number_of_sequences)])

    first_sequence = alignment.alignment[0]
    for x, column in enumerate(zip(*alignment.alignment)):
        if consensus:
            count = Counter(column[:-1])
        else:
            count = Counter(column)
        if "-" in count:
            del count["-"]
        for i, aa in enumerate(column):
            if i == 0:
                if aa == "-":
                    parsed_profile.alignment[i].append(("gap", aa))
                elif aa == "X":
                    parsed_profile.alignment[i].append(("white", aa))
                elif count[aa] > 1:
                    parsed_profile.alignment[i].append(("match", aa))
                else:
                    skip = False
                    for each_aa in count:
                        if each_aa != aa:
                            conservatives = btc.CONSERVATIVES[aa]
                            if conservatives is not None and each_aa in conservatives:
                                parsed_profile.alignment[i].append(("cons", aa))
                                skip = True
                                break
                    if not skip:
                        parsed_profile.alignment[i].append(("non-cons", aa))

            else:
                conservatives = btc.CONSERVATIVES.get(first_sequence[x])
                if aa == "X":
                    parsed_profile.alignment[i].append(("white", aa))
                elif aa == first_sequence[x] and aa != "-" and consensus and i == number_of_sequences - 1:
                    parsed_profile.alignment[i].append(("match", aa))
                elif aa == first_sequence[x] and aa != "-" and diff:       
                    parsed_profile.alignment[i].append(("match", "."))
                elif aa == first_sequence[x] and aa != "-":
                    parsed_profile.alignment[i].append(("match", aa))
                elif conservatives is not None and aa in conservatives:
                    parsed_profile.alignment[i].append(("cons", aa))
                elif aa == "-":
                    parsed_profile.alignment[i].append(("gap", aa))
                else:
                    parsed_profile.alignment[i].append(("non-cons", aa))

    return parsed_profile


def _draw_msa_number_line(
        draw:ImageDraw.ImageDraw, 
        font:ImageFont.FreeTypeFont, 
        line_num:int, 
        position_start:int,
        residue_start:int, 
        wrap:int, 
        line_height:int, 
        char_width:int, 
        name_width:int, 
        char_count:int,
        padding:int,
        alignment:list[tuple[str,str]],
        alignment_ruler: bool
    ) -> int:
    '''
    Adds a number line to the image. 
    '''
    current_residue = residue_start

    for offset, position in enumerate(range(position_start, min(position_start + wrap, char_count))):
        if alignment_ruler:
            if alignment[position][1] == "X":
                current_residue = 0
            if alignment[position][1] != "-":
                if current_residue != 0 and current_residue % 10 == 0:
                    x = padding + name_width * char_width + offset * char_width
                    draw.text((x, padding + line_num * line_height), str(current_residue), fill = btc.BLACK, font=font)
                    draw.text((x, padding + (line_num + 1) * line_height), "|", fill = btc.BLACK, font=font)
                current_residue += 1
        else:
            if alignment[position][1] == "X":
                current_residue = 0
            if current_residue != 0 and current_residue % 10 == 0:
                x = padding + name_width * char_width + offset * char_width
                draw.text((x, padding + line_num * line_height), str(current_residue), fill = btc.BLACK, font=font)
                draw.text((x, padding + (line_num + 1) * line_height), "|", fill = btc.BLACK, font=font)
            current_residue += 1

    return current_residue
                    

def _draw_msa_line(
        draw:ImageDraw.ImageDraw,
        font:ImageFont.FreeTypeFont, 
        name:str, 
        alignment:list[tuple[str, str]], 
        start:int, 
        line_num:int, 
        line_height:int, 
        char_width:int, 
        name_width:int, 
        wrap:int,
        padding: int
    ) -> None:
    '''
    Adds each sequence in the alignment to the image.
    '''

    y = padding + line_num * line_height

    for j, char in enumerate(name):
        x = padding + char_width * j
        draw.text((x,y), char, fill = btc.BLACK, font=font)

    for offset, position in enumerate(range(start, min(len(alignment), start+wrap))):           
        (status, char) = alignment[position]
        x = padding + name_width * char_width + offset * char_width
        draw.text((x, y), char, fill = btc.PROTEIN_IMAGE_COLORS[status], font=font)


