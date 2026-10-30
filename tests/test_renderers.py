import biotools.renderers as btr
import biotools.models as btm
from PIL import Image

def test_save_dna_alignment_image(tmp_path):
    output = tmp_path / "alignment.png"

    btr.save_dna_pairwise_alignment(
        [list("ACGT"), list("ACCT")],
        filename=str(output),
    )

    with Image.open(output) as image:
        assert image.format == "PNG"
        assert image.width > 0
        assert image.height > 0
        image.load()


def test_save_protein_alignment_image(tmp_path):
    output = tmp_path / "alignment.png"

    profile = btm.Profile(
        name = ["name1", "name2"],
        alignment = [["-","Q","K","P","L","E"],["N","Q","R","-","L","E"]]
    )

    image = btr.render_protein_alignment_image(profile, scale=1, wrap=100, start=None, consensus=False, diff=False)
    btr.save_protein_alignment_image(image, suffix=".png", filename=str(output))

    with Image.open(output) as image:
        assert image.format == "PNG"
        assert image.width > 0
        assert image.height > 0
        image.load()