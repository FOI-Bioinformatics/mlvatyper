from Bio.Seq import Seq

from mlvatyper.tandem import count_tandem_array, rotations

UNIT = "GATCGA"
FLANK_L = "TTTTTTTTTT"
FLANK_R = "CCCCCCCCCC"


def test_perfect_array():
    a = count_tandem_array(FLANK_L + UNIT * 5 + FLANK_R, (UNIT,))
    assert (a.copies, a.exact, a.variant, a.orientation) == (5, 5, 0, "+")
    assert a.start == len(FLANK_L) and a.span == 30


def test_variant_units_are_counted_with_tolerance():
    seq = FLANK_L + UNIT * 3 + "GATCGT" + UNIT * 4 + FLANK_R
    assert count_tandem_array(seq, (UNIT,)).copies == 8
    assert count_tandem_array(seq, (UNIT,), max_mismatch=0).copies == 4
    assert count_tandem_array(seq, (UNIT, "GATCGT"), max_mismatch=0).exact == 8


def test_reverse_orientation_and_rotation():
    seq = FLANK_L + UNIT * 6 + FLANK_R
    rc = str(Seq(seq).reverse_complement())
    a = count_tandem_array(rc, (UNIT,))
    assert a.copies == 6 and a.orientation == "-"
    rotated = UNIT[2:] + UNIT[:2]
    assert count_tandem_array(seq, (rotated,)).copies == 6
    assert len(rotations(UNIT)) == 6


def test_edge_cases():
    assert count_tandem_array("ACGT", ()).copies == 0
    assert count_tandem_array("ACG", (UNIT,)).copies == 0
    assert count_tandem_array(FLANK_L, (UNIT, "GATC")).copies == 0  # mixed unit lengths
    a = count_tandem_array(FLANK_L + "GATNGA" + FLANK_R, (UNIT,))
    assert a.copies == 1 and a.variant == 1
