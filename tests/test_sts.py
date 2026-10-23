import pytest

from mlvatyper.sts import StsError, parse_sts_line, read_sts


def test_parse_valid_line():
    pair = parse_sts_line("id1\tacgt\tTTGC\t250\tsome alias", 1)
    assert pair.id == "id1"
    assert pair.forward == "ACGT"
    assert pair.expected_size == 250
    assert pair.alias == "some alias"


def test_reject_non_acgt_primer():
    with pytest.raises(StsError, match="non-IUPAC"):
        parse_sts_line("id\t1ACGT\tTTGC\t250", 7)


def test_iupac_codes_are_accepted():
    from mlvatyper.sts import is_degenerate
    pair = parse_sts_line("id\tACGTR\tTTGCN\t250", 1)
    assert is_degenerate(pair)
    assert not is_degenerate(parse_sts_line("id\tACGT\tTTGC\t250", 1))


def test_reject_non_integer_size():
    with pytest.raises(StsError, match="not a non-negative integer"):
        parse_sts_line("id\tACGT\tTTGC\tabc", 1)


def test_read_sts_reports_invalid_and_duplicates(tmp_path):
    p = tmp_path / "x.sts"
    p.write_text("a\tACGT\tTTGC\t100\r\nb\t1ACGT\tTTGC\t100\r\na\tACGT\tTTGC\t100\r\n")
    valid, invalid = read_sts(p)
    assert [v.id for v in valid] == ["a"]
    assert len(invalid) == 2
    assert "non-IUPAC" in invalid[0][2]
    assert "duplicate" in invalid[1][2]


def test_example_sts_and_markers():
    from pathlib import Path

    from mlvatyper.markers import read_markers
    example = Path(__file__).resolve().parents[1] / "examples" / "testM"
    valid, invalid = read_sts(example / "testM.sts")
    assert len(valid) == 31 and not invalid
    defs = {d.primer_id: d for d in read_markers(example / "testM_markers.tsv")}
    assert set(defs) == {p.id for p in valid}
    assert defs["testM03"].patterns[0] == "ATCCTTATT" and len(defs["testM03"].patterns) == 5
    assert all(d.has_repeats and not d.can_call for d in defs.values())
