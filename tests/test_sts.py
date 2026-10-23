import pytest

from mlvatyper.sts import StsError, parse_sts_line, read_sts, write_sts


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


def test_real_sts_file_has_one_invalid_entry():
    from pathlib import Path
    sts = Path(__file__).resolve().parents[1] / "testID.sts"
    if not sts.exists():
        pytest.skip("testID.sts not present")
    valid, invalid = read_sts(sts)
    assert len(valid) == 207
    assert [i[0] for i in invalid] == [189]


def test_size_range_and_fix(tmp_path):
    pair = parse_sts_line("id\tACGT\tTTGC\t200-220", 1)
    assert (pair.expected_size, pair.size_low, pair.size_high, pair.size_text) == (210, 200, 220, "200-220")
    with pytest.raises(StsError, match="low-high"):
        parse_sts_line("id\tACGT\tTTGC\t220-200", 1)
    with pytest.raises(StsError, match="'1' at position 1"):
        parse_sts_line("id\t1ACGT\tTTGC\t100", 1)
    fixed = parse_sts_line("id\t1ACGT\tTT*GC\t100", 1, fix=True)
    assert (fixed.forward, fixed.reverse) == ("ACGT", "TTGC")
    p = tmp_path / "r.sts"
    p.write_text("a\tACGTACGTACGT\tTTGCTTGCTTGC\t200-220\n")
    valid, invalid = read_sts(p)
    write_sts(valid, tmp_path / "w.sts")
    assert (tmp_path / "w.sts").read_text() == "a\tACGTACGTACGT\tTTGCTTGCTTGC\t200-220\n" and not invalid
