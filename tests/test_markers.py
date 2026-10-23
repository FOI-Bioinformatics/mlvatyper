import pytest

from mlvatyper.markers import (
    MarkerError, check_against_sts, markers_from_ids, merge_definitions, parse_marker_name, read_markers,
    write_markers,
)
from mlvatyper.models import MarkerDef


def test_parse_mlvabank_style_name():
    d = parse_marker_name("Bruce06_134bp_408bp_3u")
    assert (d.unit_length, d.ref_size, d.ref_repeats) == (134, 408, 3)
    assert d.can_call and d.flanking == 408 - 3 * 134
    assert parse_marker_name("Bruce06") is None
    assert parse_marker_name("Lp03_96bp_941bp_8U").ref_repeats == 8


def test_read_markers_full_and_legacy_columns(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text(
        "marker\tunit_length\tref_size\tref_repeats\toffset\tpatterns\n"
        "m1\t6\t70\t5\t\tGATCGA;GATCGT\n"
        "m2\t9\t\t\t\tN/A\n"
        "X_12bp_300bp_4u\t\t\t\t\t\n"
    )
    defs = {d.primer_id: d for d in read_markers(p)}
    assert defs["m1"].can_call and defs["m1"].patterns == ("GATCGA", "GATCGT") and defs["m1"].flanking == 40
    assert defs["m1"].expected_size(7) == 82
    assert not defs["m2"].can_call and defs["m2"].unit_length == 9
    assert defs["X_12bp_300bp_4u"].can_call and defs["X_12bp_300bp_4u"].unit_length == 12
    legacy = tmp_path / "legacy.tsv"
    legacy.write_text("marker\tpattern\tlength\nm1\tGATCGA\t6\n")
    assert read_markers(legacy)[0].unit_length == 6


def test_read_markers_rejects_bad_values(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("id\tpatterns\tunit_length\nm1\tGATCGA\tsix\n")
    with pytest.raises(MarkerError):
        read_markers(p)
    p.write_text("marker\tunit_length\tref_size\tref_repeats\nm1\t10\t50\t6\n")
    with pytest.raises(MarkerError):
        read_markers(p)
    p.write_text("foo\tbar\nm1\tx\n")
    with pytest.raises(MarkerError):
        read_markers(p)


def test_merge_and_check(tmp_path):
    file_defs = [MarkerDef("A_6bp_60bp_3u", ("GATCGA",), 6), MarkerDef("B", (), 5)]
    from_ids = markers_from_ids(["A_6bp_60bp_3u", "C_4bp_40bp_2u", "D"])
    merged = {d.primer_id: d for d in merge_definitions(file_defs, from_ids)}
    assert merged["A_6bp_60bp_3u"].can_call and merged["A_6bp_60bp_3u"].patterns == ("GATCGA",)
    assert "C_4bp_40bp_2u" in merged and merged["B"].unit_length == 5
    kept = check_against_sts(list(merged.values()), {"A_6bp_60bp_3u", "B"})
    assert [d.primer_id for d in kept] == ["A_6bp_60bp_3u", "B"]
    out = tmp_path / "w.tsv"
    write_markers(kept, out)
    back = {d.primer_id: d for d in read_markers(out)}
    assert back == {d.primer_id: d for d in kept}


def test_read_headerless_legacy_pattern_file(tmp_path, caplog):
    from pathlib import Path

    from mlvatyper.markers import read_markers
    p = tmp_path / "pattern.txt"
    p.write_text("M1\tAAT\t3\nM2\tATCCTTATT;TTGCTCATT\t9\t1\t10\t1.0\nM3\tN/A\nM4\tACGTAC\t5\n")
    defs = {d.primer_id: d for d in read_markers(p)}
    assert list(defs) == ["M1", "M2", "M3", "M4"]  # the first line is data, not a skipped header
    assert defs["M1"].patterns == ("AAT",) and defs["M1"].unit_length == 3
    assert defs["M2"].patterns == ("ATCCTTATT", "TTGCTCATT") and defs["M2"].unit_length == 9
    assert not defs["M3"].has_repeats
    assert defs["M4"].unit_length == 6 and "stated unit length 5" in caplog.text
    repo = Path(__file__).resolve().parents[1] / "examples" / "testM" / "testM_pattern.txt"
    assert len(read_markers(repo)) == 31
