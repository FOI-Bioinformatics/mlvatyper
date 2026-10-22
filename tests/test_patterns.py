from mlvatyper.patterns import parse_pattern_line, read_patterns


def test_parse_pattern_line_variants():
    d = parse_pattern_line("m1\tgatcga;ACGT\t6\t1\t10\t0.5", 2)
    assert d.patterns == ("GATCGA", "ACGT") and d.pattern_length == 6 and d.has_repeats
    assert parse_pattern_line("m2\tN/A", 3).has_repeats is False
    assert parse_pattern_line("m3\t", 4).has_repeats is False
    d = parse_pattern_line("m4\tACGTA", 5)
    assert d.pattern_length == 5
    d = parse_pattern_line("m5\tACGTA\tx", 6)
    assert d.pattern_length == 5  # non-integer length falls back to first pattern length
    assert parse_pattern_line("only_one_column", 7) is None


def test_read_patterns_skips_header_comments_and_blank_lines(tmp_path):
    p = tmp_path / "p.tsv"
    p.write_text("marker\tpattern\n# comment\n\nm1\tACG\t3\nbad\n")
    defs = read_patterns(p)
    assert [d.primer_id for d in defs] == ["m1"]
