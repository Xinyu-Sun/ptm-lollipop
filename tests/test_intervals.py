from ptm_lollipop.intervals import compare_disorder, merge_intervals, parse_ranges


def test_parse_ranges_accepts_common_separators():
    assert parse_ranges("1-44; 449-566") == ((1, 44), (449, 566))
    assert parse_ranges("1-44 | 449-566") == ((1, 44), (449, 566))
    assert parse_ranges("no") == ()


def test_merge_intervals_merges_overlaps_but_not_gaps():
    assert merge_intervals([(5, 10), (1, 4), (10, 20), (30, 40)]) == ((1, 4), (5, 20), (30, 40))


def test_compare_disorder_flags_near_match():
    note, status = compare_disorder("10-30", ((12, 30),), 100)
    assert status == "near match"
    assert "similar" in note
