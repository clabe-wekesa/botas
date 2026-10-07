from core.scoring import score_to_mapq


def test_exact_unique_alignment_has_high_mapq():
    assert score_to_mapq(0, read_len=100, second_score=None) == 60


def test_equal_best_and_second_best_is_ambiguous():
    assert score_to_mapq(-2, read_len=100, second_score=-2) == 0


def test_score_gap_increases_mapq():
    small_gap = score_to_mapq(-2, read_len=100, second_score=-3)
    large_gap = score_to_mapq(-2, read_len=100, second_score=-12)
    assert 0 < small_gap < large_gap <= 60


def test_more_edits_reduce_mapq_when_uniqueness_is_equal():
    good = score_to_mapq(-1, read_len=100, second_score=None)
    poor = score_to_mapq(-15, read_len=100, second_score=None)
    assert 0 <= poor < good <= 60


def test_invalid_read_length_returns_zero():
    assert score_to_mapq(0, read_len=0, second_score=None) == 0
