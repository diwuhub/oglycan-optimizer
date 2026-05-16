import pytest

from oglycan.site_scoring import score_site


SCORE_CASES = [
    (0.95, 0.85, 0.30, 0.90, 0.95, 0.98),
    (0.80, 0.70, 0.20, 0.75, 0.82, 0.88),
    (0.60, 0.50, 0.10, 0.65, 0.70, 0.72),
]


def _site(difficulty: float, n_core_types: int) -> dict:
    return {
        "difficulty": difficulty,
        "core_types": [f"Core{i}" for i in range(1, n_core_types + 1)],
    }


@pytest.mark.parametrize("scores", SCORE_CASES)
@pytest.mark.parametrize("n_core_types", [1, 2, 3])
def test_harder_site_scores_no_higher_with_same_core_count(scores, n_core_types):
    easier_site = _site(0.4, n_core_types)
    harder_site = _site(0.9, n_core_types)

    assert score_site(harder_site, *scores) <= score_site(easier_site, *scores)


@pytest.mark.parametrize("scores", SCORE_CASES)
def test_difficulty_penalty_beats_extra_form_reward_in_regression_case(scores):
    easier_site = _site(0.4, 1)
    harder_site = _site(0.9, 4)

    assert score_site(harder_site, *scores) <= score_site(easier_site, *scores)
