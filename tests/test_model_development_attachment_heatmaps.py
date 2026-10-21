import numpy as np
import pytest
from matplotlib.axes import Axes

from sim_swim.analysis.model_development_evaluation import (
    SCREEN_METRICS,
    _plot_attachment_patterns,
)


def rows():
    return [
        {
            "condition_id": cid,
            "n_flagella": 2,
            "attachment_pattern": "nf02__slots02",
            "attachment_pattern_label": "slots02",
            "screen_status": status,
            **{metric: 0.0 for metric, _ in SCREEN_METRICS},
        }
        for cid, status in [
            ("nf02__slots02__bfon", "pass"),
            ("nf02__slots02__bfoff", "fail"),
        ]
    ]


def test_interaction_arms_at_same_attachment_remain_separate(tmp_path, monkeypatch):
    matrices = []
    original = Axes.imshow

    def capture(self, data, *args, **kwargs):
        matrices.append(np.asarray(data).copy())
        return original(self, data, *args, **kwargs)

    monkeypatch.setattr(Axes, "imshow", capture)
    _plot_attachment_patterns(rows(), tmp_path / "paired.png")
    # Both statuses must remain visible, regardless of which arm came last.
    np.testing.assert_array_equal(matrices[0], [[0.0, 1.0]])


def test_duplicate_condition_identity_is_rejected(tmp_path):
    pair = rows()
    pair[1]["condition_id"] = pair[0]["condition_id"]
    with pytest.raises(ValueError, match="Duplicate condition IDs"):
        _plot_attachment_patterns(pair, tmp_path / "duplicate.png")
