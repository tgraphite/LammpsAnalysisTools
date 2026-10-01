"""Unit tests for plot1d column selection / bar geometry / normalization."""

import numpy as np
import pytest

from lat.core.errors import LatError
from lat.plotting.plot1d import _bar_geometry, _normalize, _select_y


COLS = ["frame", "A_SSIP", "A_CIP", "A_AGG"]


def test_select_y_defaults_to_all_y_columns():
    assert _select_y(COLS, None) == [(1, "A_SSIP"), (2, "A_CIP"), (3, "A_AGG")]


def test_select_y_preserves_user_order_and_dedupes():
    assert _select_y(COLS, ["A_AGG", "A_SSIP", "A_AGG"]) == [(3, "A_AGG"), (1, "A_SSIP")]


def test_select_y_ignores_x_name_and_rejects_unknown():
    with pytest.raises(LatError):
        _select_y(COLS, ["frame"])          # x is always drawn, not selectable
    with pytest.raises(LatError):
        _select_y(COLS, ["nope"])


def test_normalize_gives_occupancy_fractions():
    ys = np.array([[1.0, 0.0, 2.0], [1.0, 1.0, 2.0]])
    frac = _normalize(ys)
    np.testing.assert_allclose(frac.sum(axis=0), [1.0, 1.0, 1.0])
    np.testing.assert_allclose(frac[:, 0], [0.5, 0.5])
    np.testing.assert_allclose(frac[:, 1], [0.0, 1.0])
    # an all-zero group must not divide by zero
    np.testing.assert_allclose(_normalize(np.zeros((2, 3))), 0.0)


def test_bar_geometry_no_overlap_and_single_point():
    width, off = _bar_geometry(np.arange(4), 3)
    assert off[0] < 0 < off[-1]
    assert off[-1] - off[0] == pytest.approx(2 * width)
    assert 3 * width < 1.0                    # grouped bars fit inside the spacing
    width1, off1 = _bar_geometry(np.array([5.0]), 3)
    assert width1 == width and off1.tolist() == off.tolist()
