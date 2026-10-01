"""Plotting: read lat CSV, never recompute. Compute/plot decoupled (M14)."""

from .plot1d import plot_csv as plot1d
from .plot2d import plot_grid as plot2d

__all__ = ["plot1d", "plot2d"]
