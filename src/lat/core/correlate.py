"""Time-correlation helpers: ACF, block averaging, linear fits."""

from __future__ import annotations

import numpy as np

from .errors import LatError


def autocorrelation(series, *, normalize: bool = True, unbiased: bool = False) -> np.ndarray:
    """Autocorrelation C(k) = <x(t) x(t+k)> via FFT."""
    x = np.asarray(series, float)
    if x.ndim != 1:
        raise LatError("autocorrelation expects a 1D series")
    n = len(x)
    if n == 0:
        raise LatError("autocorrelation of an empty series")
    dx = x - x.mean()
    size = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(dx, size)
    acf = np.fft.irfft(f * np.conj(f), size)[:n]
    if unbiased:
        acf = acf / np.arange(n, 0, -1)
    else:
        acf = acf / n
    if normalize:
        if acf[0] == 0:
            raise LatError("cannot normalize a zero-variance series")
        acf = acf / acf[0]
    return acf


def _acf_fft(x: np.ndarray) -> np.ndarray:
    n = len(x)
    size = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(x, size)
    return np.fft.irfft(f * np.conj(f), size)[:n]


def block_average(series, nblocks: int):
    """Return (mean, standard_error) using nblocks equal blocks."""
    x = np.asarray(series, float)
    if nblocks < 2:
        raise LatError("block averaging needs >= 2 blocks")
    if len(x) < nblocks:
        raise LatError(f"cannot split {len(x)} samples into {nblocks} blocks")
    edges = np.linspace(0, len(x), nblocks + 1).astype(int)
    means = np.array([x[edges[i]:edges[i + 1]].mean() for i in range(nblocks)])
    return float(means.mean()), float(means.std(ddof=1) / np.sqrt(nblocks))


def linear_fit(x, y):
    """Least-squares y = a*x + b; returns (slope, intercept, r2)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) < 2:
        raise LatError("linear fit needs >= 2 points")
    a, b = np.polyfit(x, y, 1)
    pred = a * x + b
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 if ss_tot == 0 else 1.0 - ss_res / ss_tot
    return float(a), float(b), r2


def integrate(y, x) -> float:
    return float(np.trapezoid(np.asarray(y, float), np.asarray(x, float)))
