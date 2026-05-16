"""Shared math utilities for sub-models."""

import math


def gauss(x: float, mu: float, sigma: float) -> float:
    """Gaussian response: 1.0 at x=mu, decaying with distance."""
    return math.exp(-0.5 * ((x - mu) / sigma) ** 2)
