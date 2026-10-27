"""Mg-alloy / PEO / HA implant vs residual bone-tumour simulation."""

from .model import OUTCOMES, OUTCOME_HELP, simulate
from .scenarios import build_params, grid, monte_carlo_params

__all__ = ["simulate", "build_params", "grid", "monte_carlo_params", "OUTCOMES", "OUTCOME_HELP"]
