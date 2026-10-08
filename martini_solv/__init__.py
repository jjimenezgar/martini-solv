"""Protein and solvent system preparation for Martini 3."""

# Install visualization refinements while keeping the existing public API.
from .visualization_overrides import apply as _apply_visualization_overrides

_apply_visualization_overrides()
del _apply_visualization_overrides
