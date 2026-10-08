"""Protein and solvent system preparation for Martini 3."""

# Install visualization refinements while keeping the existing public API.
from .visualization_overrides import apply as _apply_visualization_overrides

_apply_visualization_overrides()
del _apply_visualization_overrides

# Reconstruct proteins across periodic boundaries before trajectory display/GIF
# alignment, preventing wrapped coordinates from looking like unfolding.
from .trajectory_pbc import apply as _apply_trajectory_pbc

_apply_trajectory_pbc()
del _apply_trajectory_pbc
