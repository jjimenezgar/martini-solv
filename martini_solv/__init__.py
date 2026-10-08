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

# GōMartini CA virtual sites are force-field helpers, not physical CG beads.
from .go_virtual_visibility import apply as _apply_go_virtual_visibility

_apply_go_virtual_visibility()
del _apply_go_virtual_visibility
