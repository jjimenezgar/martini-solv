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

# Keep trajectory legend entries synchronized with components that are both
# present in the built system and enabled in the trajectory viewer.
from .trajectory_legend import apply as _apply_trajectory_legend

_apply_trajectory_legend()
del _apply_trajectory_legend

# Keep trajectory GIFs visually clean: show the connected protein backbone
# without detached side-chain dots while preserving all selected non-protein
# components and the PBC-safe coordinates from trajectory_pbc.
from .gif_backbone_only import apply as _apply_gif_backbone_only

_apply_gif_backbone_only()
del _apply_gif_backbone_only

# Accept modern protein structure uploads (PDB plus mmCIF/CIF, including
# AlphaFold downloads) while preserving the app's existing PDB-based pipeline.
from .structure_upload import install_streamlit_structure_upload_support as _install_structure_upload

_install_structure_upload()
del _install_structure_upload
