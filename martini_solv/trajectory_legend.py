"""Keep the trajectory legend in sync with the actual system and viewer toggles."""
from __future__ import annotations

import json
from pathlib import Path

from . import molecular_viewer as mv


LEGEND_SPEC = (
    ("protein", "Protein", mv.BLUE),
    ("water", "Water", "#B0BEC5"),
    ("choline", "Choline", "#7E57C2"),
    ("urea", "Urea", "#FFB74D"),
    ("sorbitol", "Sorbitol", "#FF7043"),
    ("ions", "Ions", "limegreen"),
    ("solute", "Free molecule", "yellow"),
)


def _dynamic_trajectory_legend_html(
    present: dict[str, list[str]],
    *,
    show_protein: bool,
    show_solute: bool,
    show_solvent: bool,
    show_ions: bool,
) -> str:
    enabled = {
        "protein": show_protein,
        "water": show_solvent,
        "choline": show_solvent,
        "urea": show_solvent,
        "sorbitol": show_solvent,
        "ions": show_ions,
        "solute": show_solute,
    }
    items = []
    for key, label, color in LEGEND_SPEC:
        if enabled[key] and present.get(key):
            items.append(f'<span><i style="background:{color}"></i>{label}</span>')
    return "\n        ".join(items)


def render_trajectory(
    gro_path: Path,
    xtc_path: Path,
    *,
    stride: int = 1,
    height: int = 700,
    show_protein: bool = True,
    show_solute: bool = False,
    show_solvent: bool = False,
    show_ions: bool = False,
    reline_chloride_count: int = 0,
) -> int:
    """Render a trajectory whose legend reflects only present, visible components."""
    pdb, frames = mv._trajectory_as_multimodel_pdb(
        gro_path, xtc_path, stride, reline_chloride_count=reline_chloride_count
    )
    present = mv._component_resnames(gro_path)

    # Only expose 3-character PDB viewer labels for components that are truly
    # present in the source GRO.  This prevents e.g. Urea from appearing in a
    # water-only system just because URE is a known renderer label.
    components_map = {
        "water": ["WAT"] if present.get("water") else [],
        "choline": ["CHO"] if present.get("choline") else [],
        "urea": ["URE"] if present.get("urea") else [],
        "sorbitol": ["SOR"] if present.get("sorbitol") else [],
        "reline_chloride": ["RCL"] if reline_chloride_count > 0 else [],
        "ions": [
            *(["ICL", "INA"] if present.get("ions") else []),
            *[
                name[:3]
                for name in present.get("ions", [])
                if name.upper() not in {"CL", "CL-", "NA", "NA+"}
            ],
        ],
        "solute": sorted({name[:3] for name in present.get("solute", [])}),
    }
    legend = _dynamic_trajectory_legend_html(
        present,
        show_protein=show_protein,
        show_solute=show_solute,
        show_solvent=show_solvent,
        show_ions=show_ions,
    )
    legend_block = f'<div class="viewer-legend">{legend}</div>' if legend else ""

    script = f"""
    <div class="viewer-shell short-md">
      <div id="viewer" class="viewer"></div>
      <div class="viewer-badge">Production trajectory</div>
      {legend_block}
    </div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor: "{mv.TRAJ_BG}"}});
      viewer.addModelsAsFrames({json.dumps(pdb)}, "pdb");
      viewer.setStyle({{}}, {{sphere: {{hidden: true}}}});
      const components = {json.dumps(components_map)};
      if ({json.dumps(bool(show_solvent))} && components.water.length) {{
        viewer.setStyle({{resn: components.water}}, {{sphere: {{radius: 0.462, color: "#B0BEC5", opacity: 0.62}}}});
      }}
      if ({json.dumps(bool(show_solvent))} && components.choline.length) {{
        viewer.setStyle({{resn: components.choline}}, {{sphere: {{radius: 0.74, color: "#7E57C2"}}}});
      }}
      if ({json.dumps(bool(show_solvent))} && components.urea.length) {{
        viewer.setStyle({{resn: components.urea}}, {{sphere: {{radius: 0.70, color: "#FFB74D"}}}});
      }}
      if ({json.dumps(bool(show_solvent))} && components.sorbitol.length) {{
        viewer.setStyle({{resn: components.sorbitol}}, {{sphere: {{radius: 0.72, color: "#FF7043"}}}});
      }}
      if ({json.dumps(bool(show_solvent))} && components.reline_chloride.length) {{
        viewer.setStyle({{resn: components.reline_chloride}}, {{sphere: {{radius: 0.56, color: "#EC407A"}}}});
      }}
      if ({json.dumps(bool(show_ions))} && components.ions.length) {{
        viewer.setStyle({{resn: components.ions}}, {{sphere: {{radius: 0.605, color: "limegreen"}}}});
      }}
      if ({json.dumps(bool(show_solute))} && components.solute.length) {{
        viewer.setStyle({{resn: components.solute}}, {{sphere: {{radius: 0.836, color: "yellow"}}}});
      }}
      if ({json.dumps(bool(show_protein))}) {{
        viewer.setStyle({{atom: "BB"}}, {{sphere: {{radius: 0.902, color: "{mv.BLUE}"}}}});
        viewer.setStyle({{atom: /^BB\\d+$/}}, {{sphere: {{radius: 0.902, color: "{mv.BLUE}"}}}});
        viewer.setStyle({{atom: /^SC/}}, {{sphere: {{radius: 0.902, color: "{mv.BLUE_LIGHT}"}}}});
      }}
      viewer.zoomTo();
      viewer.animate({{loop: "forward", reps: 0}});
      viewer.render();
    </script>
    <style>
      html, body {{margin:0;padding:0;overflow:hidden;background:{mv.TRAJ_BG};}}
      .viewer-shell {{
        position:relative;width:100%;height:{height}px;box-sizing:border-box;
        border:1px solid rgba(66,199,213,.28);border-radius:16px;
        overflow:hidden;background:{mv.TRAJ_BG};box-shadow:inset 0 0 46px rgba(66,199,213,.08);
      }}
      .viewer {{width:100%;height:{height}px;overflow:hidden;}}
      .viewer-badge {{
        position:absolute;left:14px;bottom:14px;padding:8px 10px;
        border:1px solid rgba(66,199,213,.35);border-radius:999px;
        background:rgba(14,13,17,.78);color:#F5FAFA;
        font:700 12px/1.2 sans-serif;pointer-events:none;
      }}
      .viewer-legend {{
        position:absolute;right:14px;bottom:14px;display:flex;gap:10px;flex-wrap:wrap;
        max-width:70%;padding:7px 10px;border:1px solid rgba(66,199,213,.22);
        border-radius:12px;background:rgba(14,13,17,.78);color:#F5FAFA;
        font:600 11px/1.2 sans-serif;pointer-events:none;
      }}
      .viewer-legend span {{display:flex;align-items:center;gap:5px;}}
      .viewer-legend i {{width:9px;height:9px;border-radius:50%;display:inline-block;}}
    </style>
    """
    mv._components_html(script, height=height + 2)
    return frames


def apply() -> None:
    mv.render_trajectory = render_trajectory
