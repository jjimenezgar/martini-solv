from martini_solv.trajectory_legend import _dynamic_trajectory_legend_html


def test_trajectory_legend_shows_only_present_and_enabled_components():
    present = {
        "protein": ["ALA"],
        "water": ["W"],
        "choline": [],
        "urea": [],
        "sorbitol": [],
        "ions": ["NA"],
        "solute": [],
    }
    html = _dynamic_trajectory_legend_html(
        present,
        show_protein=True,
        show_solute=False,
        show_solvent=True,
        show_ions=False,
    )
    assert "Protein" in html
    assert "Water" in html
    assert "Urea" not in html
    assert "Choline" not in html
    assert "Sorbitol" not in html
    assert "Ions" not in html
    assert "Free molecule" not in html


def test_trajectory_legend_hides_present_component_when_viewer_toggle_is_off():
    present = {
        "protein": ["ALA"],
        "water": [],
        "choline": [],
        "urea": ["UREA"],
        "sorbitol": [],
        "ions": [],
        "solute": ["LIG"],
    }
    html = _dynamic_trajectory_legend_html(
        present,
        show_protein=True,
        show_solute=False,
        show_solvent=False,
        show_ions=False,
    )
    assert "Protein" in html
    assert "Urea" not in html
    assert "Free molecule" not in html
