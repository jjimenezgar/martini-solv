from pathlib import Path

from martini_solv import molecular_viewer as mv


def test_go_virtual_ca_sites_are_hidden_from_trajectory_component_map(tmp_path: Path):
    gro = tmp_path / "system.gro"
    gro.write_text(
        "test\n"
        "4\n"
        "    1ALA     BB    1   0.100   0.100   0.100\n"
        "    1ALA    SC1    2   0.200   0.100   0.100\n"
        "    1ALA     CA    3   0.100   0.100   0.100\n"
        "    2W        W    4   0.500   0.500   0.500\n"
        "   2.00000   2.00000   2.00000\n"
    )

    components = mv._gif_component_map(gro)

    assert components[:3] == ["protein", "protein", "hidden"]
    assert components[3] == "solvent"
