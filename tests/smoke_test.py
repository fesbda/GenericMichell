"""Build-level smoke test using only the Python standard library."""

import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples"))
sys.path.insert(0, str(ROOT / "python"))

import michell  # noqa: E402
from wigley import make_wigley_hull  # noqa: E402
from paper_configuration import make_paper_solver  # noqa: E402


def main():
    hull = make_wigley_hull(michell, nx=41, nz=21)
    assert 3.9 < hull.get_length_wl() < 4.1
    assert 0.39 < hull.get_beam_wl() < 0.41
    assert 0.24 < hull.get_draft() < 0.26
    assert hull.get_displacement() > 0.0

    solver = michell.Resistance(hull)
    solver.set_grid_precision(michell.ResistanceGridPrecision.Coarse)
    solver.prepare()
    speed = hull.get_speed(0.30)
    components = (
        solver.get_friction_resistance(speed, michell.ResistanceMethod.Auto),
        solver.get_wave_resistance(speed, michell.ResistanceMethod.Auto),
        solver.get_pressure_resistance(speed, michell.ResistanceMethod.Auto),
        solver.get_total_resistance(speed, michell.ResistanceMethod.Auto),
    )
    assert all(math.isfinite(value) and value >= 0.0 for value in components)
    assert components[-1] > 0.0
    paper_solver = make_paper_solver(michell, hull)
    assert paper_solver.is_ok()
    print("smoke test passed:", ", ".join(f"{value:.6g}" for value in components))


if __name__ == "__main__":
    main()
