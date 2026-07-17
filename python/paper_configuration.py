"""Fixed solver configuration used in the associated manuscript."""


def make_paper_solver(michell, hull):
    """Return a prepared Resistance instance with the manuscript settings."""
    solver = michell.Resistance(hull)
    solver.set_grid_precision(michell.ResistanceGridPrecision.Fine)
    solver.set_transom_method(
        michell.TransomMethod.DoctorsDynamic,
        michell.TransomMethod.DoctorsDynamic,
    )
    solver.set_params(-0.8, -2.0, 2.0, 3.8, 0.0)
    solver.set_pdyn_projection(True, 1.6, 0.0, 2.5, 4.5)
    solver.set_whisker_warp_only(True)
    solver.set_friction_width(True)
    solver.set_friction_cv_ramp(99.0, 100.0)
    solver.set_savitsky_trim_weight(0.3)
    solver.prepare()
    return solver
