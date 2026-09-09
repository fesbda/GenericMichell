# Changelog

## 1.2.0 — 2026-09-09

Solver state behind the final revision of the Ocean Engineering manuscript. Both items
are numerical repairs found by auditing the implementation against the manuscript's own
equations; neither is a modelling choice, and no closure was retuned.

- **The sinkage solve is repaired and its result is accepted or refused.** Its Newton
  step took the derivative from the hydrostatic waterplane stiffness alone. That omits
  the speed-dependent wave and bottom-pressure lift and under-estimates
  `|dFz/dh|` by about a factor of two at planing, so each step overshot the root and the
  iterate oscillated about it, the bracket shrinking only microscopically because the
  overshoot always landed inside it. On one comparator point the iterate ping-ponged
  between 38.9 and 41.3 mm for its whole budget while plain bisection on the identical
  residual found the root in twelve evaluations. The solve now takes the secant slope of
  the two most recent evaluations, forces a bisection step whenever the bracket fails to
  halve over two iterations, and returns its pose, torque and residual from one final
  evaluation. It is then tested against the same 0.1%-of-weight force tolerance the
  search uses, and the outcome is reported through the new `get_equilibrium_status`
  rather than returned silently. Over the 119-point comparator set the maximum vertical
  residual falls from 3.88% to 0.67% of weight, no point exceeds 1%, and 105 of 119 are
  accepted. Where the residual is a step function of sinkage — the centre-plane integral
  loses a mesh row at once as the hull rises — the solve returns the smaller-imbalance
  side of the step, which is the attainable minimum there.
- **The pressure-drag law reads the attitude balance's own support ledger.** It
  re-derived `L_P,bal` from the raw Michell lift, while the attitude balance fades that
  lift by `(1 − w)` wherever the sectional water-entry closure replaces the planing lift,
  and the displacement-side local-flow lift `(1 − w) L_loc` was in neither. The two
  supports agreed on prismatic hulls but differed by up to 32% of weight on strongly
  warped ones. `L_P,bal = max(0, W − B_pose − L̂_W − (1 − w) L_loc)` now holds on every
  branch. The new `set_pdyn_ledger_support` / `get_pdyn_ledger_support` default to on;
  the change is inert at a prescribed attitude and wherever the planing weight has
  reached unity, so it moves only warped hulls, and it does not touch predicted attitude
  at all. Pooled planing resistance MAPE over the comparator set is 8.37% and pooled trim
  MAE 0.68 degrees, both unchanged at the manuscript's precision; warped-hull resistance
  MAPE improves from 6.5% to 5.8%.
- **`python/paper_configuration.py` states the ledger support explicitly** rather than
  inheriting it from the compiled default, as it already does for the Blount–Fox
  amplitude.
- **`CITATION.cff` now carries the concept DOI** in its `doi` field, which always
  resolves to the latest archived version, with the 1.1.0 version DOI kept under
  `identifiers`. The previous file pinned a version DOI that went stale on release.

## 1.1.0 — 2026-09-07

Solver state behind the revised Ocean Engineering manuscript. Two of the items below
entered on 20 July 2026 with the submitted version and reach the public repository for
the first time here; the rest come from the revision, and one from verifying this
release on Windows.

- **Blount–Fox hump multiplier withdrawn.** The default amplitude is 0.0 rather than
  0.5, and the hump is carried by the pre-planing residuary shortfall alone. The term
  duplicated a deficit that closure already supplies: a crossed ablation puts the
  interaction at −1.66 percentage points, and pooled planing resistance MAPE is 8.41%
  at zero amplitude against 8.76% at half amplitude over the 119-point comparator set.
  `set_blount_fox_amplitude` and `get_blount_fox_amplitude` are new, so the amplitude
  ladder that justifies the removal stays reproducible. Running attitude is unaffected
  in either setting: the multiplier scales planing friction, and shear carries no
  vertical force or moment, so trim and heave are bit-identical across the change.
- **Lift-partition correction.** Pressure drag now uses the bottom-pressure support
  `L_P,bal = max(0, W − B_pose − L_W)`, the same quantity the attitude balance uses,
  instead of projecting the wave-supported share of lift.
- **Pose-stable transom detection** (from the submitted version). A short look-ahead
  stops a one-column positive fringe aft of a full-width transom from being read as a
  pointed stern, and an immersed transom spanning the deepest grid row is no longer
  reported as zero-draft.
- **Doctors transom closure on one representation** (from the submitted version). Beam
  and draft are both taken from the source grid; mixing a mesh-waterplane breadth with
  a grid draft let the closure toggle under sub-millimetre pose changes.
- **The CMake build defaults to OpenMP off.** OpenMP makes the equilibrium solve
  nondeterministic at the 1e-4 level, so a serial build is the configuration the
  manuscript's numbers were produced with; it also broke the MSVC build outright,
  since MSVC implements OpenMP 2.0 and rejects the unsigned loop indices in the field
  loops. 1.0.0 defaulted it on and therefore did not build on Windows.
  `-DMICHELL_USE_OPENMP=ON` restores the previous behaviour.
- **`python/paper_configuration.py` now sets robust chine extraction**, which the
  manuscript configuration requires and 1.0.0 omitted, and states the Blount–Fox
  amplitude explicitly rather than inheriting it from the compiled default.

## 1.0.0 — 2026-07-17

- First curated public solver release.
- Added a pinned CMake build, procedural Wigley example and smoke test.
- Added the fixed configuration used in the associated manuscript.
- Excluded proprietary CFD artifacts, third-party hull meshes and experimental data.
