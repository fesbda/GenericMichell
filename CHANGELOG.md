# Changelog

## 1.1.0 — 2026-09-07

Solver state behind the revised Ocean Engineering manuscript. Two of the four items
below entered on 20 July 2026, with the submitted version, and are released publicly
for the first time here; the other two come from the revision.

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
