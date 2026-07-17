# GenericMichell

GenericMichell is a reduced-order resistance and running-attitude solver for
speed sweeps spanning displacement, semi-planing and planing operation. The
Python module combines a C++17 computational core with pybind11 bindings.

This repository is the compact solver release associated with the manuscript
“Continuous prediction of resistance and running attitude across the
displacement-to-planing transition.” It intentionally contains no proprietary
CFD files, third-party hull geometries or experimental datasets.

## What is included

- displacement and planing resistance components;
- free-sinkage and trim equilibrium;
- hydrostatic and dynamic force diagnostics;
- Python bindings for mesh construction and solver access;
- a procedural Wigley-hull example requiring no external geometry files.

The public API is research software. Treat results outside the validation
domain described in the paper as extrapolations, and independently verify
predictions used for design or safety decisions.

## Build

Requirements are CMake 3.18+, a C++17 compiler, Python 3 with development
headers, and internet access during the first configuration. CMake downloads
the pinned pybind11 2.13.6 and GLM 1.0.1 source archives.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Run the example from the repository root:

```bash
PYTHONPATH=build python examples/wigley_sweep.py
```

The output columns are length Froude number, speed, total resistance, friction
resistance, wave resistance and pressure resistance. Forces are in newtons.

## Minimal use

```python
import michell
from examples.wigley import make_wigley_hull

hull = make_wigley_hull(michell)
solver = michell.Resistance(hull)
solver.prepare()
speed = hull.get_speed(0.30)
resistance = solver.get_total_resistance(speed, michell.ResistanceMethod.Auto)
```

See `examples/wigley.py` for mesh conventions and
`python/paper_configuration.py` for the fixed configuration used by the
associated manuscript.

## Coordinate and input conventions

- +x points forward and +z points upward;
- a half-hull surface is sufficient; pass `1` to `michell.Hull` for the +y side;
- meshes must be triangulated and should cross the requested waterline;
- SI units are used throughout;
- pitch is in radians in the API, with negative pitch representing bow-up trim.

## Reproducibility boundary

The numerical core and the manuscript configuration are provided here. The
validation measurements and some hull surfaces are not redistributed because
their publication rights are separate from the solver license. The repository
records the private development source revision in `SOURCE_COMMIT`.

## License and third-party code

GenericMichell is distributed under GPL-2.0-or-later; see `LICENSE`.
`ext/xfoil/` contains a GPL-2.0-or-later C++ translation of XFoil and retains
the notices of Mark Drela and Andre Deperrois. pybind11 and GLM are downloaded
at build time under their own licenses.

## Citation

Citation metadata are provided in `CITATION.cff`. Replace the unpublished
article entry with its DOI after acceptance, and cite the archived software
release DOI once assigned.
