# Third-party software

## XFoil C++ translation

Files in `ext/xfoil/` are a C++ translation of XFoil. Their retained source
headers identify the following copyright holders:

- Copyright (C) 2000 Mark Drela
- Copyright (C) 2003 Andre Deperrois

Those files are distributed under GPL version 2 or, at the recipient's option,
any later version. The repository-level `LICENSE` contains the GPL version 2
text.

## Build-time dependencies

CMake downloads the following pinned sources at configuration time; they are
not copied into this repository:

- pybind11 2.13.6, BSD-3-Clause
- GLM 1.0.1, MIT

Consult each dependency's source distribution for its complete notices.
