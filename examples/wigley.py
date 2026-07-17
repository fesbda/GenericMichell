"""Procedural half-hull mesh for a classical Wigley form."""


def make_wigley_hull(
    michell,
    length=4.0,
    beam=0.4,
    draft=0.25,
    freeboard=0.25,
    nx=81,
    nz=33,
):
    """Construct a triangulated +y Wigley half hull and float it at z=0."""
    if nx < 3 or nz < 3:
        raise ValueError("nx and nz must both be at least 3")

    mesh = michell.Mesh()
    for ix in range(nx):
        x = length * ix / (nx - 1)
        longitudinal = 1.0 - (2.0 * x / length - 1.0) ** 2
        for iz in range(nz):
            z = -draft + (draft + freeboard) * iz / (nz - 1)
            scale = draft if z <= 0.0 else freeboard
            vertical = max(0.0, 1.0 - (z / scale) ** 2)
            y = 0.5 * beam * longitudinal * vertical
            mesh.add_vertex(x, y, z)

    for ix in range(nx - 1):
        for iz in range(nz - 1):
            a = ix * nz + iz
            b = (ix + 1) * nz + iz
            c = b + 1
            d = a + 1
            mesh.add_face(a, b, c)
            mesh.add_face(a, c, d)
    mesh.update()

    environment = michell.Environment()
    environment.density = 999.3
    hull = michell.Hull(mesh, environment, 1)
    hull.set_waterline(0.0)
    hull.set_cg(hull.get_cb().x, 0.0)
    hull.capture_reference_displacement()
    return hull
