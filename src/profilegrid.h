#ifndef PROFILEGRID_H
#define PROFILEGRID_H

#include "hull.h"

// Auxiliary calls for managing the uniform calculation grid
class ProfileGrid {
private:

    std::vector<double> beams, slopes, slopes0;
    std::vector<glm::dvec3> areas;
    std::shared_ptr<Hull> hull;
    double dx = 0;
    double dz = 0;
    double x_start = 0;
    double z_start = 0;
    std::size_t nx = 0;
    std::size_t nz = 0;

public:

    ProfileGrid(std::shared_ptr<Hull> ship_hull) : hull(ship_hull) {}
    ~ProfileGrid() {}

    std::shared_ptr<Hull> get_hull() {
        return hull;
    }

    double get_dx() const {
        return dx;
    }

    double get_dz() const {
        return dz;
    }

    std::size_t size_x() const {
        return nx;
    }

    std::size_t size_z() const {
        return nz;
    }

    std::size_t index(std::size_t xi, std::size_t zi) const {
        return zi * size_x() + xi;
    }

    std::vector<double> get_xs() const {
        std::vector<double> ret;
        for (std::size_t xi = 0; xi < size_x(); xi++) ret.push_back(x_start + xi * dx);
        return ret;
    }

    std::vector<double> get_zs() const {
        std::vector<double> ret;
        for (std::size_t zi = 0; zi < size_z(); zi++) ret.push_back(z_start + zi * dz);
        return ret;
    }

    std::size_t get_x_transom(std::size_t zi) const {

        double const eps = hull->get_beam_wl() * 0.1; // everything wider than this is still transom

        for (std::size_t xi = 0; xi < size_x(); xi++) { // TODO allows only +X for now

            double const y = get_beam(xi, zi);

            if (y > eps) {
                return xi;
            } else if (y > 0.0) {
                return nx;
            }
        }

        return nx;
    }

    bool has_transom() const {
        return get_beam_transom(0) > 1e-2;
    }

    double get_beam_transom(std::size_t zi) const {

        auto const xi = get_x_transom(zi);

        if (xi >= size_x()) {   // get_x_transom returns nx (an x-index) as the "no transom" sentinel
            return 0;
        }

        return get_beam(xi, zi);
    }

    double get_draft_transom_approx() const {

        double const eps = hull->get_beam_wl() * 0.2; // everything wider than this is still transom
        for (std::size_t zi = 0; zi < size_z(); zi++) {
            if (get_beam_transom(zi) < eps) {
                return std::abs(get_z(0) - get_z(zi));
            }
        }

        return 0;
    }

    double get_x(std::size_t xi) const {
        return x_start + xi * dx;
    }

    double get_z(std::size_t zi) const {
        return z_start + zi * dz;
    }

    double get_beam(std::size_t xi, std::size_t zi) const {
        return beams[index(xi, zi)];
    }

    double get_slope(std::size_t xi, std::size_t zi) const {
        return slopes[index(xi, zi)];
    }

    double get_slope0(std::size_t xi, std::size_t zi) const {
        return slopes0[index(xi, zi)];
    }

    const std::vector<double>& get_beams() const { return beams; }
    const std::vector<double>& get_slopes() const { return slopes; }
    const std::vector<double>& get_slopes0() const { return slopes0; }

    glm::dvec3 get_area(std::size_t xi, std::size_t zi) const {
        return areas[index(xi, zi)];
    }

    void set_slope(std::size_t xi, std::size_t zi, double new_slope) {
        slopes[index(xi, zi)] = new_slope;
    }

    void setup(std::size_t n_x, std::size_t n_z) {

        nx = n_x;
        nz = n_z;

        // x capture whole hull
        double const eps_x = 0.02 * hull->get_size(true).x;
        x_start = hull->get_min(true).x - eps_x;
        dx = (hull->get_max(true).x + eps_x - x_start) / nx;

        // z from waterline to bottom
        double const eps_z = 0.02 * (hull->get_waterline() - hull->get_min(true).z);
        z_start = hull->get_waterline();
        dz = (hull->get_min(true).z - eps_z - z_start) / nz;

        // prepare heightfield. Use assign (not resize): map_mesh_to_grid only writes the
        // cells a hull triangle covers, so cells left un-hit must be zeroed every rebuild —
        // resize keeps stale values when the size is unchanged, which made the grid (and
        // hence the wave/Noblesse forces) depend on the attitude it was FIRST built at.
        auto const N = nx * nz;
        beams.assign(N, 0);
        slopes.assign(N, 0);
        slopes0.assign(N, 0);
        areas.assign(N, {0, 0, 0});

        // construct heightfield
        map_mesh_to_grid();
    }

    bool is_valid() const {
        return nx >= 3 && nz >= 3;
    }

private:

    void map_mesh_to_grid() {

        auto mesh = hull->get_mesh();
        if (!mesh) {
            std::cout << "Cannot make the calculation grid without the mesh" << std::endl;
            return;
        }

//        std::cout << "Making hydrodynamic calculation grid" << std::endl;

        auto verts = mesh->get_vertices();
        auto faces = mesh->get_faces();
        double const hx = get_dx();
        double const hz = std::abs(get_dz());
        double const eps_x = 1e-2 * hx;
        double const eps_z = 1e-2 * hz;

        glm::dvec3 y_dir(0, hull->get_side(), 0);

        // find beams
        for (std::size_t fi = 0; fi < faces.size(); fi++) {

            auto const f = faces[fi];
            glm::dvec3 tv[3] = {verts[f.x], verts[f.y], verts[f.z]};
            glm::dvec3 tv_alt[3] = {verts[f.z], verts[f.y], verts[f.x]};

            // xz bbox
            int start_x = int((std::min({tv[0].x, tv[1].x, tv[2].x}) - x_start - eps_x) / hx);
            int end_x = int((std::max({tv[0].x, tv[1].x, tv[2].x}) - x_start + eps_x) / hx);
            int start_z = int((z_start - std::max({tv[0].z, tv[1].z, tv[2].z}) - eps_z) / hz);
            int end_z = int((z_start - std::min({tv[0].z, tv[1].z, tv[2].z}) + eps_z) / hz);

            for (int xi = start_x+1; xi <= end_x; xi++) {

                if (xi < 0 || xi >= int(size_x())) {
                    continue;
                }

                for (int zi = start_z+1; zi <= end_z; zi++) {

                    if (zi < 0 || zi >= int(size_z())) {
                        continue;
                    }

                    auto const id = index(xi, zi);
                    auto from = glm::dvec3(get_x(xi), 0.0, get_z(zi));
                    double dist;
                    if (intersect_ray_triangle(from, y_dir, tv, dist)) {
                        beams[id] = std::abs(dist);
                    } else if (intersect_ray_triangle(from, y_dir, tv_alt, dist)) {
                        beams[id] = std::abs(dist);
                    }
                }
            }
        } // beams

        // calculate slopes from beams
        for (std::size_t xi = 1; xi < size_x()-1; xi++) { // skip boundary nodes
            for (std::size_t zi = 0; zi < size_z(); zi++) {

                auto const id = index(xi, zi);
                double const beam_xz = beams[id];

                if (beam_xz == 0.0) {
                    continue;
                }

                std::size_t const x_bck = get_beam(xi - 1, zi) == 0.0 ? xi : xi - 1;
                std::size_t const x_fwd = get_beam(xi + 1, zi) == 0.0 ? xi : xi + 1;

                if (x_bck != x_fwd) {

                    slopes0[id] = slopes[id] = (-hull->get_fwd()) * (get_beam(x_fwd, zi) - get_beam(x_bck, zi)) / (get_x(x_fwd) - get_x(x_bck));

                    // get normal from the grid, not from the mesh
                    std::size_t const z_up = std::max(0, int(zi) - 1);
                    std::size_t const z_down = std::min(size_z() - 1, zi + 1);

                    glm::dvec3 const edge1(get_dx(), (get_beam(x_fwd, zi) - get_beam(x_bck, zi)) / (x_fwd - x_bck), 0.0);
                    glm::dvec3 const edge2(0.0, (get_beam(xi, z_up) - get_beam(xi, z_down)) * 0.5, get_dz());

                    areas[id] = glm::cross(edge1, edge2);

                    if (get_beam(x_fwd, zi) == 0.0 || get_beam(x_bck, zi) == 0.0 || get_beam(xi, z_up) == 0.0 || get_beam(xi, z_down) == 0.0)
                        areas[id].x = 0.0;

                }

            }
        }

//        std::cout << "Converted hull to XZ grid of size " << size_x() << " x " << size_z() << std::endl;
        return;
    }

};

#endif // PROFILEGRID_H
