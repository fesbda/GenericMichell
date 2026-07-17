#pragma once

#include "profilegrid.h"
#include "streamline.h"

/// Computation of local flow
class LocalFlow {

private:
    std::shared_ptr<Hull> hull;
    std::shared_ptr<Environment> env;
    std::unique_ptr<ProfileGrid> grid;
    PchipInterpolator deflection;

    double param1 = -0.8; // bow coef
    double param2 = -2; // fore coef
    double param3 = 1.0; // aft coef

    bool prepared = false;

    std::vector<double> deflection_field;

    // buttocks
    std::vector<double> ys;
    std::vector<Streamline> streamlines;

public:

    LocalFlow(std::shared_ptr<Hull> ship_hull)
    : hull(ship_hull)
    {
        env = hull->get_environment();
        grid = std::make_unique<ProfileGrid>(hull);
        deflection = PchipInterpolator(
            {0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.08, 1.28, 1.47, 1.57, 1.67, 1.87, 2.06, 2.26, 2.45, 2.65, 2.85, 3.04, 3.142},
            {0.0, 2.58787224e-02,    -2.07058187e-02,     -4.11867727e-02,    -5.41260291e-02,   -4.65911214e-02,    -2.29935673e-02,    -5.64612831e-03,    -7.46950397e-04,     0.0,            -2.56595981e-02,    -1.69933303e-01,    -4.37234271e-01,    -5.45162121e-01,    -6.96964561e-01,    -7.83925538e-01,    -9.04849779e-01,  -1.11388172,        -1.35}
        );
    }

    ~LocalFlow() {}

    ProfileGrid* get_grid() const {
        return grid.get();
    }

    void prepare(int n_slices) {

        // use minimum 3 buttocks as streamlines
        n_slices = std::max(3, n_slices);

        double const y_max = 0.97 * (hull->get_beam_wl() * 0.5);
        ys = linspace(1e-3, y_max, n_slices);

        streamlines.resize(n_slices);

        for (std::size_t i = 0; i < ys.size(); i++) {

            // get buttock geometry
            std::vector<double> xs;
            std::vector<double> zs;
            hull->slice_buttock(ys[i], xs, zs, true);

            // assume that the buttock is a streamline
            streamlines[i].setup_one_side(xs, zs, hull->get_fwd());
        }

        prepared = true;
    }

    void set_params(double p1, double p2, double /*p3*/, double /*p4*/, double /*p5*/) {
        param1 = p1;
        param2 = p2;
    }

    void set_deflection_table(const std::vector<double>& xs, const std::vector<double>& ys) {
        deflection = PchipInterpolator(xs, ys);
    }

    void set_deflection_field(const std::vector<double>& field) {
        deflection_field = field;
    }

    void set_grid(std::size_t nx, std::size_t nz) {

        grid->setup(nx, nz);

        bool use_field = !deflection_field.empty() && deflection_field.size() == grid->size_x() * grid->size_z();

        for (std::size_t xi = 0; xi < grid->size_x(); xi++) {
            for (std::size_t zi = 0; zi < grid->size_z(); zi++) {
                double const old_slope = grid->get_slope0(xi, zi);
                double new_slope;
                
                if (use_field) {
                    double delta_beta = deflection_field[grid->index(xi, zi)];
                    new_slope = std::max(std::tan(std::atan(old_slope * hull->get_fwd()) + delta_beta), -0.55);
                } else {
                    new_slope = get_bl_slope(old_slope);
                }
                
                grid->set_slope(xi, zi, new_slope);
            }
        }

    }

    glm::dvec3 get_drag_lift_torque_noblesse(double speed, std::vector<double>* out_pressure_field = nullptr) {

        double const F = hull->get_fn(speed);
        auto const cg = hull->get_cg(true);
//        double const dS = std::abs(grid->get_dx() * grid->get_dz());

        // output
        double drag = 0;
        double lift = 0;
        double torque = 0;
//        double p_min = 1e8;
//        double p_max = -1e8;
        double sum_area = 0;

        if (out_pressure_field) {
            // assign (not resize): NaN/skipped cells below are left unwritten, so a reused buffer of
            // the same size would keep stale values from a previous call. Zero the whole field.
            out_pressure_field->assign(grid->size_x() * grid->size_z(), 0.0);
        }

        #pragma omp parallel for reduction(+:sum_area,drag,lift,torque) collapse(2) //reduction(max:p_max) reduction(min:p_min)
        for (std::size_t xi = 0; xi < grid->size_x(); xi++) {
            for (std::size_t zi = 0; zi < grid->size_z(); zi++) {

//                if (grid->get_beam(xi,zi) < 1e-4) {
//                    continue;
//                }

                auto [p, p0] = get_thin_local_pressure(speed, F, xi, zi);

                if (std::isnan(p) || std::isnan(p0)) {
                    continue; // avoid NaN
                }

                if (out_pressure_field) {
                    (*out_pressure_field)[grid->index(xi, zi)] = p;
                }

//                p_min = std::min(p, p_min);
//                p_max = std::max(p, p_max);

                auto const dA = grid->get_area(xi, zi);

                drag += p * dA.x;
                lift += p * dA.z;
                // bow-up pitch moment = Δx·Fz − Δz·Fx (Fz=p·dA.z, Fx=p·dA.x); the Δz·Fx
                // (horizontal-force × vertical-lever) term was + (wrong rotation sense).
                torque += p * (dA.z * (grid->get_x(xi) - cg.x) - dA.x * (grid->get_z(zi) - cg.z));

                sum_area += std::abs(dA.z);
            }
        }

        // account for two hull sides
        drag *= 2;
        lift *= 2;// * (0.5*hull->get_waterplane_area() / (sum_area + 1e-5));
        torque *= 2;

//        std::cout << "for speed " << speed << " local drag is " << drag << " N, lift is " << lift/9.81 << "kg, and torque is " << torque/9.81 << " kg.m" << std::endl;
//        std::cout << "pmin " << p_min << ", pmax " << p_max  << ", pinf " << 0.5 * env->get_density() * sq(speed) << std::endl;
//        std::cout << "Area_WL " << hull->get_waterplane_area() << " vs " << 2*sum_area << " ratio " << (0.5*hull->get_waterplane_area() / sum_area) << std::endl;

        return {std::abs(drag), lift, torque};
    }


    /// Local-flow centerline elevation E_L(x, 0, 0).
    /// Same Noblesse simplified Green function as get_thin_local_pressure but with
    /// flo at the free surface (z = wl, so nondim flo.z = 0). At z_flo = 0 we have
    /// r == r1 exactly, so the 1/r^3 - 1/r1^3 term in M (paper Eq. 30b) vanishes;
    /// only the F^2 group survives. Result is in metres (pressure / rho*g, doubled
    /// for port + starboard).
    std::vector<double> get_local_elevation_centerline(double speed, const std::vector<double>& xs_eval) {

        std::vector<double> result(xs_eval.size(), 0.0);

        if (speed <= 0 || !grid->is_valid()) {
            return result;
        }

        double const F = hull->get_fn(speed);
        if (F < 1e-2) {
            return result;
        }

        auto const linv = 1.0 / std::abs(grid->get_x(0) - grid->get_x(grid->size_x()));
        auto const wl = grid->get_z(0);
        auto const hp = std::abs(grid->get_dx()) * linv;
        auto const hp3 = cub(hp);
        auto const hp4 = sq(sq(hp));
        auto const hp6 = sq(hp3);
        auto const F2 = sq(F);

        double const dx = std::abs(grid->get_dx());
        double const dz = std::abs(grid->get_dz());
        double const dS = dx * dz * sq(linv);

        double const rho = env->get_density();
        double const rho_g = rho * env->get_gravity();

        std::size_t const Nx = grid->size_x();
        std::size_t const Nz = grid->size_z();
        int const fwd = hull->get_fwd();

        #pragma omp parallel for
        for (std::size_t e = 0; e < xs_eval.size(); e++) {

            // flo at centerplane y=0 and free surface z=wl, in nondim coords.
            glm::dvec3 const flo(xs_eval[e] * linv, 0.0, 0.0);

            double sum = 0.0;

            for (std::size_t xi = 0; xi < Nx; xi++) {
                for (std::size_t zi = 0; zi < Nz; zi++) {

                    glm::dvec3 const src(grid->get_x(xi) * linv, 0.0, (grid->get_z(zi) - wl) * linv);

                    double const dxi = flo.x - src.x;
                    double const r1 = std::sqrt(sq(dxi) + sq(src.z));  // = r at z_flo=0
                    if (r1 > 3.0) {
                        continue;
                    }
                    if (r1 < 1e-9) { continue; }   // P5: skip the singular self source term (psi=-0/0)

                    double const psi = -src.z / (r1 + std::abs(dxi));  // src.z<=0, so psi>=0
                    double const sigma = fwd * grid->get_slope(xi, zi);
                    // Smoothed sign function: as x_eval sweeps the hull, the bare sgn(Δx)
                    // flips at every source position, producing per-grid-cell ripple in the
                    // discrete sum. Replacing sgn with Δx/√(Δx²+δ²) spreads each flip over
                    // ~one grid cell. The continuous integral is unchanged in the δ→0 limit.
                    // dxi = x' - x, and sign(x - x') multiplies the WHOLE bracket (Eq. 30a),
                    // so the smoothed sign is -dxi/|dxi|.
                    double const delta = hp;
                    double const sgn = -dxi / std::sqrt(sq(dxi) + sq(delta));   // sign(x - x')

                    double const eps_ = hp6 / (hp4 + sq(sq(r1)));
                    double const M = 2.0 / cub(F2 + r1) *
                        (1.0 + F2 / std::sqrt(sq(r1) + eps_) *
                            (1.0 + 2.0*psi - 2.3*sq(psi) * (F2 - 2.0*r1) / (F2 + r1)));

                    double const second = 2.0 * psi * F2 / cub(F2 + r1) *
                        (1.0 + 4.6*psi + F2 / std::sqrt(sq(r1) + eps_));

                    sum += sigma * sgn * (std::abs(dxi) * M + second);
                }
            }

            sum *= dS / (2.0 * M_PI) * rho * sq(speed);  // single-side pressure (Pa)
            if (std::isnan(sum)) {
                sum = 0.0;
            }
            result[e] = 2.0 * sum / rho_g;  // both sides, pressure → elevation (m)
        }

        return result;
    }


    /// Near-field axial-deficit disturbance velocity at world point (X, Y, Z),
    /// from the Noblesse simplified Green function. Only the axial component is
    /// returned (linearized Bernoulli u_x = -p/(rho*V)); the transverse and
    /// vertical disturbance components are well-modelled by the ThinShip Rankine
    /// + image gradient already used in ThinShip::get_velocity, and re-deriving
    /// them from the Noblesse pressure-only kernel would over-count.
    glm::dvec3 get_velocity_local(double speed, double X, double Y, double Z) const {

        if (speed <= 0 || !grid->is_valid()) {
            return {0.0, 0.0, 0.0};
        }

        double const p = eval_local_pressure_at(speed, X, Y, Z);
        double const denom = env->get_density() * speed;
        double const u_x = (denom > 0.0 && !std::isnan(p)) ? (-p / denom) : 0.0;

        return {u_x, 0.0, 0.0};
    }

    glm::dvec3 get_drag_lift_torque_xflow(double speed) {

        if (speed <= 1e-2 || hull->get_fn(speed) < 1e-2 || !prepared) {
            return {0, 0, 0};
        }

        auto const cg = hull->get_cg(true);
        double drag = 0, lift = 0, torque = 0;

        for (std::size_t i = 0; i < ys.size(); i++) {

            streamlines[i].calculate(hull->get_environment(), speed);
            auto const D = streamlines[i].get_drag(2, false);
            auto const L = streamlines[i].get_lift(2, false);

            // simplest form of integration
            drag += D.first;
            lift += L.first;
            torque += L.first * (L.second - cg.x)
                    + D.first * (cg.z - D.second); // cos here drag is positive in direction of movement
        }

        double const dy = 2 * std::abs(ys[1] - ys[0]);
        return {drag*dy, lift*dy, torque*dy};
    }

private:

    double get_bl_slope(double geometry_slope) const {

        if (std::abs(geometry_slope) < 1e-16) {
            return 0.0;
        }

        double const slope_angle = std::atan(geometry_slope * hull->get_fwd());
        double delta_beta = deflection(slope_angle);

        if (true)
        {
            double const multiplier = slope_angle < 0.0 ? param1 : param2;
            delta_beta = deflection(-slope_angle) * multiplier;
        }
//        else
//        {

//            double ang = M_PI_2 + std::atan(geometry_slope * hull->get_fwd());
//            if (ang < M_PI_2)
//                delta_beta = (-1.6211*sq(ang) + 2.5465*ang) * param1;
//            else
//                delta_beta = sq(ang - M_PI_2) * param2;
//        }

        double const new_slope = std::max(std::tan(slope_angle + delta_beta), -0.55);
        return new_slope;

//        double const tan_beta = geometry_slope;
//        double const beta_from_slope = M_PI_2 - std::atan(tan_beta);
//        double const delta_beta = deflection(beta_from_slope);
//        double const multiplier = delta_beta > 0.0 ? param1 : (beta_from_slope < M_PI_2 ? param2 : param3);
//        double const new_slope = std::max(std::tan(M_PI_2 - beta_from_slope - multiplier * delta_beta), -0.5);
//        return new_slope;
    }

    /// Evaluate the Noblesse simplified-Green-function pressure at an arbitrary
    /// world point (X, Y, Z). Mirrors get_thin_local_pressure but takes a world
    /// point rather than a hull grid index, and returns the BL-corrected pressure
    /// only (single-side, doubled for port+starboard symmetry, in Pa).
    double eval_local_pressure_at(double speed, double X, double Y, double Z) const {

        if (speed <= 0 || !grid->is_valid()) {
            return 0.0;
        }

        double const F = hull->get_fn(speed);
        if (F < 1e-2) {
            return 0.0;
        }

        double const wl = grid->get_z(0);
        double const linv = 1.0 / std::abs(grid->get_x(0) - grid->get_x(grid->size_x()));
        double const hp = std::abs(grid->get_dx()) * linv;
        double const hp3 = cub(hp);
        double const hp4 = sq(sq(hp));
        double const hp6 = sq(hp3);
        double const F2 = sq(F);

        double const dx = std::abs(grid->get_dx());
        double const dz = std::abs(grid->get_dz());
        double const dS = dx * dz * sq(linv);

        glm::dvec3 const flo(X * linv, Y * linv, (Z - wl) * linv);

        double sum = 0.0;

        #pragma omp parallel for reduction(+:sum) collapse(2)
        for (std::size_t xi = 0; xi < grid->size_x(); xi++) {
            for (std::size_t zi = 0; zi < grid->size_z(); zi++) {

                glm::dvec3 const src(grid->get_x(xi) * linv, 0.0, (grid->get_z(zi) - wl) * linv);

                double const dxi = flo.x - src.x;
                double const dyi = flo.y - src.y;
                double const r = std::sqrt(sq(dxi) + sq(dyi) + sq(flo.z - src.z));
                if (r > 3.0) continue;

                double const r1 = std::sqrt(sq(dxi) + sq(dyi) + sq(flo.z + src.z));
                if (r1 < 1e-9) { continue; }   // P5: skip the singular self/mirror source term (psi=-0/0)
                double const psi = -(flo.z + src.z) / (r1 + std::abs(dxi));
                double const sigma = hull->get_fwd() * grid->get_slope(xi, zi);
                // Smoothed sign(x - x') for numerical stability on/near the centerplane.
                // dxi = x' - x (flo - src), so sign(x - x') = -dxi/|dxi|. Multiplies the WHOLE
                // bracket per Eq. (30a) (kernel is odd in x - x').
                double const delta = hp;
                double const sgn = -dxi / std::sqrt(sq(dxi) + sq(delta));   // sign(x - x')

                double const eps  = hp6 / (hp3 + cub(r));
                double const eps1 = hp6 / (hp3 + cub(r1));
                double const eps_ = hp6 / (hp4 + sq(sq(r1)));
                double const M = 1.0 / (cub(r) + eps) - 1.0 / (cub(r1) + eps1)
                    + 2.0 / cub(F2 + r1) * (1.0 + F2 / std::sqrt(sq(r1) + eps_)
                        * (1.0 + 2.0 * psi - 2.3 * sq(psi) * (F2 - 2.0 * r1) / (F2 + r1)));

                double const second_term = 2.0 * psi * F2 / cub(F2 + r1)
                    * (1.0 + 4.6 * psi + F2 / std::sqrt(sq(r1) + eps_));

                sum += sigma * sgn * (std::abs(dxi) * M + second_term);
            }
        }

        sum *= dS / (2.0 * M_PI) * (env->get_density() * sq(speed));
        if (std::isnan(sum)) sum = 0.0;
        // Single-side pressure (Pa), matching the convention of get_thin_local_pressure
        // and the BiLinearInterpolator field returned by Resistance::get_local_pressure_field.
        return sum;
    }

    std::pair<double, double> get_thin_local_pressure(double speed, double F, std::size_t floxi, std::size_t flozi) {

        auto const linv = 1.0 / std::abs(grid->get_x(0) - grid->get_x(grid->size_x())); // hull->get_length_wl();
        auto const wl = grid->get_z(0);
        // Noblesse Sec. 9 recommends regularizer h_p on the order of the panel size sigma.
        auto const hp = std::abs(grid->get_dx()) * linv;
        auto const hp2 = sq(hp);
        auto const hp3 = cub(hp);
        auto const hp4 = sq(hp2);
        auto const hp6 = sq(hp3);
        auto const F2 = sq(F);

        auto dx = std::abs(grid->get_dx());
        auto dz = std::abs(grid->get_dz());
        // Noblesse simplified Green function (Eq. 50) is derived for y - y' = 0;
        // rely on the eps/eps1/eps_ regularizers for the on-source singularity.
        auto const flo = glm::dvec3(grid->get_x(floxi), 0.0, grid->get_z(flozi) - wl) * linv;
        auto const dS = dx * dz * sq(linv);
        double sum = 0;
        double sum0 = 0;

//        #pragma omp parallel for reduction(+:sum) collapse(2)
        for (std::size_t xi = 0; xi < grid->size_x(); xi++) {
            for (std::size_t zi = 0; zi < grid->size_z(); zi++) {

//                if (grid->get_beam(xi,zi) < 1e-4) {
//                    continue;
//                }

                auto const src = glm::dvec3(grid->get_x(xi), 0.0, grid->get_z(zi) - wl) * linv;
                double const r = std::sqrt(sq(flo.x - src.x) + sq(flo.y - src.y) + sq(flo.z - src.z));

                if (r > 3.0) {
                    continue;
                }

                double const r1 = std::sqrt(sq(flo.x - src.x) + sq(flo.y - src.y) + sq(flo.z + src.z)); // differs by z addition
                // P5: at a waterline-row eval cell (flo.z=0) the coincident waterline source (src.z=0,
                // same x) gives r1=0 -> psi=-0/0=NaN, which used to poison the whole cell (dropping the
                // top row, which carries the largest dA.x). Skip ONLY this singular source term.
                if (r1 < 1e-9) { continue; }
                double const psi = -(flo.z + src.z) / (r1 + std::abs(flo.x - src.x));
                // Paper Eq. (30a): sign(x - x') multiplies the WHOLE bracket (both the |x-x'|*M
                // term and the second term). The kernel equals d g^SL/dx (Eq. 23a), which is odd
                // in (x - x'); applying the sign to the second term only left the M term with the
                // wrong (even) parity. Verified against Eq. (23a) by finite difference.
                double const sigma0 = hull->get_fwd() * grid->get_slope0(xi,zi);
                double const sigma  = hull->get_fwd() * grid->get_slope(xi,zi);
                double const sgn = sign(src.x - flo.x);   // sign(x - x')

                double const eps = hp6 / (hp3 + cub(r));
                double const eps1 = hp6 / (hp3 + cub(r1));
                double const eps_ = hp6 / (hp4 + sq(sq(r1)));
                double const M = 1 / (cub(r) + eps) - 1 / (cub(r1) + eps1) + 2 / cub(F2 + r1) * (1 + F2 / std::sqrt(sq(r1) + eps_) * (1 + 2 * psi - 2.3 * sq(psi) * (F2 - 2 * r1) / (F2 + r1)));

                double const second_term = 2.0 * psi * F2 / cub(F2 + r1) * (1 + 4.6 * psi + F2 / std::sqrt(sq(r1) + eps_));
                sum  += sigma  * sgn * (std::abs(src.x - flo.x) * M + second_term);
                sum0 += sigma0 * sgn * (std::abs(src.x - flo.x) * M + second_term);
            }
        }

        sum *= dS / (2*M_PI) * (env->get_density() * sq(speed));
        sum0 *= dS / (2*M_PI) * (env->get_density() * sq(speed));

        return std::make_pair(sum, sum0);
    }

};
