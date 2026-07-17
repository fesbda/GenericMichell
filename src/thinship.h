#ifndef THINSHIP_H
#define THINSHIP_H

#include "profilegrid.h"
#include "interpolation.h"

enum class DemihullMethod {
    Tuck,
    Yeung,
    Mixed
};

enum class TransomShapeMethod {
    DoctorsStatic,
    DoctorsDynamic,
    Robards,
    Mixed
};

class ThinShip {

private:
    std::shared_ptr<Hull> hull;
    std::shared_ptr<Environment> env;
    std::unique_ptr<ProfileGrid> grid;

    bool fix_phase = false;
    bool fix_bl = false;
    bool fix_depth = false;

    // Cached hull spectrum (P, Q) tabulated on a theta grid for the far-field wave field.
    // P, Q only depend on lam = sec(theta) which is even in theta, so we tabulate theta in [0, theta_max].
    struct FarFieldSpectrumEntry {
        double theta;
        double secT;       // sec(theta)
        double sec2T;      // sec^2(theta)
        double sec4T;      // sec^4(theta)
        double secTanT;    // sec(theta) * tan(theta) = sec^2(theta) * sin(theta)
        double P, Q;
    };
    std::vector<FarFieldSpectrumEntry> ff_spectrum_;
    double ff_spectrum_speed_ = -1.0;
    double ff_spectrum_k_ = -1.0;

    // The spectrum cache above is keyed on speed only. Invalidate it whenever an input the
    // spectrum depends on (grid/pose, phase params, deflection, transom method) changes, so a
    // same-speed re-tabulation cannot return a stale spectrum (the fixed-speed multi-pose
    // diagnostic path; production tabulates once per speed at the converged pose, so unaffected).
    void invalidate_ff_spectrum() { ff_spectrum_.clear(); ff_spectrum_speed_ = -1.0; }

    // Upper bound of the theta-grid shared by get_wave_resistance,
    // tabulate_farfield_spectrum, get_farfield_wave_elevation and get_velocity.
    // 85 deg keeps sec(theta) bounded (~11.5) while reaching well past the
    // diverging-wave decay knee; the Filon-induced |PQ|^2 falloff handles the tail.
    static constexpr double kFarFieldThetaMax = 85.0 * M_PI / 180.0;

    LinearInterpolator deflection;
    double k_multiplier = 1.0;
    double param1 = -0.8; // bow coef
    double param2 = -2;   // fore coef
    double param3 = 1.0;  // aft coef
    double param4 = 3.8;  // wl arc coef
    double param5 = 0;    // wl arc sq coef

    // compact phase model (replaces param4/param5 when use_compact_phase=true)
    bool use_compact_phase = false;
    double phase_alpha = 3.8;
    double phase_gamma = 0.0;

    std::vector<double> deflection_field;

    double catamaran_separation = 0.0;

    static constexpr bool use_filon = true;
    DemihullMethod demihull_method = DemihullMethod::Yeung;
    TransomShapeMethod transom_hollow_method = TransomShapeMethod::DoctorsDynamic;
    TransomShapeMethod transom_dry_method = TransomShapeMethod::DoctorsDynamic;



public:

    ThinShip(std::shared_ptr<Hull> ship_hull)
    : hull(ship_hull)
    {
        env = hull->get_environment();

        grid = std::make_unique<ProfileGrid>(hull);

        deflection = LinearInterpolator(
            {0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.08, 1.28, 1.47, 1.57, 1.67, 1.87, 2.06, 2.26, 2.45, 2.65, 2.85, 3.04, 3.142},
            {0.0, 2.58787224e-02,    -2.07058187e-02,     -4.11867727e-02,    -5.41260291e-02,   -4.65911214e-02,    -2.29935673e-02,    -5.64612831e-03,    -7.46950397e-04,     0.0,            -2.56595981e-02,    -1.69933303e-01,    -4.37234271e-01,    -5.45162121e-01,    -6.96964561e-01,    -7.83925538e-01,    -9.04849779e-01,  -1.11388172,        -1.35}
        );

        for (auto& x : deflection.xi) x -= M_PI_2;

//        std::cout << "Made Basic-Michell calculator for the hull" << std::endl;
    }

    ~ThinShip() {}

    void set_transom_method(TransomShapeMethod hollow, TransomShapeMethod dryness) {
        transom_hollow_method = hollow;
        transom_dry_method = dryness;
        invalidate_ff_spectrum();
    }

    void set_params(double p1, double p2, double p3, double p4, double p5) {
        param1 = p1;
        param2 = p2;
        param3 = p3;
        param4 = p4;
        param5 = p5;
        use_compact_phase = false; // legacy params take over
        invalidate_ff_spectrum();

//        deflection = LinearInterpolator({-M_PI_2, -M_PI_4, 0, M_PI_4, M_PI_2}, {0, p1, 0, p2, p3});
    }

    void set_deflection_table(const std::vector<double>& xs, const std::vector<double>& ys) {
        deflection = LinearInterpolator(xs, ys);
        for (auto& x : deflection.xi) x -= M_PI_2;
        invalidate_ff_spectrum();
    }

    void set_phase_model(double alpha, double gamma) {
        phase_alpha = alpha;
        phase_gamma = gamma;
        use_compact_phase = true;
        invalidate_ff_spectrum();
    }

    void set_deflection_field(const std::vector<double>& field) {
        deflection_field = field;
        invalidate_ff_spectrum();
    }

    void set_grid(std::size_t nx, std::size_t nz) {

        grid->setup(nx, nz);
        invalidate_ff_spectrum();
    }

    ProfileGrid* get_grid() const {
        return grid.get();
    }

    void set_demihull(double separation) {
        catamaran_separation = separation;
    }

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
    }

    void fix_stuff(bool phase, bool bl, bool depth) {

        fix_phase = phase;
        fix_bl = bl;
        fix_depth = depth;
        invalidate_ff_spectrum();   // k_multiplier / grid slopes below feed the far-field spectrum

        if (fix_phase) {

            double arc_ratio = 0; // ignore for immersed bodies
            double waterplane_ratio = 0;

            if (hull->get_length_wl() != 0.0) {
                arc_ratio = std::max(hull->get_arc_shoulder() / std::abs(hull->get_x_shoulder() - hull->get_fp()) - 1.0, 0.0);
                waterplane_ratio = (1 - hull->get_waterplane_coef());
            }

            double phase_base;
            if (use_compact_phase) {
                // compact: k = 1 + alpha * arc_ratio * (1 + gamma * (1-CWP))
                phase_base = phase_alpha * arc_ratio * (1.0 + phase_gamma * waterplane_ratio);
            } else {
                // legacy: k = arc_ratio*p4 + (1-CWP)*p5 + 1
                phase_base = arc_ratio * param4 + waterplane_ratio * param5;
            }

            if (is_demihull()) {
                auto x_in = std::max(0.0, catamaran_separation - hull->get_beam_wl()) / catamaran_separation;
                auto temp = 6.8499340193314706E-01 + 3.4147198623042684E+01 * std::pow(6.5110604379930515E-03, x_in);
                phase_base *= std::max(temp, 1.0);
            }

            k_multiplier = 1.0 + phase_base;
//            std::cout << "Fixed phase by ratio " << k_multiplier << std::endl;
        }

        if (fix_bl) {

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

//            std::cout << "Fixed boundary layer" << std::endl;
        }
    }

    bool is_demihull() const {
        return catamaran_separation > hull->get_beam_wl();
    }

    double get_wave_resistance(double speed) {

        if (speed <= 0 || !grid->is_valid()) {
            return 0;
        }

        auto k = env->get_gravity() / sq(speed);

        if (fix_phase) {
            k *= k_multiplier;
        }

        // Michell's far-field integral in theta-parametrization (lam = sec theta,
        // dlam = sec theta tan theta dtheta) gives
        //   CW = (4/pi) * integral_0^{theta_max} sec^3(theta) * |PQ|^2 * M(theta) dtheta
        // where M(theta) is the catamaran-separation modulation factor (1 for monohulls).
        // The lam=1 transverse-wave singularity vanishes in this form (sec^3(0) = 1);
        // the diverging-wave behaviour at theta -> pi/2 is controlled by the Filon-induced
        // |PQ|^2 falloff, so a uniform Simpson sum is appropriate.
        // 501 points => 500 even subintervals, so composite Simpson 1/3 applies cleanly.
        std::size_t constexpr N = 501;
        tabulate_farfield_spectrum(speed, static_cast<int>(N));
        if (ff_spectrum_.size() != N) {
            return 0;
        }

        double const dtheta = kFarFieldThetaMax / static_cast<double>(N - 1);

        auto modulation = [&](double secT, double secTanT) -> double {
            if (!is_demihull()) {
                return 1.0;
            }
            double const phase = k * catamaran_separation * secTanT;  // k * sep * lam * tan(theta)
            (void)secT;
            switch (demihull_method) {
                case DemihullMethod::Yeung: return 2.0 * (1.0 + std::cos(phase));
                case DemihullMethod::Tuck:  return sq(2.0 * std::cos(0.5 * phase));
                default:                    return 1.0;
            }
        };

        double s_odd = 0.0;
        double s_even = 0.0;
        double f_first = 0.0;
        double f_last = 0.0;

        for (std::size_t i = 0; i < N; i++) {
            auto const& e = ff_spectrum_[i];
            double const PQ_sq = sq(e.P) + sq(e.Q);
            double const fi = e.secT * e.sec2T * PQ_sq * modulation(e.secT, e.secTanT);

            if (i == 0)            f_first = fi;
            else if (i == N - 1)   f_last  = fi;
            else if (i & 1u)       s_odd  += fi;
            else                   s_even += fi;
        }

        double const integral = (dtheta / 3.0) * (f_first + f_last + 4.0 * s_odd + 2.0 * s_even);
        double const Rwp = integral * (4.0 / M_PI) * env->get_density() * sq(speed * k);

        if (std::isnan(Rwp)) {
            return 0;
        }

        return Rwp;
    }

    /// Near-field / far-field wave-drag ratio diagnostic (Batch 3 acceptance gate): the
    /// causal near-field pressure drag integrates to the far-field Havelock result up to
    /// the documented quadrature/truncation differences (t-grid vs theta-grid, d_angle
    /// windows). Expect ~2-3% deviation on a Wigley, ~5-10% on transom hulls.
    /// Returns 0 when either side is zero/invalid.
    double get_nearfar_drag_ratio(double speed) {
        double const near = get_drag_lift_torque(speed).x;
        double const far = get_wave_resistance(speed);
        if (far == 0.0 || std::isnan(far) || std::isnan(near)) {
            return 0.0;
        }
        return near / far;
    }

    glm::dvec3 get_drag_lift_torque(double speed, std::vector<double>* out_pressure_field = nullptr) {

        if (speed <= 0 || !grid->is_valid() || hull->get_fn(speed) < 0.1) {
            return {0, 0, 0};
        }

        auto k = env->get_gravity() / sq(speed);

        if (fix_phase) {
            k *= k_multiplier;
        }

        //
        auto const cg = hull->get_cg(true);
        // Paper Eq. (29c) dimensional prefactor: -(4/pi) * rho * V^2 * k^2.
        // The previous factor of 4 above was spurious and overshot the near-field drag by 4x.
        // Batch 3 sign anchor: the prefactor keeps its sign for the causal kernel — the
        // legacy accumulation equals +causal on symmetric hulls (test_causal_filon.cpp
        // check 3), so -4/pi is what preserves the validated symmetric-hull outputs and
        // the Havelock far-field drag cross-check.
        const double amp = (-4.0/M_PI) * env->get_density() * sq(k * speed);

        // integration specifics
        double constexpr eps = 1e-5;
        double constexpr dt = 5e-2;
        std::size_t N = 5.0 / dt;

        // frames influences
        std::vector<double> Fs(grid->size_x(), 0.0); // frame influence (over the depth)
        std::vector<double> Ss(grid->size_x(), 0.0); // wavemaking influence (fwd from the current frame)
        std::vector<double> As, Bs;                  // causal prefix quadratures

        // output
        double drag = 0;
        double lift = 0;
        double torque = 0;

        if (out_pressure_field) {
            // assign (not resize): a reused same-size buffer would keep stale values in cells that the
            // loop below leaves unwritten. Zero the whole diagnostic field each call.
            out_pressure_field->assign(grid->size_x() * grid->size_z(), 0.0);
        }

        // integration via summing
        for (std::size_t i = 0; i < N; i++) {

            double const t_prev = i*dt;
            double const t_next = (i+1)*dt;
            double const t = (t_prev + t_next) * 0.5;
            double const lam = std::sqrt(1 + sq(t));

            double const d_angle = k * lam * grid->get_dx() * grid->size_x();
            double constexpr d_limit = 1.0 / (2*5.0*M_PI);
            if (d_angle * d_limit > 1.0) {
                continue;
            }

            // get frames influence
            for (std::size_t xi = 0; xi < grid->size_x(); xi++) {
                Fs[xi] = get_frame(xi, k, lam);
            }

            // Causal accumulation (Batch 3): a station is influenced by the sources
            // forward of it only. The transom hollow lies aft of every hull station,
            // so it contributes no hull pressure (hence no lift/torque leak); its
            // drag enters as the far-field-consistent spectral delta below.
            accumulate_causal_influence(k * lam, Fs, As, Bs, Ss);

            if (grid->has_transom()) {
                // Hollow drag from the same closed forms as get_PQ, so the near/far
                // pairing cannot drift. get_PQ integrates stern→bow, so its hull term
                // is (-A_full, -B_full); the per-λ near-field drag of the hull rows is
                // -amp·dt·λ·0.5·(A_full²+B_full²) per side, and the hollow adds the
                // remainder of the far-field drag identity
                //   0.5·[(cosF_hull+dC)² + (sinF_hull+dS)²] = 0.5·|PQ|².
                auto const dCS = get_transom_spectrum_delta(k, lam);
                double const cosF_hull = -As[0];
                double const sinF_hull = -Bs[0];
                double const spec_delta = sq(cosF_hull + dCS.first) + sq(sinF_hull + dCS.second)
                                        - (sq(cosF_hull) + sq(sinF_hull));
                drag += -amp * dt * lam * 0.5 * spec_delta;
            }

            double const lift_old = lift;

            for (std::size_t xi = 0; xi < grid->size_x(); xi++) {

                for (std::size_t zi = 0; zi < grid->size_z(); zi++) {

                    double const depth = grid->get_z(zi) - grid->get_z(0); // negative
                    double const dp = amp * dt * lam * Ss[xi] * std::exp(depth * k * sq(lam));
                    auto const dA = grid->get_area(xi, zi);
                    // Sign matches the Havelock far-field result; cross-checked vs get_wave_resistance.
                    drag += dp * dA.x;
                    lift += dp * dA.z;
                    // bow-up pitch moment = Δx·Fz − Δz·Fx (the Δz·Fx term was + — wrong sense)
                    torque += dp * (dA.z * (grid->get_x(xi) - cg.x) - dA.x * (grid->get_z(zi) - cg.z));
                    if (out_pressure_field) {
                        (*out_pressure_field)[grid->index(xi, zi)] += dp;
                    }
                }
            }

//            if (t > 3.0 && std::abs(lift - lift_old) < std::abs(lift * eps)) {
//                break;
//            }

        }

        // account for port + starboard sides (paper Eq. 29c is per-side, matching LocalFlow's convention)
        drag *= 2.0;
        lift *= 2.0;
        torque *= 2.0;

        if (std::isnan(drag) || std::isnan(lift) || std::isnan(torque)) {
            drag = 0;
            lift = 0;
            torque = 0;
        }

        // add the friction force contribution to the torque, approximate
        torque += hull->get_force(speed, get_friction_coef(speed)) * (hull->get_wetted_centroid().z - cg.z);

        return {drag, lift, torque};
    }


    /// Wave-component centerline elevation E_W(x, 0, 0).
    /// Noblesse 2009 Eqs. 5-6: at z=0 the hull-pressure integral equals the free-surface
    /// elevation, divided by rho*g to convert pressure (Pa) to height (m). The integrand
    /// matches get_drag_lift_torque with exp(depth*k*lam^2) collapsing to 1 at z=0.
    /// Result is in metres. Reuses the per-lam Filon spectrum machinery; Ss is linearly
    /// interpolated in xi-space for arbitrary x_eval.
    std::vector<double> get_wave_elevation_centerline(double speed, const std::vector<double>& xs_eval) {

        std::vector<double> result(xs_eval.size(), 0.0);

        if (speed <= 0 || !grid->is_valid() || hull->get_fn(speed) < 0.1) {
            return result;
        }

        auto k = env->get_gravity() / sq(speed);
        if (fix_phase) {
            k *= k_multiplier;
        }

        const double amp = (-4.0/M_PI) * env->get_density() * sq(k * speed);

        double constexpr dt = 5e-2;
        std::size_t const N = static_cast<std::size_t>(5.0 / dt);

        std::size_t const Nx = grid->size_x();
        std::vector<double> Fs(Nx, 0.0);
        std::vector<double> Ss(Nx, 0.0);
        std::vector<double> As, Bs;   // causal prefix quadratures

        auto const fp = hull->get_fp();
        auto const fwd = hull->get_fwd();
        double const dx_grid = grid->get_dx();
        double const x0 = grid->get_x(0);

        for (std::size_t i = 0; i < N; i++) {

            double const t_prev = i*dt;
            double const t_next = (i+1)*dt;
            double const t = (t_prev + t_next) * 0.5;
            double const lam = std::sqrt(1 + sq(t));

            double const d_angle = k * lam * grid->get_dx() * Nx;
            double constexpr d_limit = 1.0 / (2*5.0*M_PI);
            if (d_angle * d_limit > 1.0) {
                continue;
            }

            for (std::size_t xi = 0; xi < Nx; xi++) {
                Fs[xi] = get_frame(xi, k, lam);
            }

            // Hull-only causal envelopes (the transom hollow enters the drag identity
            // only; the near-field elevation stays a hull-pressure quantity).
            accumulate_causal_influence(k * lam, Fs, As, Bs, Ss);

            for (std::size_t e = 0; e < xs_eval.size(); e++) {

                double const x_e = xs_eval[e];
                double const xi_real = (x_e - x0) / dx_grid;

                // Interpolate the SMOOTH A/B envelopes and apply the exact phase at the
                // evaluation point — interpolating the oscillatory Ss itself aliases.
                // Downstream of the last station the full-hull envelope propagates with
                // the true phase.
                double const s_e = (fp - x_e) * fwd;
                double A_e = 0.0, B_e = 0.0;

                if (xi_real <= 0.0) {
                    // fwd>0 ⇒ downstream of stern: full-hull envelope. fwd<0 ⇒ forward of bow: 0.
                    if (fwd > 0) { A_e = As[0]; B_e = Bs[0]; }
                } else if (xi_real >= static_cast<double>(Nx - 1)) {
                    if (fwd < 0) { A_e = As[Nx - 1]; B_e = Bs[Nx - 1]; }
                } else {
                    std::size_t const xi0 = static_cast<std::size_t>(xi_real);
                    double const frac = xi_real - static_cast<double>(xi0);
                    A_e = (1.0 - frac) * As[xi0] + frac * As[xi0 + 1];
                    B_e = (1.0 - frac) * Bs[xi0] + frac * Bs[xi0 + 1];
                }

                double const Ss_at_xe = std::cos(k * lam * s_e) * A_e + std::sin(k * lam * s_e) * B_e;

                result[e] += amp * dt * lam * Ss_at_xe;  // single-side pressure-like contribution (Pa)
            }
        }

        // Both hull sides contribute equally at y=0; then pressure → elevation height.
        double const rho_g = env->get_density() * env->get_gravity();
        double const scale = 2.0 / rho_g;
        for (auto& r : result) {
            if (std::isnan(r)) r = 0.0;
            r *= scale;
        }

        return result;
    }


    /// Tabulate the far-field hull spectrum (P, Q) on a uniform theta grid in
    /// [0, kFarFieldThetaMax] (85°). Cached per speed; shared by get_wave_resistance,
    /// get_farfield_wave_elevation and get_velocity so they integrate over the same
    /// grid. Re-tabulates if the requested speed or sample count differs from the cache.
    void tabulate_farfield_spectrum(double speed, int n_theta = 401) {

        if (speed <= 0 || !grid->is_valid()) {
            ff_spectrum_.clear();
            ff_spectrum_speed_ = -1.0;
            return;
        }

        auto k = env->get_gravity() / sq(speed);
        if (fix_phase) {
            k *= k_multiplier;
        }

        // re-use cache when speed (and therefore k) is unchanged
        if (!ff_spectrum_.empty() && std::abs(speed - ff_spectrum_speed_) < 1e-9
            && static_cast<int>(ff_spectrum_.size()) == n_theta) {
            return;
        }

        ff_spectrum_.clear();
        ff_spectrum_.reserve(n_theta);

        double const dtheta = kFarFieldThetaMax / static_cast<double>(n_theta - 1);

        for (int i = 0; i < n_theta; i++) {
            double const theta = dtheta * i;
            double const cT = std::cos(theta);
            double const sT = std::sin(theta);
            double const secT = 1.0 / cT;
            double const sec2T = sq(secT);
            double const sec4T = sq(sec2T);
            // sec^2(theta) * sin(theta) = sec(theta) * tan(theta); used in y-phase.
            double const secTanT = secT * (sT / cT);
            auto pq = get_PQ(k, secT);
            ff_spectrum_.push_back({theta, secT, sec2T, sec4T, secTanT, pq.first, pq.second});
        }

        ff_spectrum_speed_ = speed;
        ff_spectrum_k_ = k;
    }

    /// Far-field wave elevation at world point (x_field, y_field) on the free surface,
    /// using the Tuck-Scullen-Lazauskas Fourier-Kochin form (their Eq. 5) with the
    /// hull spectrum cached by tabulate_farfield_spectrum. Eddy-viscosity damping
    /// (TSL Eq. 7) suppresses very short diverging waves; pass viscosity=0 to disable.
    ///
    /// LIMITATIONS:
    ///   - far-field only: no near-field correction Z_L (Newman's G_L approximation
    ///     is out of scope), so the field is inaccurate within roughly 1 ship length
    ///     of the hull.
    ///   - returns 0 forward of the bow (no causally connected wave sources).
    double get_farfield_wave_elevation(double speed, double x_field, double y_field,
                                       double viscosity = 5e-3) {

        if (speed <= 0 || ff_spectrum_.empty()) {
            return 0.0;
        }
        // Speed-consistency guard: the cached spectrum (and k0 below) is keyed on the tabulated speed.
        // If asked for a different speed, refresh it rather than silently using the wrong-k spectrum.
        if (std::abs(speed - ff_spectrum_speed_) >= 1e-9) {
            tabulate_farfield_spectrum(speed);
            if (ff_spectrum_.empty()) { return 0.0; }
        }

        double const k0 = ff_spectrum_k_;
        double const fp = hull->get_fp();
        double const L = hull->get_length_wl();
        int const fwd = hull->get_fwd();

        // Downstream distance from FP. >0 aft of bow, ∈[0,L] on hull, >L aft of stern.
        double const s_field = (fp - x_field) * fwd;

        if (s_field < 0.0) {
            return 0.0;  // forward of bow: no Kelvin wake reaches here
        }

        // Damping distance: how far the waves have propagated from the stern to the observer.
        // For observers on the hull (s_field < L), use 0; downstream of stern, use s_field - L.
        double const s_eff = std::max(0.0, s_field - L);
        double const damping_pre = -2.0 * viscosity * sq(k0) * s_eff / speed;

        double const dtheta = kFarFieldThetaMax / static_cast<double>(ff_spectrum_.size() - 1);

        double sum = 0.0;
        for (std::size_t i = 0; i < ff_spectrum_.size(); i++) {
            auto const& e = ff_spectrum_[i];
            double const alpha = k0 * e.secT * s_field;
            double const beta = k0 * e.secTanT * y_field;
            double const damping = (viscosity > 0.0) ? std::exp(damping_pre * e.sec4T) : 1.0;
            double const integrand = e.sec4T * std::cos(beta)
                                     * (e.P * std::sin(alpha) - e.Q * std::cos(alpha))
                                     * damping;
            // trapezoidal weights: 0.5 at endpoints, 1 elsewhere
            double const w = (i == 0 || i == ff_spectrum_.size() - 1) ? 0.5 : 1.0;
            sum += w * integrand;
        }
        sum *= dtheta;

        return -(4.0 * sq(k0) / M_PI) * sum;
    }

    /// Bulk evaluation of the far-field wave field on the cartesian outer product xs × ys.
    /// Returns a flattened (ny × nx) array (x varies fastest). Cache the spectrum once,
    /// then parallelise over field points with OpenMP.
    std::vector<double> get_farfield_wave_field(double speed,
                                                const std::vector<double>& xs,
                                                const std::vector<double>& ys,
                                                double viscosity = 5e-3,
                                                int n_theta = 401) {

        std::size_t const N = xs.size() * ys.size();
        std::vector<double> result(N, 0.0);

        if (speed <= 0 || !grid->is_valid() || xs.empty() || ys.empty()) {
            return result;
        }

        tabulate_farfield_spectrum(speed, n_theta);

        #pragma omp parallel for collapse(2)
        for (std::size_t j = 0; j < ys.size(); j++) {
            for (std::size_t i = 0; i < xs.size(); i++) {
                result[j * xs.size() + i] = get_farfield_wave_elevation(speed, xs[i], ys[j], viscosity);
            }
        }

        return result;
    }


    double get_transom_resistance(double speed) {

        if (!grid->has_transom()) {
            return 0;
        }

        // calculation based on dryness coef
        double const beam_0 = 2 * grid->get_beam_transom(0);
        double const draft_0 = grid->get_draft_transom_approx();
        double const dry = get_transom_dryness(speed);
        double const draft_dynamic = draft_0 * (1 - dry);
        double R_transom_simple = 0.9 * env->get_density() * env->get_gravity() * beam_0 * (0.5*(sq(draft_0) - sq(draft_dynamic)));
        if (is_demihull()) {
            R_transom_simple *= 2;
        }

        // calculation based on each beam slice
        double const transom_wl = get_transom_wl(speed);
        std::vector<double> dfs;

        for (std::size_t zi = 0; zi < grid->size_z() - 1; zi++) {
            double const y = grid->get_beam_transom(zi);
            if (y >= 1e-4 && grid->get_z(zi) > transom_wl) {

                // loss of hydrostatic pressure, p = rho * g * depth
                double const pressure = env->get_density() * env->get_gravity() * std::abs(grid->get_z(zi+1) - grid->get_z(0));
                // dF/dz = p * B(Z)
                double const df = pressure * (2*y);

                dfs.push_back(df);
            } else {
                break;
            }
        }

        double R_transom = std::abs(trapz(grid->get_dz(), dfs));

        if (is_demihull()) {
            R_transom *= 2;
        }

        return (R_transom + R_transom_simple) * 0.5;
    }


    /*
    double get_viscous_pressure_resistance_alt(double speed) {

        auto F = hull->get_fn(speed);
        auto fp = hull->get_fp();
        auto linv = 1.0 / hull->get_length_wl();
        auto wl = zs[0];
        auto coef = 1.0 /(4*M_PI) * (env->density * sq(speed));
        auto hp = std::abs(xs[1]-xs[0]) * 1e-2;
        auto hp3 = cub(hp);
        auto hp6 = sq(cub(hp3));

        // gradient of the potential gives velocity
        double R = 0.0;
        double dS = std::abs((xs[1]-xs[0])*(zs[1]-zs[0]));
        double dSn = dS * sq(linv);
        double sum_area = 0;
        for (std::size_t xi1 = 0; xi1 < xs.size(); xi1++) {
            for (std::size_t zi1 = 0; zi1 < zs.size(); zi1++) {

                if (get_beam(xi1,zi1) < 1e-4) {
                    continue;
                }

                auto index1 = zi1 * xs.size() + xi1;
                auto pt1 = glm::dvec3(fp - xs[xi1], 0, zs[xi1] - wl) * linv;

                for (std::size_t xi2 = 0; xi2 < xs.size(); xi2++) {
                    for (std::size_t zi2 = 0; zi2 < zs.size(); zi2++) {

                        if (get_beam(xi2,zi2) < 1e-4) {
                            continue;
                        }

                        auto index2 = zi2 * xs.size() + xi2;
                        if (index2 <= index1) {
                            continue;
                        }

                        auto pt2 = glm::dvec3(fp - xs[xi2], 0, zs[xi2] - wl) * linv;
                        auto r = std::sqrt(sq(pt1.x - pt2.x) + sq(pt1.y - pt2.y) + sq(pt1.z - pt2.z));
//                        if (r <= 1e-5 || r > 0.4) {
//                            continue;
//                        }
                        auto r1 = std::sqrt(sq(pt1.x - pt2.x) + sq(pt1.y - pt2.y) + sq(pt1.z + pt2.z));
                        auto eps = hp6 / (hp3 + cub(r));
                        auto eps1 = hp6 / (hp3 + cub(r1));
                        auto eps_ = hp6 / (hp3*hp + sq(sq(r1)));
                        auto psi = -(pt1.z + pt2.z) / (r1 + std::abs(pt1.x - pt2.x) + 1e-8);
                        auto sigma1 = sign(pt1.x - pt2.x) * get_slope(xi1,zi1).x;
                        auto sigma2 = sign(pt2.x - pt1.x) * get_slope(xi2,zi2).x;
                        auto M = 1/(cub(r)+eps) - 1/(cub(r1)+eps1) + 2/cub(sq(F)+r1) * (1 + sq(F)/std::sqrt(sq(r1)+eps_)*(1 + 2*psi - 2.3*sq(psi) * (sq(F)-2*r1)/(sq(F)+r1)));
                        auto MM = coef * dSn * (std::abs(pt1.x - pt2.x)*M + 2*psi*sq(F)/cub(sq(F)+r1) * (1+4.6*psi+sq(F)/(r1+eps1)));
//                        pressures[index1] += sigma1 * MM;
//                        pressures[index2] += sigma2 * MM;

//                        auto beta_from_slope1 = M_PI_2 - std::atan(get_slope(xi1,zi1).x);
                        R += (sigma1 * MM) * get_normal(xi1,zi1).x * dS; //std::cos(beta_from_slope1); //get_nx(xi1,zi1);

//                        auto beta_from_slope2 = M_PI_2 - std::atan(get_slope(xi2,zi2).x);
                        R += (sigma2 * MM) * get_normal(xi2,zi2).x * dS; //std::cos(beta_from_slope2); //get_nx(xi2,zi2);

                        sum_area += 2*dS;
                    }
                }
            }
        }

        R = std::abs(R/sum_area*hull->get_wetted_area());
        return R;
    }
*/
    double get_friction_coef(double speed) {

        // initial waterline contribution
        double const Rn0 = hull->get_rn(speed);
        double const Cf0 = 0.075 / sq(std::log10(Rn0) - 2);
        int sum_form = int(4*grid->size_x());
        double sum_coef = Cf0 * sum_form;

        // if we have large variations of form accros draft, this calculates friction for each depth
        // it should be for streamlines with included pressure effect, but nobody cares
        for (std::size_t zi = 1; zi < grid->size_z(); zi++) {

            int form_start = -1;

            for (std::size_t xi = 0; xi < grid->size_x(); xi++) {

                double const y = grid->get_beam(xi, zi);

                if (y > 1e-4 && form_start < 0) {

                    form_start = int(xi);

                } else if (y <= 1e-4 && form_start >= 0) {

                    double const len = std::abs(grid->get_x(xi) - grid->get_x(form_start));
                    auto Rn = hull->get_rn(speed, len);
                    auto Cf = 0.075 / sq(std::log10(Rn) - 2);
                    sum_coef += Cf * (xi-form_start);
                    sum_form += int(xi) - form_start;
                    form_start = -1;

                }
            }
        }

        return (sum_coef / sum_form);
    }

private:

    // Doctors, L. J. (2007). A Numerical Study of the Resistance of Transom-Stern Monohulls. Ship Technology Research, 54(3), 134–144. doi:10.1179/str.2007.54.3.005
    double get_transom_dryness(double speed) {

        double const transom_draft = grid->get_draft_transom_approx();

        if (speed <= 0 || !grid->has_transom() || transom_draft <= 0.0) {
            return 0.0;
        }

        // Doctors static:
        double constexpr C1s = 0.002472;
        double constexpr C2s = 1.862;
        double constexpr C3s = 0.2859;
        double constexpr C4s = 0.3588;

        // Doctors dynamic:
        auto constexpr C1d = 0.004856;
        auto constexpr C2d = 1.821;
        auto constexpr C3d = 0.1990;
        auto constexpr C4d = 0.3126;

        double const Fn_T = speed / std::sqrt(env->get_gravity() * transom_draft);
        double const Fn_B = speed / std::sqrt(env->get_gravity() * hull->get_beam_transom());
        double const BT = hull->get_beam_transom() / transom_draft;
        double const Rn = std::sqrt(env->get_gravity() * cub(transom_draft)) / env->get_viscosity();

        double const eta_doctors_static = clamp(C1s * std::pow(Fn_T, C2s) * std::pow(BT, C3s) * std::pow(Rn, C4s), 0.0, 1.0);
        double const eta_doctors_dynamic = clamp(C1d * std::pow(Fn_T, C2d) * std::pow(BT, C3d) * std::pow(Rn, C4d), 0.0, 1.0);

//        double const eta_maki = clamp(0.1578 * std::pow(Fn_T, 1.966), 0.0, 1.0);
//        return eta_maki; // TEMP SOLUTION

         // Robards:
        double const a0 = -0.0002444*cub(BT) + 0.003303*sq(BT) - 0.01494*(BT) + 0.02260;
        double const a1 = -3.381*(BT) - 2.609;
        double const a2 = 0.7594*sq(sq(BT)) - 5.981*cub(BT) + 19.10*sq(BT) - 29.35*(BT) + 45.70;
        double const a3 = -0.002399*cub(BT) + 0.03279*sq(BT) - 0.1486*(BT) + 0.2251;
        double const a4 = -0.1111*sq(BT) + 0.9967*(BT) + 0.07370;
        double eta_robards = (0.5 * (1 + std::tanh(3.765 + BT*(0.145*Fn_T - 0.377) - 1.77*Fn_T))
                                         + BT * (a0 * std::sin(a1*Fn_B + a2*std::sqrt(Fn_T)) / (a3 + sq(sq(a4 - Fn_T))))
                                         - 0.008*(Rn/50000 - Fn_T) / (0.1*std::pow(BT,0.333) + sq(sq(1.5*std::pow(BT,0.333) - Fn_T)))
                                          );
        eta_robards = clamp(1 - eta_robards, 0.0, 1.0);

        switch (transom_dry_method) {
        case TransomShapeMethod::DoctorsStatic: return eta_doctors_static;
        case TransomShapeMethod::DoctorsDynamic: return eta_doctors_dynamic;
        case TransomShapeMethod::Robards: return eta_robards;
        default:  return (eta_doctors_dynamic + eta_robards) * 0.5;
        }
    }

    double get_transom_wl(double speed) {

        double const transom_draft = grid->get_draft_transom_approx();
        double const eta_dry = get_transom_dryness(speed);
        auto z_transom = grid->get_z(0) - eta_dry * transom_draft;
        return z_transom;
    }

    // Doctors, L. J. (2007). A Numerical Study of the Resistance of Transom-Stern Monohulls. Ship Technology Research, 54(3), 134–144. doi:10.1179/str.2007.54.3.005
    double get_transom_hollow_length(double speed) {

        double const transom_draft = grid->get_draft_transom_approx();

        if (speed <= 0 || !grid->has_transom() || transom_draft <= 0.0) {
            return 0;
        }

        // Doctors static:
        auto constexpr C1s = 0.6095;
        auto constexpr C2s = 2.733;
        auto constexpr C3s = 0.3468;
        auto constexpr C4s = -0.1514;

        // Doctors dynamic:
        auto constexpr C1d = 0.2491;
        auto constexpr C2d = 3.107;
        auto constexpr C3d = 0.1598;
        auto constexpr C4d = -0.1225;

        auto Fn = speed / std::sqrt(env->get_gravity() * transom_draft);
        auto BT = hull->get_beam_transom() / transom_draft;
        auto Rn = std::sqrt(env->get_gravity() * cub(transom_draft)) / env->get_viscosity();
        auto L_doctors_static = transom_draft * C1s * std::pow(Fn, C2s) * std::pow(BT, C3s) * std::pow(Rn, C4s);
        auto L_doctors_dynamic = transom_draft * C1d * std::pow(Fn, C2d) * std::pow(BT, C3d) * std::pow(Rn, C4d);

        auto b0 = -0.4574*cub(BT) + 3.377*sq(BT) - 8.585*(BT) + 11.51;
        auto b1 = -0.5960*sq(sq(BT)) + 6.174*cub(BT) - 21.77*sq(BT) + 30.17*(BT) + 9.929;
        auto b2 = -0.03840*cub(BT) + 0.3986*sq(BT) - 1.420*(BT) + 2.529;
        auto L_robards = transom_draft * 0.0113 * std::exp(1.9*Fn - 1.1223) + 1.2 + BT*0.03*std::sin(b0*Fn + b1)/(0.08 +sq(sq(b2-Fn)));

        switch (transom_hollow_method) {
        case TransomShapeMethod::DoctorsStatic: return L_doctors_static;
        case TransomShapeMethod::DoctorsDynamic: return L_doctors_dynamic;
        case TransomShapeMethod::Robards: return L_robards;
        default:
            return L_doctors_dynamic * 0.7 + L_robards * 0.3;
            //return (L_doctors_dynamic * Fn + L_robards * std::max(0.0, 6 - Fn)) / (Fn + std::max(0.0, 6 - Fn));
        }
    }

    // Causal near-field influence of the hull sources for one wave component kl = k*lam.
    // s = (fp - x)*fwd is the downstream distance from the FP; the station with the
    // smallest s (grid index size_x()-1 under the +X-forward convention asserted by
    // Resistance) is the bow. Waves trail aft, so a station is influenced by the sources
    // forward of it:
    //   Ss[xi] = ∫_{s_bow}^{s_i} F(s') cos(kl (s_i - s')) ds'
    //          = cos(kl s_i) A[xi] + sin(kl s_i) B[xi],
    //   A[xi]  = ∫_{s_bow}^{s_i} F cos(kl s') ds',  B[xi] = ∫_{s_bow}^{s_i} F sin(kl s') ds'.
    // A and B are smooth (non-oscillatory in s_i) prefix quadratures, accumulated in one
    // bow→stern march from the exact linear-interpolant segment integrals (tools.h
    // filon_cos_seg/filon_sin_seg). O(N), no odd-count constraint, no parity fill.
    // Validated against a high-res reference in src/tests/test_causal_filon.cpp.
    void accumulate_causal_influence(double kl, const std::vector<double>& Fs,
                                     std::vector<double>& A, std::vector<double>& B,
                                     std::vector<double>& Ss) const {

        std::size_t const Nx = grid->size_x();
        A.assign(Nx, 0.0);
        B.assign(Nx, 0.0);
        Ss.assign(Nx, 0.0);

        if (Nx < 2) {
            return;
        }

        double const fp = hull->get_fp();
        int const fwd = hull->get_fwd();

        for (std::size_t xi = Nx - 1; xi-- > 0; ) {
            double const s1 = (fp - grid->get_x(xi + 1)) * fwd;  // forward neighbor (smaller s)
            double const s0 = (fp - grid->get_x(xi)) * fwd;
            A[xi] = A[xi + 1] + filon_cos_seg(Fs[xi + 1], Fs[xi], s1, s0, kl);
            B[xi] = B[xi + 1] + filon_sin_seg(Fs[xi + 1], Fs[xi], s1, s0, kl);
            Ss[xi] = std::cos(kl * s0) * A[xi] + std::sin(kl * s0) * B[xi];
        }
    }

    // Closed-form spectrum contribution (dC, dS) of the transom hollow: the Doctors/Robards
    // hollow modeled as a constant frame strength Ft over the hollow length aft of the
    // transom, integrated against cos/sin(kl·s). Shared by get_PQ (far field) and the
    // near-field hollow drag delta in get_drag_lift_torque so the two cannot drift.
    std::pair<double, double> get_transom_spectrum_delta(double k, double lam) {

        if (!grid->has_transom()) {
            return {0.0, 0.0};
        }

        double const Ft = get_frame_transom(k, lam);
        double const fp = hull->get_fp();
        double const x_transom = grid->get_x(grid->get_x_transom(0));
        double const speed = std::sqrt(env->get_gravity() / k);
        double const hollow = get_transom_hollow_length(speed);
        double const kl = k * lam;
        double const s_t = (fp - x_transom) * hull->get_fwd();

        double const dC = ( std::sin(kl * (s_t + hollow)) - std::sin(kl * s_t)) * Ft / kl;
        double const dS = (-std::cos(kl * (s_t + hollow)) + std::cos(kl * s_t)) * Ft / kl;

        return {dC, dS};
    }

    double get_frame(std::size_t xi, double k, double lam) {

        std::vector<double> vals(grid->size_z());
        for (std::size_t zi = 0; zi < vals.size(); zi++) {
            double const slope = grid->get_slope(xi, zi);
            double const depth = grid->get_z(zi) - grid->get_z(0); // negative
            vals[zi] = slope * std::exp(depth * k * sq(lam));
        }
        return trapz(grid->get_dz(), vals);
    }

    double get_frame_transom(double k, double lam) {

        if (k <= 0.0) {
            return 0;
        }

        double const speed = std::sqrt(env->get_gravity() / k);
        double const transom_draft = grid->get_draft_transom_approx();
//        double const transom_wl = get_transom_wl(speed);
        double  hollow = get_transom_hollow_length(speed);
        //hollow = std::max(hollow, 0.05*hull->get_length_wl());

        std::vector<double> vals(grid->size_z(), 0.0);
        for (std::size_t zi = 0; zi < vals.size(); zi++) {

            int xi = grid->get_x_transom(zi);
            if (xi < 0 || static_cast<std::size_t>(xi) >= grid->size_x()) {
                break;
            }

            double slope = -1.0 * std::max(transom_draft, grid->get_beam_transom(zi)) / hollow;
            //double slope_geo = grid->get_slope(xi, zi);
            //slope = (slope + slope_geo) * 0.5;
            double const depth = grid->get_z(0) - grid->get_z(zi);
            vals[zi] = std::max(-0.33, slope) * std::exp(-depth * k * sq(lam));
        }

        while (vals.size() < 3) {
            vals.push_back(vals[0]);
        }

        return trapz(grid->get_dz(), vals);
    }

    std::pair<double, double> get_PQ(double k, double lam) {

        double const d_angle = k * lam * grid->get_dx() * grid->size_x();
        double constexpr d_limit = 1.0 / (2*5.0*M_PI);
        double const d_ratio = d_angle * d_limit;
        if (d_ratio > 1.2) {
            return {0.0, 0.0};
        }

        double const amp = d_ratio <= 1.0 ? 1.0 : std::exp(-10.0 * (d_ratio - 1.0)); //lerp(1.0, 0.0, (d_ratio - 1.0)/0.4);
        double const fp = hull->get_fp();
        double const len_start = (fp - grid->get_x(0)) * hull->get_fwd();
        double const len_end = (fp - grid->get_x(grid->size_x() - 1)) * hull->get_fwd();
        std::vector<double> Fs(grid->size_x());

        for (std::size_t xi = 0; xi < Fs.size(); xi++) {
            Fs[xi] = get_frame(xi, k, lam);
        }

        double cosF = filon_cos(Fs, len_start, len_end, k*lam);
        double sinF = filon_sin(Fs, len_start, len_end, k*lam);

        if (grid->has_transom()) {
            // Doctors/Robards hollow closed forms, shared with the near-field drag delta
            // (get_transom_spectrum_delta) — bit-identical refactor of the historical
            // inline cosF_transom_2/sinF_transom_2 expressions.
            auto const dCS = get_transom_spectrum_delta(k, lam);
            cosF += dCS.first;
            sinF += dCS.second;
        }

        return {cosF*amp, sinF*amp};

    }

    double get_PQ_squared(double k, double lam) {
        auto PQ = get_PQ(k, lam);
        return sq(PQ.first) + sq(PQ.second);
    }

public:

    /// Depth of the free-surface depression aft of a ventilated transom, below the calm
    /// waterline, at the world point (X, Y). A dry transom leaves the water resting on the
    /// wetted part of the transom edge, so the surface starts eta_dry * T_transom below the
    /// calm level and recovers over Doctors' hollow length. Zero forward of the transom,
    /// beyond the hollow, and outside the transom beam.
    ///
    /// len_scale stretches (>1) or shortens (<1) the recovery length, for A/B calibration.
    double get_transom_hollow_drop(double speed, double X, double Y, double len_scale = 1.0) {

        if (speed <= 0 || !grid->has_transom()) {
            return 0.0;
        }

        double const T = grid->get_draft_transom_approx();
        double const eta = get_transom_dryness(speed);
        double const L = get_transom_hollow_length(speed) * len_scale;

        if (T <= 0.0 || eta <= 0.0 || L <= 1e-9) {
            return 0.0;
        }

        std::size_t const xi_t = grid->get_x_transom(0);

        if (xi_t >= grid->size_x()) {
            return 0.0;
        }

        // Downstream distance aft of the transom plane.
        double const s = (grid->get_x(xi_t) - X) * hull->get_fwd();

        if (s < 0.0 || s > L) {
            return 0.0;
        }

        // Half-cosine recovery from the transom edge back to the calm level.
        double const f = 0.5 * (1.0 + std::cos(M_PI * s / L));

        // Transverse extent: full inside the transom half-beam, cosine-tapered over its outer
        // fifth so a disc wider than the transom does not see a step.
        double const b = grid->get_beam_transom(0);
        double const a = std::abs(Y);
        double g_y = 1.0;

        if (b > 1e-9) {
            double const b0 = 0.8 * b;
            if (a >= b) {
                g_y = 0.0;
            } else if (a > b0) {
                g_y = 0.5 * (1.0 + std::cos(M_PI * (a - b0) / (b - b0)));
            }
        }

        return eta * T * f * g_y;
    }

    /// Diagnostics for the transom head-release term, in grid/world coordinates:
    /// {transom draft, dryness eta, hollow length, transom station x, calm waterline z}.
    std::vector<double> get_transom_hollow_params(double speed) {

        if (!grid->has_transom()) {
            return {0.0, 0.0, 0.0, 0.0, grid->get_z(0)};
        }

        std::size_t const xi_t = grid->get_x_transom(0);
        double const x_t = (xi_t < grid->size_x()) ? grid->get_x(xi_t) : 0.0;

        return {grid->get_draft_transom_approx(), get_transom_dryness(speed),
                get_transom_hollow_length(speed), x_t, grid->get_z(0)};
    }

    /// Axial speed excess, as a fraction of the ship speed, from the hydrostatic head released
    /// by the transom hollow. A particle under the depression has lost rho*g*dz of head relative
    /// to the same depth far upstream, so Bernoulli gives u = sqrt(speed^2 + 2*g*dz). The excess
    /// is independent of depth — the long-wave limit, in which the depression is wide compared
    /// with the depth of interest — which is why it is nearly uniform over a propeller disc.
    /// Returns (u - speed)/speed, and zero above the depressed surface (that point is in air).
    ///
    /// The same head loss drives get_transom_resistance; this is its velocity-field counterpart,
    /// which the thin-ship kernel and the far-field wave spectrum do not carry.
    double get_transom_head_release(double speed, double X, double Y, double Z,
                                    double len_scale = 1.0) {

        double const dz = get_transom_hollow_drop(speed, X, Y, len_scale);

        if (dz <= 0.0 || Z > grid->get_z(0) - dz) {
            return 0.0;
        }

        return std::sqrt(1.0 + 2.0 * env->get_gravity() * dz / sq(speed)) - 1.0;
    }

    /// Disturbance velocity (u_x, u_y, u_z) at an arbitrary world point (X, Y, Z),
    /// in the ship-fixed frame. The total fluid velocity at the point is
    /// (-speed * fwd, 0, 0) + this vector; for a propeller plane the "axial inflow"
    /// is therefore speed - u_x * fwd (always positive when flow enters from ahead).
    ///
    /// Combines per-cell Rankine source + free-surface mirror image (analytic
    /// gradient, accurate in the near field) with the Havelock wave term from the
    /// cached far-field spectrum (Tuck-Scullen-Lazauskas form). The wave term is
    /// causal: zero forward of the bow.
    glm::dvec3 get_velocity(double speed, double X, double Y, double Z,
                            double viscosity = 5e-3, int n_theta = 121) {

        if (speed <= 0 || !grid->is_valid()) {
            return {0.0, 0.0, 0.0};
        }

        double const k0 = env->get_gravity() / sq(speed);
        int const fwd = hull->get_fwd();
        double const wl = grid->get_z(0);
        double const fp = hull->get_fp();
        double const L = hull->get_length_wl();
        double const dx = std::abs(grid->get_dx());
        double const dz = std::abs(grid->get_dz());
        double const dS = dx * dz;

        double const eps2 = 1e-6 * sq(L);

        // ===== 1. Per-cell Rankine + free-surface image (analytic gradient) =====
        // Centerplane source m(x,z) = -(V/pi) * fwd * (dY/dx); dimensional 1/s.
        // Disturbance potential phi = - sum m*dS/(4 pi) * (1/r - 1/r').
        // Velocity = grad phi = sum m*dS/(4 pi) * (r_vec/r^3 - r_mirror_vec/r'^3).

        double u_x_db = 0.0, u_y_db = 0.0, u_z_db = 0.0;

        #pragma omp parallel for reduction(+:u_x_db,u_y_db,u_z_db) collapse(2)
        for (std::size_t xi = 0; xi < grid->size_x(); xi++) {
            for (std::size_t zi = 0; zi < grid->size_z(); zi++) {

                double const beam = grid->get_beam(xi, zi);
                if (beam < 1e-6) continue;

                double const slope = grid->get_slope(xi, zi);
                double const x_src = grid->get_x(xi);
                double const z_src = grid->get_z(zi) - wl;  // <= 0

                double const m = -fwd * slope * (speed / M_PI);

                double const dxf = X - x_src;
                double const dyf = Y;
                double const dzf_real = (Z - wl) - z_src;
                double const dzf_mir  = (Z - wl) + z_src;

                double const r2    = sq(dxf) + sq(dyf) + sq(dzf_real);
                double const r1_2  = sq(dxf) + sq(dyf) + sq(dzf_mir);
                double const r3    = std::pow(r2   + eps2, 1.5);
                double const r1_3  = std::pow(r1_2 + eps2, 1.5);

                double const coef = m * dS / (4.0 * M_PI);

                u_x_db += coef * (dxf      / r3 - dxf      / r1_3);
                u_y_db += coef * (dyf      / r3 - dyf      / r1_3);
                u_z_db += coef * (dzf_real / r3 - dzf_mir  / r1_3);
            }
        }

        // ===== 2. Havelock wave term from cached far-field spectrum =====
        tabulate_farfield_spectrum(speed, n_theta);

        double u_x_w = 0.0, u_y_w = 0.0, u_z_w = 0.0;

        if (!ff_spectrum_.empty()) {

            double const s_field = (fp - X) * fwd;

            if (s_field >= 0.0) {

                double const Z_dep = std::min(Z - wl, -1e-4);  // depth, <= 0
                double const s_eff = std::max(0.0, s_field - L);
                double const damping_pre = -2.0 * viscosity * sq(k0) * s_eff / speed;
                double const dtheta = kFarFieldThetaMax / static_cast<double>(ff_spectrum_.size() - 1);

                double Sx = 0.0, Sy = 0.0, Sz = 0.0;

                for (std::size_t i = 0; i < ff_spectrum_.size(); i++) {

                    auto const& e = ff_spectrum_[i];
                    double const alpha = k0 * e.secT * s_field;
                    double const beta  = k0 * e.secTanT * Y;
                    double const gamma = k0 * e.sec2T * Z_dep;
                    double const env_exp = std::exp(gamma);
                    double const damping = (viscosity > 0.0) ? std::exp(damping_pre * e.sec4T) : 1.0;
                    double const w_trap = (i == 0 || i == ff_spectrum_.size() - 1) ? 0.5 : 1.0;
                    double const com = w_trap * env_exp * damping;

                    double const cA = std::cos(alpha), sA = std::sin(alpha);
                    double const cB = std::cos(beta),  sB = std::sin(beta);
                    double const PQ_inphase    = e.P * sA - e.Q * cA;  // matches elevation form
                    double const PQ_quadrature = e.P * cA + e.Q * sA;

                    // sec^5(theta) * sin(theta) = sec^3(theta) * secTanT
                    double const sec3T = e.sec2T * e.secT;
                    double const sec5T_sinT = sec3T * e.secTanT;
                    double const sec5T = e.sec4T * e.secT;

                    Sx += com * e.sec4T   * cB * PQ_inphase;
                    Sy += com * sec5T_sinT * sB * PQ_quadrature;
                    Sz += com * sec5T     * cB * PQ_quadrature;
                }

                double const A = (4.0 * speed * sq(k0) / M_PI) * dtheta;

                u_x_w = -A * fwd * Sx;
                u_y_w =  A       * Sy;
                u_z_w = -A       * Sz;
            }
        }

        double const ux = u_x_db + u_x_w;
        double const uy = u_y_db + u_y_w;
        double const uz = u_z_db + u_z_w;

        return {std::isnan(ux) ? 0.0 : ux,
                std::isnan(uy) ? 0.0 : uy,
                std::isnan(uz) ? 0.0 : uz};
    }

    /*
    glm::dvec3 get_velocity(double speed, double F, std::size_t floxi, std::size_t flozi) {

        auto const F2 = sq(F);
        auto const F4 = sq(F2);
        auto const F6 = F2 * F4;
        auto const linv = 1.0 / hull->get_length_wl();
        auto const wl = zs[0];
        auto const hp = std::abs(xs[1]-xs[0]) * linv * 1e-2;
        auto const hp2 = sq(hp);
        auto const hp3 = cub(hp);
        auto const hp4 = sq(hp2);
        auto const hp6 = sq(hp3);

//        auto const flo = glm::dvec3((xs[floxi] + xs[floxi+1])*0.5,
//        (get_beam(floxi,flozi) + get_beam(floxi+1,flozi) + get_beam(floxi,flozi+1) + get_beam(floxi+1,flozi+1)) * 0.25,
//        (zs[flozi] + zs[flozi+1])*0.5 - wl) * linv;

        auto const flo = glm::dvec3(xs[floxi]+hp, get_beam(floxi,flozi), std::min(zs[flozi] - wl - hp, -0.02)) * linv;

        glm::dvec3 local_vel(0,0,0);
        double sum_area = 0;

        auto const dS = std::abs((xs[1] - xs[0]) * (zs[1] - zs[0])) * sq(linv);

        for (std::size_t xi = 0; xi < xs.size(); xi++) {
            for (std::size_t zi = 0; zi < zs.size(); zi++) {

                if (get_beam(xi,zi) < 1e-4) {
                    continue;
                }

                auto src = glm::dvec3(xs[xi], get_beam(xi,zi), zs[xi] - wl) * linv;
                auto mirror = src;
                mirror.z *= -1;
                auto r = std::sqrt(sq(flo.x - src.x) + sq(flo.y - src.y) + sq(flo.z - src.z));
                auto r1 = std::sqrt(sq(flo.x - src.x) + sq(flo.y - src.y) + sq(flo.z + src.z));
                auto r1_sq = sq(r1);
                auto r1_cub = cub(r1);
                auto F2_r1 = F2 + r1;
                auto nazivnik = 1 / std::sqrt(r1_sq + hp4/(hp2 + r1_sq));
                auto ee = std::abs(flo.x - src.x) * nazivnik;
                auto eta = (flo.y - src.y) * nazivnik;
                auto zeta = -(flo.z + src.z) * nazivnik;
                auto etazeta = std::sqrt(sq(eta) + sq(zeta)) / (sq(eta) + sq(zeta) + hp4/(hp2 + sq(eta) + sq(zeta)));
                auto A_star = 4*F6 - 4*F4*r1 + 60*F2*r1_sq - 51*r1_cub;
                auto B_star = F6 + 74*F4*r1 - 189*F2*r1_sq + 48*r1_cub;
                auto C_star = 4*F6 - 14*F4*r1 + 3*F2*r1_sq - 15*r1_cub;
                auto Q = 2*F2*nazivnik / (sq(F2_r1)*(1+ee));
                auto P = 2/sq(F2_r1) + 4*F2*zeta/(cub(F2_r1)*(1+ee)) - 0.4*F2/sq(cub(F2_r1)) * (
                    (A_star + (B_star*zeta)*etazeta)*(1-ee) - (F2*C_star*ee)/(F2_r1)
                );
                auto D = 0.4*F2*(8*F6+13*F4*r1+37*F2*r1_sq+26*r1_cub)/sq(cub(F2_r1));
                auto E = 0.4*F2*(F4+39*F2*r1-24*r1_sq)/(cub(F2_r1)*sq(F2_r1));
                auto S = P + Q*zeta/(1+ee) - (D + E*zeta*etazeta)*ee;
                auto T = E/(1+ee) * eta*etazeta;
                auto RLx = sign(flo.x - src.x) * (P*ee + Q*zeta + (D/etazeta + E*zeta) / etazeta);
                auto RLy = S*eta + T*zeta;
                auto RLz = Q + T*eta - S*zeta;

                local_vel += ((flo - src) * (1/std::sqrt(cub(r) + hp6/(hp3+cub(r)))) - (flo - mirror) * (1/std::sqrt(cub(r1) + hp6/(hp3+cub(r1))))
                            + glm::dvec3(RLx, RLy, RLz) * double(get_strength(xi,zi)) // U_L_regularized
                            ) * dS / get_normal(xi,zi).y; // singular

                sum_area += dS;
            }
        }

//        local_vel *= (0.5*hull->get_wetted_area()*sq(linv))/sum_area; // fix area

        local_vel += get_normal(floxi, flozi) * (2 * M_PI * get_strength(floxi, flozi));

        local_vel *= (speed/(4*M_PI));
        local_vel.x = speed - local_vel.x;

        return local_vel;
    }
*/

public:

#ifdef USE_BEM
    glm::dvec3 get_viscous_pressure_resistance_bem(double speed, bool mirror_y) {

        auto N = xs.size() * zs.size();
        auto dx = xs[1] - xs[0];

        std::vector<PanelVector3*> panel_vertices(N);
        std::vector<Panel3*> panels;

        for (std::size_t xi = 0; xi < xs.size(); xi++) {
            for (std::size_t zi = 0; zi < zs.size(); zi++) {

                auto id = zi * xs.size() + xi;
                panel_vertices[id] = new PanelVector3(xs[xi], dx+get_beam(xi,zi), zs[zi]);
            }
        }

//        for (std::size_t xi = 0; xi < xs.size(); xi++) {
//            (*panel_vertices[xi])[1] = 0;
//        }

        for (std::size_t xi = 1; xi < xs.size(); xi++) {
            for (std::size_t zi = 1; zi < zs.size(); zi++) {

                auto id1 = (zi - 1) * xs.size() + xi - 1;
                auto id2 = (zi) * xs.size() + xi - 1;
                auto id3 = (zi) * xs.size() + xi;
                auto id4 = (zi - 1) * xs.size() + xi;

                if (std::abs(panel_vertices[id1]->y() + panel_vertices[id2]->y() + panel_vertices[id3]->y() + panel_vertices[id4]->y()) * 0.25 < 1e-3) {
                    continue;
                }

                Panel3* pnl = new Panel3();
                pnl->append_point(panel_vertices[id1]);
                pnl->append_point(panel_vertices[id2]);
                pnl->append_point(panel_vertices[id3]);
                pnl->append_point(panel_vertices[id4]);
                if (mirror_y) pnl->set_symmetric();
                panels.push_back(pnl);
            }
        }

        std::vector<PanelVector3*> trailing_edge;
//        for (std::size_t zi = 0; zi < zs.size(); zi++) {
//            trailing_edge.push_back(new PanelVector3(xs[xs.size()-1], 0, zs[zi]));
//        }

        DirichletDoublet0Source0Case3 bem(panels, trailing_edge);
        bem.v_inf = Vector3(speed, 0, 0);
        bem.A_ref = hull->get_wetted_area();
        bem.farfield = 1e5;
        auto res = bem.polar();

        /*
        Vector3 F();
        for (auto p : panels) {
            F += p->get_force();
        }
        F *= 0.5 * env->density * sq(speed);
        if (double_body) {
            F *= 0.5; // get rid of upper half
        }
*/
        // cleanup
        for (auto p : panels) {
            delete p;
        }

        for (auto pv: panel_vertices) {
            delete pv;
        }

        return {res.values[0].cD, res.values[0].cL, res.values[0].cS};
    }
#endif

};


#endif // THINSHIP_H
