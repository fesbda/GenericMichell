#ifndef PLANING_H
#define PLANING_H

#include "hull.h"
#include "interpolation.h"

// --- Savitsky root solves with un-bracketed fallbacks (SOLVER_FIX_PLAN.md H3) ---
// bisection() silently returns the far bracket end when the root is not bracketed. These wrappers
// detect the two physical corners where that happens and fall back sensibly. Both are bit-identical
// to the old inline bisection call whenever the root is properly bracketed (the common case).

/// Flat-plate lift coefficient from the deadrise lift coefficient: solve
/// C_L0 - 0.0065 β_deg C_L0^0.6 = C_Lβ. In the high-β / light-load corner the root can lie above
/// the initial upper bracket 4·C_Lβ, where bare bisection returns that end (under-predicting C_L0);
/// expand the upper bracket (×4, a few steps) until it brackets.
inline double savitsky_clo(double beta_deg, double C_Lbeta, bool* found = nullptr) {
    auto g = [=](double CL0) { return CL0 - 0.0065 * beta_deg * std::pow(std::max(CL0, 1e-12), 0.60); };
    double const lo = C_Lbeta;
    double hi = 4.0 * C_Lbeta + 1e-3;
    for (int k = 0; k < 8 && g(hi) < C_Lbeta; ++k) { hi = lo + 4.0 * (hi - lo); }
    return bisection(g, C_Lbeta, lo, hi, 1e-6, std::max(1e-8, 1e-4 * C_Lbeta), found);
}

/// Mean wetted-length ratio λ from the load coefficient: solve
/// C_L0 = τ^1.1 (0.012 √λ + 0.0055 λ^2.5 / Cv²) on [0.05, 15]. When the root is below the low end
/// (g(0.05) > C_L0), return the low end so the caller's 0.2 floor applies — rather than bare bisection
/// railing to the far end (which the caller would then min() to the λ cap).
inline double savitsky_lambda(double tau11, double cv2, double C_L0, double min_step, double abs_tol) {
    auto g = [=](double lam) { double const le = std::max(lam, 1e-3);
                               return tau11 * (0.012 * std::sqrt(le) + 0.0055 * std::pow(le, 2.5) / cv2); };
    double const target = std::max(C_L0, 1e-9);
    if (g(0.05) > target) { return 0.05; }
    return bisection(g, target, 0.05, 15.0, min_step, abs_tol);
}

/// Closed-form planing total-resistance from Savitsky's empirical equations,
/// with Blount-Fox hump multiplier. Returns the friction and wave-resistance
/// components separately so the caller can substitute its own friction model.
struct PlaningResistance {
    double friction = 0.0;
    double wave = 0.0;
};

struct TrueWettedProfile {
    std::vector<double> lambda_y;
    double mean = 0.0;
    double max = 0.0;
    int count = 0;

    bool valid() const { return count > 0; }
};

/// Savitsky-2012 warp-guard floor (rad) on the fore-aft deadrise spread. Prismatic hulls
/// measure up to ~1e-3 rad of pure section-read noise at the two fixed stations (a 1e-6
/// guard FLICKERS on them per-eval — observed as ±0.1 pt drift on the Fridsma/Taunton
/// suites); genuine warp starts near 0.09 rad (Begovic warp1). 0.02 sits two decades above
/// the noise and 5x below the mildest real warp.
static constexpr double S12_WARP_SPREAD_FLOOR = 0.02;

/// Static waterline wetted planform in beam units, sampled on half-beam buttocks.
///
/// Hull::get_wetted_buttock_profile returns each strip's local x-range. For a planing
/// planform the longitudinal coordinate must be measured from one common transom
/// reference, not from each strip's own aft-most intersection; otherwise flared/warped
/// sections lose the aft planform area and the pressure CoP/friction length shift.
inline TrueWettedProfile get_true_wetted_profile(std::shared_ptr<Hull> const& hull,
                                                 int n, double beam, double scale = 1.0) {
    TrueWettedProfile out;
    n = std::max(1, n);
    out.lambda_y.assign(std::size_t(n), 0.0);

    double const b = std::max(beam, 1e-6);
    double const fwd = double(hull->get_fwd());
    std::vector<double> const wp = hull->get_wetted_buttock_profile(n);

    bool any = false;
    double x_ref = (fwd > 0.0) ? std::numeric_limits<double>::max()
                               : std::numeric_limits<double>::lowest();
    for (int j = 0; j < n; j++) {
        double const xm = wp[2 * j], xM = wp[2 * j + 1];
        if (xM <= xm) {
            continue;
        }
        any = true;
        x_ref = (fwd > 0.0) ? std::min(x_ref, xm) : std::max(x_ref, xM);
    }
    if (!any) {
        return out;
    }

    double sum = 0.0;
    for (int j = 0; j < n; j++) {
        double const xm = wp[2 * j], xM = wp[2 * j + 1];
        if (xM <= xm) {
            continue;
        }
        double const x_forward = (fwd > 0.0) ? xM : xm;
        double const dist = (x_forward - x_ref) * fwd;
        double const lam = std::max(0.0, scale * dist / b);
        out.lambda_y[std::size_t(j)] = lam;
        out.max = std::max(out.max, lam);
        sum += lam;
        out.count++;
    }
    if (out.count > 0) {
        out.mean = sum / double(out.count);
    }
    return out;
}

/// Savitsky's planing resistance method (Savitsky 1964 + Blount & Fox 1976,
/// 1996). Treats the hull as a prismatic planing surface and computes lift,
/// drag, and friction in closed form. Mean β and chine beam are extracted from
/// the hull mesh (see Hull::get_deadrise, Hull::get_beam_chine).
class Savitsky {

private:
    std::shared_ptr<Hull> hull;
    std::shared_ptr<Environment> env;
    // ---- Tunable levers (grouped; mutated only via the set_* forwarders below) ----
    struct Config {
        double lambda_cap = 12.0;        // hard upper bound on the wetted λ (see Morabito)
        bool use_true_wetted = true;     // friction on the TRUE wetted length (physics-first default 2026-06-25)
        double friction_factor = 1.0;    // λ_f = friction_factor × geometric mean wetted length
        double planing_ca = 0.0006981;   // planing-bottom additive correlation allowance (planing practice).
                                         // Exposed only so the MaxSurf benchmark asymmetry (that tool run at
                                         // zero allowance) can be quantified by a Ca=0 A/B; default unchanged.
        // Speed-ramped spray-root friction: blend the true (calm-water) wetted length toward
        // Savitsky's spray-inclusive λ as Fn_vol rises through [lo,hi] (the spray-root wets the
        // bottom forward of the calm waterline, and grows with speed). Set hi<=lo to disable.
        // DEFAULT [5.5, 6.5] (adopted 2026-07-05, fixes_ab/taunton_datafix): the imposed-pose
        // probe on the corrected Taunton data proved the very-high-Fn deficit is the calm-plane
        // wetting collapsing at flat trim (the real wetted area is set by the spray-root
        // pile-up, which Savitsky's λ embeds); this window engages ABOVE every other suite's
        // range (Begović max Fn_vol 4.52, Fridsma < 4.5, displacement ~1 -> all byte-inert)
        // and takes Taunton 20.0 -> 14.8% (pooled 20.7 -> 15.4) with Series 62 9.4 -> 9.1%.
        double fric_ramp_lo = 5.5, fric_ramp_hi = 6.5;
        // Derived deadrise-scaled spray-root friction (default OFF). Blend the static true-wetted
        // length toward the Savitsky load-coefficient λ by the lift-efficiency deficit (1 - C_Lβ/C_L0):
        // a deadrise V-bottom lifts less per unit area, so it carries the load over a LONGER wetted
        // length -> more spray-root friction; the factor -> 0 at β=0 (flat plate) so it is do-no-harm
        // on flat hulls by construction. Parameter-free: the β-dependence is Savitsky's own lift law.
        bool deadrise_scaled_fric = false;
        // Savitsky (2012) warped-hull effective-attitude corrections (default OFF; see Morabito's
        // Config for the full note). In this class they re-express the running-trim friction
        // closure (the Auto path) at the effective deadrise + quarter-beam buttock trim; the
        // standalone predictive method is left as-is. Warp-guarded: constant-deadrise hulls take
        // the production path byte-identically.
        int s12_mode = 0;
        double s12_station_frac = 1.0;
        // pdyn Stage 3b friction-width correction (see the use_true_wetted block in calc).
        bool fric_width = false;
    } cfg;
    // s12 warp-guard decision, cached once per instance (see Morabito; pose-flicker hazard).
    int s12_guard_ = -1;

public:
    Savitsky(std::shared_ptr<Hull> ship_hull, std::shared_ptr<Environment> environment)
        : hull(ship_hull), env(environment) {}

    ~Savitsky() {}

    void set_lambda_cap(double cap) { cfg.lambda_cap = std::max(0.2, cap); }
    double get_lambda_cap() const { return cfg.lambda_cap; }

    /// Warp-resolved-trim productionization Step 2 (opt-in, default OFF): take the planing
    /// friction wetted length from the TRUE geometric per-strip profile (× friction_factor)
    /// instead of the Savitsky-inverted λ, so the friction is consistent with the
    /// true-planform pressure rather than the over-long load-coefficient λ.
    void set_true_wetted(bool on) { cfg.use_true_wetted = on; }
    void set_friction_factor(double f) { cfg.friction_factor = std::max(0.1, f); }
    void set_friction_cv_ramp(double lo, double hi) { cfg.fric_ramp_lo = lo; cfg.fric_ramp_hi = hi; }
    void set_friction_width(bool on) { cfg.fric_width = on; }
    void set_planing_ca(double ca) { cfg.planing_ca = std::max(0.0, ca); }
    void set_deadrise_scaled_friction(bool on) { cfg.deadrise_scaled_fric = on; }
    /// Savitsky (2012) effective-attitude corrections (running-trim friction path; see Morabito).
    void set_savitsky2012(int mode, double station_frac) {
        cfg.s12_mode = std::max(0, std::min(2, mode));
        cfg.s12_station_frac = std::max(0.1, std::min(2.0, station_frac));
        s12_guard_ = -1; }

    /// use_running_trim: evaluate friction/induced at the hull's actual (equilibrium)
    /// trim instead of the predictive Savitsky trim. The Auto blend uses this so the
    /// friction is consistent with the converged attitude (and the Morabito pressure);
    /// the standalone Savitsky method keeps the predictive 2-DOF trim (default).
    PlaningResistance calc(double speed, bool use_running_trim = false) {

        PlaningResistance out;

        if (speed <= 1e-2) {
            return out;
        }

        // Use the fixed design displacement, not the running one: a planing hull
        // rises out of the water as it speeds up, so get_displacement() collapses,
        // but the weight the planing surface must support stays constant.
        double const vol = hull->get_reference_displacement();
        double const vol13 = std::pow(std::max(vol, 1e-12), 1.0/3.0);
        double const mass = vol * env->get_density();
        double const g = env->get_gravity();
        double const chine_beam = std::max(hull->get_beam_chine(), 1e-3);
        // beta/tau re-expressed at the Savitsky-2012 effective attitude below when active.
        double beta = hull->get_deadrise();                // radians, in [0, 40°]
        double beta_deg = degrees(beta);
        double const Fn_vol = hull->get_fn(speed, vol13);
        double const Fn_beam = hull->get_fn(speed, chine_beam);

        // LCG expressed in beam units forward of the transom (Savitsky's convention).
        double const fwd = double(hull->get_fwd());
        double const x_cg = hull->get_cg(true).x;
        double const x_transom = (fwd > 0) ? hull->get_min(false).x : hull->get_max(false).x;
        double const LCG = std::max(0.1, (x_cg - x_transom) * fwd / chine_beam);

        // Lift coefficient of deadrise planing surface.
        double const C_Lbeta = g * vol / (0.5 * sq(speed) * sq(chine_beam));

        // Equivalent flat-plate lift coefficient (Savitsky's iterative correction).
        // Solve C_L0 - 0.0065 β_deg C_L0^0.60 - C_Lbeta = 0.
        bool clo_found = false;
        double const C_L0 = savitsky_clo(beta_deg, C_Lbeta, &clo_found);

        double const Fnb2 = std::max(1e-6, sq(Fn_beam));

        // Trim: the actual running (equilibrium) trim when use_running_trim — so the
        // Auto blend's friction is consistent with the converged attitude — otherwise the
        // predictive Savitsky 2-DOF trim where the centre of pressure matches the LCG
        // (Savitsky 1964 Eq 31: C_p = λ (0.75 - 1 / (5.21 (Cv/λ)^2 + 2.39))).
        double tau_deg;
        if (use_running_trim) {
            tau_deg = std::max(0.25, -degrees(hull->get_pitch()));
        } else {
            bool tau_found = false;
            tau_deg = bisection(
                [=](double guess_tau_deg) {
                    double const t11 = std::pow(std::max(guess_tau_deg, 0.1), 1.1);
                    double const lam = savitsky_lambda(t11, Fnb2, C_L0, 1e-6, 1e-5);
                    double const lam_safe = std::max(lam, 0.01);
                    double const cp_over_lam = 0.75 - 1.0 / (5.21 * (Fnb2 / sq(lam_safe)) + 2.39);
                    return lam_safe * cp_over_lam;  // C_p in beam units forward of transom
                },
                LCG, 0.5, 15.0, 1e-4, 1e-4, &tau_found);
            tau_deg = std::max(0.25, tau_deg);  // floor for well-posedness
            // Mirror find_equilibrium's predict_trim guard (resistance.h): the CoP=LCG solve returns
            // its far bracket end (~15 deg) when there is no physical root; don't trust a >10 deg
            // predictive trim, cap the un-bracketed rail at the validity ceiling.
            if (!tau_found && tau_deg > 10.0) { tau_deg = 10.0; }
        }
        double tau = radians(tau_deg);

        // Mean wetted-length-to-beam ratio λ from Savitsky's lift equation inverted
        // at the *running* trim: C_L0 = τ^1.1 (0.012 λ^0.5 + 0.0055 λ^2.5 / Cv^2).
        // This load-coefficient λ is independent of LCG, unlike the centre-of-pressure
        // relation it replaces, so it does not collapse for the aft LCG of warped
        // hulls (where CoP-λ under-predicts the wetted length and hence the friction).
        double const tau11 = std::pow(tau_deg, 1.1);
        double const lambda = savitsky_lambda(tau11, Fnb2, C_L0, 1e-7, 1e-6);

        double const lam_sav = std::min(cfg.lambda_cap, std::max(lambda, 0.2));
        double lam_eff = lam_sav;
        // POSE COHERENCE (savblend Stage -1 validation, 2026-07-03): the true-wetted profile is
        // read from the mesh AS POSED. In the Auto path (use_running_trim) the hull sits at the
        // converged planing pose, so λ_geom is the physical wetted length. In the STANDALONE
        // predictive method the hull is at design float — reading the static profile there pairs
        // full-length friction with predictive planing trim (incoherent; Begović standalone
        // over-predicted 1.5–3.4×, growing with v² on the mono, while MaxSurf's Savitsky tracks
        // EFD ±3–16%; fixes_ab/savblend/). Predictive mode therefore keeps Savitsky's own
        // load-coefficient λ, as published.
        if (cfg.use_true_wetted && use_running_trim) {
            // True per-strip wetted length from the common transom reference
            // (posed mesh) -> mean λ_geom × friction_factor.
            TrueWettedProfile const tw = get_true_wetted_profile(hull, 50, chine_beam, cfg.friction_factor);
            if (tw.valid()) {
                lam_eff = std::max(0.2, tw.mean);
                // Friction-width correction (pdyn Stage 3b, DEFAULT OFF): the friction area
                // below is lam_eff * b^2, i.e. the mean strip length times the FULL chine
                // beam as the wetted width. The DYNAMIC wetted width is the Wagner pile-up
                // spread of the calm-plane width, capped at the chine:
                //     w_dyn = min(b, (pi/2) * beam_wl(pose)).
                // Chines-wet hulls clamp to 1 and are untouched (GPPH: the EFD wet-chine
                // column and the CFD shear both certify full-width friction; the raw
                // calm-plane ratio without the pi/2 factor wrongly shrank it, 3.2 -> 7.4%).
                // Genuinely chines-dry states (slender prisms at very high CvB: Taunton,
                // warp3 top band) keep a sub-1 width -- the physics the legacy pose-collapsed
                // beam readout was accidentally supplying, whose removal by robust_chine blew
                // those bands up +18..27%. Applied ONLY to the geometric lam -- the Savitsky
                // load-coefficient lam keeps its own published b^2 convention.
                if (cfg.fric_width) {
                    double const w_wet = hull->get_beam_wl();
                    if (w_wet > 1e-6) {
                        double const w_dyn = 1.5707963267948966 * w_wet;   // Wagner pi/2
                        lam_eff *= std::max(0.05, std::min(1.0, w_dyn / chine_beam));
                        lam_eff = std::max(0.2, lam_eff);
                    }
                }
            }
        }
        // Derived deadrise-scaled spray-root friction (default off): blend the static true-wetted
        // length toward the Savitsky load λ by the lift-efficiency deficit (1 - C_Lβ/C_L0). f_beta=0
        // at β=0 (flat, do-no-harm) and grows with deadrise. Only meaningful when use_true_wetted
        // (otherwise lam_eff is already lam_sav).
        if (cfg.deadrise_scaled_fric && cfg.use_true_wetted) {
            double const f_beta = std::max(0.0, std::min(1.0, 1.0 - C_Lbeta / std::max(C_L0, 1e-9)));
            lam_eff = lam_eff + f_beta * (lam_sav - lam_eff);
        }
        // Speed-ramped spray-root friction (default off): as Fn_vol rises through
        // [fric_ramp_lo, fric_ramp_hi], blend the true (calm-water) wetted length toward the
        // spray-inclusive Savitsky λ -- the spray-root wets the bottom forward of the static
        // waterline and that wetted area grows with speed (the high-Fn friction deficit).
        if (cfg.fric_ramp_hi > cfg.fric_ramp_lo) {
            double t = (Fn_vol - cfg.fric_ramp_lo) / (cfg.fric_ramp_hi - cfg.fric_ramp_lo);
            t = std::max(0.0, std::min(1.0, t));
            double const g = t * t * (3.0 - 2.0 * t);   // smoothstep
            lam_eff = (1.0 - g) * lam_eff + g * lam_sav;
        }

        // Savitsky (2012) effective attitude (default OFF; running-trim/Auto path only): evaluate
        // the friction closure below at the effective deadrise + quarter-beam buttock trim sampled
        // at the mean-wetted-length station, consistent with the Morabito pressure kernel. The
        // wetted length keeps the geometric lam_eff. Warp-guarded => prisms byte-identical.
        if (cfg.s12_mode > 0 && use_running_trim) {
            if (s12_guard_ < 0) {   // decided once per instance (pose-flicker hazard; see Morabito)
                double const Lwl_s = std::max(1e-3, hull->get_length_wl());
                double const bA_s = hull->get_deadrise_at_x(x_transom + fwd * 0.15 * Lwl_s);
                double const bF_s = hull->get_deadrise_at_x(x_transom + fwd * 0.60 * Lwl_s);
                s12_guard_ = (bA_s > 0.0 && bF_s > 0.0 &&
                              std::abs(bF_s - bA_s) > S12_WARP_SPREAD_FLOOR) ? 1 : 0;
            }
            if (s12_guard_ == 1) {
                double const x_s = x_transom
                                   + fwd * cfg.s12_station_frac * std::max(lam_eff, 0.2) * chine_beam;
                if (cfg.s12_mode == 1) {
                    double const bl_s = hull->get_deadrise_about(x_s, 0.25 * chine_beam);
                    if (bl_s > 0.0) {
                        beta = std::max(1e-3, std::min(bl_s, radians(50.0)));
                        beta_deg = degrees(beta);
                    }
                }
                double const delta = hull->get_quarter_buttock_delta(x_s, 0.25 * chine_beam);
                tau_deg = std::max(0.25, tau_deg + degrees(delta));
                tau = radians(tau_deg);
            }
        }

        // Mean velocity over the wetted bottom (Savitsky's velocity reduction).
        double const cl0_corr = 0.012 * std::sqrt(lam_eff) * std::pow(tau_deg, 1.1);
        double const beta_corr = 0.0065 * beta_deg * std::pow(std::max(cl0_corr, 1e-12), 0.6);
        double const speed_sq_ratio = std::max(0.05,
            1.0 - (cl0_corr - beta_corr) / (lam_eff * std::cos(tau)));
        double const speed_mean = speed * std::sqrt(speed_sq_ratio);

        double const length_wetted = lam_eff * chine_beam;
        double const Rn = std::max(1.0, length_wetted * speed_mean / env->get_viscosity());
        double const Cf = 0.075 / sq(std::log10(Rn) - 2.0);
        double const Ca = cfg.planing_ca;  // additive correlation allowance (planing practice; settable)

        double const cos_beta = std::cos(beta);
        double const cos_tau = std::cos(tau);
        double const Rf = 0.5 * env->get_density() * sq(speed_mean) * lam_eff * sq(chine_beam)
                          * (Cf + Ca) / std::max(cos_beta * cos_tau, 1e-6);

        // Wave / induced-drag component from the trim angle.
        double const Rw = mass * g * std::tan(tau);

        // Blount-Fox hump multiplier (recommended in Blount's High Speed Propulsion Design).
        // Applied symmetrically to both components so the proportion is preserved.
        double const LCB_over_b = LCG;
        double const M = 0.98
                         + 2.0 * std::pow(LCB_over_b, 1.45) * std::exp(-2.0 * (Fn_vol - 0.85))
                         - 3.0 * LCB_over_b * std::exp(-3.0 * (Fn_vol - 0.85));
        double const MP = 1.0 + std::max(0.0, 0.5 * (M - 1.0));

        out.friction = Rf * MP;
        out.wave = Rw * MP;
        return out;
    }

    /// Convenience: total Savitsky resistance.
    double calc_total(double speed) {
        auto r = calc(speed);
        return r.friction + r.wave;
    }

    /// Savitsky's predictive 2-DOF running trim (deg): the trim at which the
    /// centre of pressure coincides with the LCG (Savitsky 1964, Eq. 31), with the
    /// wetted length λ(τ) from the lift equation. This is the coupled λ–τ planing
    /// attitude — it falls with speed as the data require — and is used by the Auto
    /// blend to set the planing trim, since the bottom-pressure moment alone puts the
    /// centre of pressure too far forward at high Cv and over-predicts the trim.
    /// Evaluated against the design CG (call with the hull at the reference attitude).
    /// beta_override (rad, >=0) replaces the hull's area-mean deadrise — pass the
    /// effective (wetted-bottom, section-local) planing deadrise so the predictive trim
    /// of a warped hull is referenced to the lower aft deadrise it actually planes on,
    /// not the area mean (which over-predicts the trim). <0 keeps the area mean.
    double predict_trim(double speed, double beta_override = -1.0) {
        if (speed <= 1e-2) return 0.0;
        double const vol = hull->get_reference_displacement();
        double const g = env->get_gravity();
        double const b = std::max(hull->get_beam_chine(), 1e-3);
        double const Fnb2 = std::max(1e-6, sq(hull->get_fn(speed, b)));
        double const C_Lbeta = g * vol / (0.5 * sq(speed) * sq(b));
        double const fwd = double(hull->get_fwd());
        double const x_cg = hull->get_cg(true).x;
        double const x_transom = (fwd > 0) ? hull->get_min(false).x : hull->get_max(false).x;
        double const LCG = std::max(0.1, (x_cg - x_transom) * fwd / b);

        // Savitsky (2012): the CoP=LCG relation solved below is prismatic, so on a warped hull
        // it yields the EFFECTIVE (quarter-beam buttock) trim at the effective deadrise; the
        // keel pitch reported to the equilibrium blend is tau_eff - delta. Warp-guarded.
        bool s12 = false;
        if (cfg.s12_mode > 0) {
            if (s12_guard_ < 0) {   // decided once per instance (pose-flicker hazard; see Morabito)
                double const Lwl_s = std::max(1e-3, hull->get_length_wl());
                double const bA_s = hull->get_deadrise_at_x(x_transom + fwd * 0.15 * Lwl_s);
                double const bF_s = hull->get_deadrise_at_x(x_transom + fwd * 0.60 * Lwl_s);
                s12_guard_ = (bA_s > 0.0 && bF_s > 0.0 &&
                              std::abs(bF_s - bA_s) > S12_WARP_SPREAD_FLOOR) ? 1 : 0;
            }
            s12 = (s12_guard_ == 1);
        }

        double beta_deg = degrees((beta_override >= 0.0) ? beta_override : hull->get_deadrise());
        double tau_deg = 0.25, delta_s = 0.0;
        int const passes = s12 ? 2 : 1;   // pass 0 locates the station, pass 1 re-solves there
        for (int pass = 0; pass < passes; pass++) {
            double const C_L0 = savitsky_clo(beta_deg, C_Lbeta);
            bool found = false;
            tau_deg = bisection(
                [=](double guess_tau_deg) {
                    double const t11 = std::pow(std::max(guess_tau_deg, 0.1), 1.1);
                    double const lam = savitsky_lambda(t11, Fnb2, C_L0, 1e-6, 1e-5);
                    double const ls = std::max(lam, 0.01);
                    return ls * (0.75 - 1.0 / (5.21 * (Fnb2 / sq(ls)) + 2.39));
                },
                LCG, 0.5, 15.0, 1e-4, 1e-4, &found);
            tau_deg = std::max(0.25, tau_deg);
            if (!s12 || pass == passes - 1) break;
            double const t11 = std::pow(std::max(tau_deg, 0.1), 1.1);
            double const lam0 = std::max(0.2, savitsky_lambda(t11, Fnb2, C_L0, 1e-6, 1e-5));
            double const x_s = x_transom + fwd * cfg.s12_station_frac * lam0 * b;
            // mode 1 resamples the effective deadrise (unless the caller pinned one); mode 2
            // keeps the caller's/mean deadrise and applies only the trim conversion.
            if (cfg.s12_mode == 1 && beta_override < 0.0) {
                double const bl_s = hull->get_deadrise_about(x_s, 0.25 * b);
                if (bl_s > 0.0) {
                    beta_deg = degrees(std::max(1e-3, std::min(bl_s, radians(50.0))));
                }
            }
            delta_s = hull->get_quarter_buttock_delta(x_s, 0.25 * b);
        }
        if (s12) {
            tau_deg = std::max(0.25, tau_deg - degrees(delta_s));
        }
        return tau_deg;
    }

};


/// Morabito (2014) empirical bottom-pressure distribution for prismatic
/// planing hulls. Integrates the pressure field over (X, Y) in non-dimensional
/// beam units to produce (drag, lift, torque about CG) at a given speed.
/// Matches LocalFlow::get_drag_lift_torque_noblesse so it slots directly into
/// Resistance equilibrium iteration and pressure-aware dispatch.
class Morabito {

private:
    std::shared_ptr<Hull> hull;
    std::shared_ptr<Environment> env;
    // ---- Tunable levers (grouped; mutated only via the set_* forwarders below) ----
    struct Config {
        int nx = 200;
        int ny = 50;
        double lambda_cap = 12.0;   // hard upper bound on the Savitsky-inverted wetted λ
        bool use_true_planform = true;    // warp-resolved true wetted planform (physics-first default 2026-06-25)
        double pileup = 1.0;              // manual spray-root pile-up factor
        bool auto_pileup = false;         // compute the pile-up f(β_aft,τ) automatically
        bool use_hydrostatic = false;     // Pstat lift+moment handed to hull buoyancy (keep drag); physics-first default 2026-06-25
        double dynamic_cop_fraction = 0.70;  // dynamic resultant at frac*lambda fwd of transom (just aft of the
                                             // Wagner/Savitsky 3/4-length CoP). 0.70 selected 2026-07-01 as the
                                             // do-no-harm sweep optimum (improves all 4 series; see sweep_results.db
                                             // run 4). Was 0.75. <0 = integral CoP.
        bool cop_savitsky = false;           // true = use Savitsky 1964 Eq31 speed-dependent Cp/lambda for the
                                             // dynamic CoP fraction instead of the constant 0.75 high-Cv asymptote.
        // Cv-ramp of the dynamic CoP fraction, GATED by the hull's fore-aft deadrise spread (warp).
        // frac_eff = frac + cop_cv_slope * warp * (cop_cv_ref - Cv), clamped. Default slope 0 = OFF
        // (production byte-identical). The warp gate (fwd-aft deadrise, rad) is ~0 on prismatic hulls
        // (Series62/Fridsma/monohedral) so the ramp activates only on warped forms. Physically: the
        // constant 3/4-length CoP over-trims warped hulls at high Cv and under-trims them through the
        // hump; a positive slope moves the CoP forward (more bow-up) at low Cv and aft at high Cv,
        // matching the measured monotonic trim fall. See diag_wave_climb_lead.py / ATTITUDE notes.
        double cop_cv_slope = 0.0;    // 0 = OFF (byte-identical)
        double cop_cv_ref = 2.0;      // neutral Cv where the warp ramp crosses zero
        double cop_warp_floor = 0.13; // rad; fore-aft deadrise spread below this -> ramp inert.
                                      // Excludes prismatic (Fridsma) AND mildly-warped hulls whose
                                      // literature 3/4 CoP already validates (Series62 spread ~0.09,
                                      // warp1 ~0.09); activates only the strongly-warped warp2/warp3
                                      // (spread ~0.19/0.27). The gate samples FIXED Lwl-fractions so it
                                      // is a stable hull-shape property, not a running quantity.
        // Warped-transition lift closure (default-OFF). On warped hulls only, through a Cv transition
        // window, raise the aft (low-deadrise) wetted area / spray-root pile-up so the dynamic LIFT and
        // the induced DRAG rise together, with an INDEPENDENT control on the bow-up MOMENT. Unlike the
        // CoP ramp above (which relocates lift, is anti-coupled with sinkage, and never moves the sinkage
        // phase). The lift and moment amplitudes are DECOUPLED: lift fixes the squat phase; moment is the
        // CoP shift toward the enlarged planform's integral centroid, dialled separately so the squat fix
        // does not drag the trim past measured through the upper hump. Shares the cop_warp_floor hinge, so
        // it stays inert on prismatic / mildly-warped hulls. Both amps 0 (default) = OFF (byte-identical).
        double warp_lift_amp = 0.0;     // LIFT amplitude (Hook 1: pile-up/wetted-area boost). 0 = OFF
        double warp_moment_amp = 0.0;   // MOMENT amplitude (Hook 2: CoP blend toward integral centroid). 0 = OFF
        double warp_lift_cv_lo = 1.0;   // beam-Froude window lower Cv (ramp up over [lo, lo+0.3])
        double warp_lift_cv_hi = 2.4;   // window upper Cv (ramp down over [hi-0.4, hi]); dies in full planing
        // Savitsky (2012) warped-hull effective-attitude corrections (default OFF). The 3rd
        // Chesapeake Power Boat Symposium result: prismatic planing equations remain valid on a
        // warped bottom if the EFFECTIVE deadrise (local deadrise at the mean-wetted-length
        // station) and the EFFECTIVE trim (quarter-beam buttock angle there, NOT the keel angle)
        // are substituted. Mode 1 = as published (constant beta_eff replaces the local-beta
        // profile, tau_eff everywhere); mode 2 = tau-only hybrid (keep the kernel's own per-cell
        // local-beta resolution, substitute tau_eff only). Guarded by the fore-aft deadrise
        // spread, so constant-deadrise hulls take the production path byte-identically.
        int s12_mode = 0;               // 0 = OFF (production byte-identical)
        double s12_station_frac = 1.0;  // station = transom + frac*lambda*b ("forward edge of
                                        // the mean wetted length" = 1.0; settable pending PDF)
    } cfg;
    // s12 warp-guard decision, cached once per instance (-1 = undecided). NOT config: state.
    int s12_guard_ = -1;

public:
    Morabito(std::shared_ptr<Hull> ship_hull, std::shared_ptr<Environment> environment)
        : hull(ship_hull), env(environment) {}

    ~Morabito() {}

    /// Morabito recommends 200×50; bump nx for very low deadrise (sharp stagnation peak).
    void set_grid(int new_nx, int new_ny) {
        cfg.nx = std::max(20, new_nx);
        cfg.ny = std::max(10, new_ny);
    }

    int get_nx() const { return cfg.nx; }
    int get_ny() const { return cfg.ny; }

    /// Cap the Savitsky-inverted mean wetted length-to-beam ratio λ at a physical bound.
    /// The unconstrained inversion extrapolates to λ > L/b (longer than the hull) at light
    /// load / high deadrise; a physical cap (e.g. L/b or the geometric wetted length) keeps
    /// the bottom-pressure planform on the actual hull. Set ≥12 to disable.
    void set_lambda_cap(double cap) { cfg.lambda_cap = std::max(0.2, cap); }
    double get_lambda_cap() const { return cfg.lambda_cap; }

    /// Warp-resolved-trim productionization (opt-in, default OFF so production is unchanged):
    /// integrate the bottom pressure over the TRUE per-strip wetted planform (geometric wetted
    /// length from the posed mesh × a spray-root pile-up factor) instead of the Savitsky-λ
    /// triangle. The pile-up bridges the static-intersection length to the physical one;
    /// f∈[1,2.3], ~1 at low (warped-aft) deadrise. See proto_phase*.py for the characterisation.
    void set_true_planform(bool on) { cfg.use_true_planform = on; }
    void set_pileup(double f) { cfg.pileup = std::max(0.1, f); }
    void set_auto_pileup(bool on) { cfg.auto_pileup = on; }
    /// Drop the hydrostatic (Pstat) part of the bottom pressure. Default ON (keep it). When the
    /// hull's own buoyancy moment is in the force balance, the planing Pstat double-counts the
    /// bottom hydrostatic — set OFF to let the hull buoyancy carry it (Stage-C CoP experiment).
    void set_hydrostatic(bool on) { cfg.use_hydrostatic = on; }
    /// Place the DYNAMIC planing resultant at frac*lambda forward of the transom (the Savitsky/Wagner
    /// centre of pressure that marches forward with Cv), instead of the integral's own dynamic
    /// centroid (which sits too far aft at high Cv). <0 = off (default). Stage-C2 CoP experiment.
    void set_dynamic_cop_fraction(double frac) { cfg.dynamic_cop_fraction = frac; }
    /// Use Savitsky's speed-dependent CoP (Eq31 Cp/lambda) for the dynamic resultant instead of the
    /// constant 0.75 asymptote. Parameter-free (cited 1964 relation). off (default) = byte-identical.
    void set_cop_savitsky(bool on) { cfg.cop_savitsky = on; }
    /// Cv-ramp the dynamic-CoP fraction, gated by the hull's fore-aft deadrise spread (warp).
    /// slope=0 (default) is OFF and leaves the constant-fraction CoP byte-identical. ref is the
    /// neutral Cv (beam Froude V/sqrt(g*b)) where the ramp crosses zero. Warped hulls only.
    void set_cop_cv_ramp(double slope, double ref, double warp_floor) {
        cfg.cop_cv_slope = slope; cfg.cop_cv_ref = std::max(0.1, ref); cfg.cop_warp_floor = std::max(0.0, warp_floor); }

    /// Warped-transition lift closure (Stage-B prototype). Raise the aft wetted area / spray-root
    /// pile-up through the Cv transition window on warped hulls so lift + induced drag move together,
    /// with the bow-up moment dialled INDEPENDENTLY: lift_amp drives Hook 1 (the squat/phase fix),
    /// moment_amp drives Hook 2 (the CoP blend toward the integral centroid). Both 0 (default) OFF.
    /// cv_lo/cv_hi bound the beam-Froude window; it dies out below the hump and in full planing. Warp
    /// gate shared with the CoP ramp (cop_warp_floor), so prismatic / mildly-warped hulls stay inert.
    void set_warp_lift_closure(double lift_amp, double moment_amp, double cv_lo, double cv_hi) {
        cfg.warp_lift_amp = std::max(0.0, lift_amp);
        cfg.warp_moment_amp = std::max(0.0, moment_amp);
        cfg.warp_lift_cv_lo = std::max(0.1, cv_lo);
        cfg.warp_lift_cv_hi = std::max(cfg.warp_lift_cv_lo + 1e-3, cv_hi); }

    /// Savitsky (2012) warped-hull effective-attitude corrections. mode 0 = OFF (default,
    /// byte-identical), 1 = as-published equivalent prism (beta_eff + tau_eff substituted into
    /// the prismatic relations), 2 = tau-only hybrid (local-beta profile kept). station_frac
    /// places the sampling station at transom + frac*lambda*b.
    void set_savitsky2012(int mode, double station_frac) {
        cfg.s12_mode = std::max(0, std::min(2, mode));
        cfg.s12_station_frac = std::max(0.1, std::min(2.0, station_frac));
        s12_guard_ = -1; }
    int get_savitsky2012_mode() const { return cfg.s12_mode; }

    /// Geometric mean wetted length-to-beam ratio from the keel/chine waterline
    /// intersections (the value the Savitsky-λ replaced because it collapses as the hull
    /// rises). Exposed to bracket the truth against the over-long Savitsky λ.
    double get_geometric_lambda() {
        return compute_lambda(std::max(hull->get_beam_chine(), 1e-3));
    }

    /// Returns {drag, lift, torque-about-CG} in N, N, N·m. If out_pressure_field
    /// is non-null, fills a flattened ny × nx (X varying fastest) array of P in Pa.
    glm::dvec3 get_drag_lift_torque(double speed, std::vector<double>* out_pressure_field = nullptr,
                                    double* out_lambda = nullptr, double* out_tau = nullptr,
                                    glm::dvec4* out_split = nullptr) {

        if (out_pressure_field) {
            out_pressure_field->assign(std::size_t(cfg.nx) * std::size_t(cfg.ny), 0.0);
        }

        if (speed <= 1e-2) {
            return {0.0, 0.0, 0.0};
        }

        double const rho = env->get_density();
        double const g = env->get_gravity();
        double const V = speed;
        double const q = 0.5 * rho * sq(V);
        double const b = std::max(hull->get_beam_chine(), 1e-3);
        double const Cv = V / std::sqrt(g * b);
        double const Cv2 = std::max(sq(Cv), 1e-12);

        // β: trim is the pitch angle (sign convention matches find_equilibrium).
        // beta/tau and their derived values are re-expressed at the Savitsky-2012 effective
        // attitude below (s12_mode active + warped hull); the default path leaves them untouched.
        double beta = std::max(1e-3, std::min(hull->get_deadrise(), radians(50.0)));
        double const tau_raw = -hull->get_pitch();
        // Morabito's equations are derived for positive bow-up trim. Clamp to validity range.
        double tau = std::max(radians(0.05), std::min(tau_raw, radians(30.0)));
        double tau_deg = degrees(tau);
        double sin_tau = std::sin(tau);
        double cos_tau = std::cos(tau);

        // λ = mean wetted-length-to-beam ratio from Savitsky's lift equation inverted
        // at the running trim (the load-coefficient λ, consistent with the friction
        // term). This replaces the geometric keel/chine wetted length, which collapses
        // as the hull rises — taking the bottom pressure, lift and CG moment down with
        // it, and pinning the equilibrium trim too flat. With the load-coefficient λ the
        // pressure spans the correct planform, so the centre of pressure sits forward
        // enough to support the planing trim.
        double const vol = hull->get_reference_displacement();
        double const C_Lbeta = g * vol / (0.5 * sq(V) * sq(b));
        double const beta_deg_l = degrees(beta);
        double const C_L0 = savitsky_clo(beta_deg_l, C_Lbeta);
        double lambda = std::min(cfg.lambda_cap, std::max(0.2,
            savitsky_lambda(std::pow(tau_deg, 1.1), Cv2, C_L0, 1e-7, 1e-6)));

        // Warped-transition lift closure (default-OFF): build the shared gate g_tr = warp × Cv-window.
        // warp reuses the fixed-Lwl-fraction fore-aft deadrise spread + cop_warp_floor hinge (zero on
        // prismatic/mildly-warped hulls); the Cv window is a bump that rises over [cv_lo, cv_lo+0.3],
        // holds through the semi-planing transition, and dies over [cv_hi-0.4, cv_hi] so it does not
        // touch full-planing behaviour. Both hooks (lift+drag below, moment at the CoP block) use it.
        double warp_lift_gate = 0.0;
        if (cfg.warp_lift_amp > 0.0 || cfg.warp_moment_amp > 0.0) {
            auto const smooth01 = [](double t) {
                t = std::min(1.0, std::max(0.0, t)); return t * t * (3.0 - 2.0 * t); };
            double const fwd_w = double(hull->get_fwd());
            double const xT_w = (fwd_w > 0) ? hull->get_min(false).x : hull->get_max(false).x;
            double const Lwl_w = std::max(1e-3, hull->get_length_wl());
            double const bA_w = hull->get_deadrise_at_x(xT_w + fwd_w * 0.15 * Lwl_w);
            double const bF_w = hull->get_deadrise_at_x(xT_w + fwd_w * 0.60 * Lwl_w);
            double warp_w = (bF_w > 0.0 && bA_w > 0.0) ? (bF_w - bA_w) : 0.0;
            warp_w = std::max(0.0, warp_w - cfg.cop_warp_floor);
            double const up = smooth01((Cv - cfg.warp_lift_cv_lo) / 0.3);
            double const dn = 1.0 - smooth01((Cv - (cfg.warp_lift_cv_hi - 0.4)) / 0.4);
            warp_lift_gate = warp_w * up * dn;          // rad × [0,1]; 0 unless warped AND in transition
        }

        // True-planform mode (opt-in): replace the Savitsky-λ triangle with the geometric
        // per-strip wetted length (posed mesh) × the pile-up factor. beta_scale is the forward
        // extent used to sample the deadrise profile so it spans the actual wetted planform.
        std::vector<double> true_lamY;
        double beta_scale = lambda;
        double lambda_used = lambda;
        if (cfg.use_true_planform) {
            double const fwd_t = double(hull->get_fwd());
            double const xT_t = (fwd_t > 0) ? hull->get_min(false).x : hull->get_max(false).x;
            double f_pile = cfg.pileup;
            if (cfg.auto_pileup) {
                // Warp-resolved pile-up f(β_aft, τ), fit to proto_phase3_pileup.py (β_aft = the
                // deadrise ~0.3 beam fwd of the transom, the low aft deadrise on a warped hull).
                double const ba = degrees(hull->get_deadrise_at_x(xT_t + 0.3 * b * fwd_t));
                double const b_aft = (ba > 0.0) ? ba : degrees(beta);
                f_pile = std::min(3.0, std::max(0.5, 0.1624 + 0.0589 * b_aft + 0.2221 * tau_deg));
            }
            // Hook 1 (warped-transition lift closure): extend the spray-root pile-up forward through
            // the transition so the enlarged planform raises dynamic lift AND its induced drag
            // together. f_pile unchanged when OFF (warp_lift_gate = 0) ⇒ byte-identical.
            f_pile *= (1.0 + cfg.warp_lift_amp * warp_lift_gate);
            TrueWettedProfile const tw = get_true_wetted_profile(hull, cfg.ny, b, f_pile);
            if (tw.valid()) {
                true_lamY = tw.lambda_y;
                beta_scale = std::max(tw.max, 1e-3);
                lambda_used = std::max(0.2, tw.mean);
            } else {
                true_lamY.assign(cfg.ny, 0.0);
                beta_scale = 1e-3;
                lambda_used = 0.0;
            }
        }
        // Savitsky (2012) effective attitude (default OFF): on a warped bottom the prismatic
        // planing relations stay valid when evaluated at the EFFECTIVE deadrise (local deadrise
        // at the mean-wetted-length station) and the EFFECTIVE trim (the quarter-beam buttock
        // angle there — on a warped hull the aft body runs steeper than the keel pitch reports,
        // so the same pressure resolves into MORE drag per unit lift). Guard: fore-aft deadrise
        // spread ~0 => constant-deadrise hull => production path, byte-identical. The guard is
        // decided ONCE per instance and cached (s12_guard_): a per-eval read of the posed
        // geometry can transiently mis-read the spread at extreme equilibrium-scan poses (the
        // posed-Lwl station lands on the stem rake), flipping the toggle for a single residual
        // eval and shifting the root — the Radojcic pose-dependence lesson. The decision has
        // two decades of margin at any sane pose (prisms ~0.002 rad, real warp >= 0.09).
        bool s12_active = false;
        if (cfg.s12_mode > 0) {
            double const fwd_s = double(hull->get_fwd());
            double const xT_s = (fwd_s > 0) ? hull->get_min(false).x : hull->get_max(false).x;
            if (s12_guard_ < 0) {
                double const Lwl_s = std::max(1e-3, hull->get_length_wl());
                double const bA_s = hull->get_deadrise_at_x(xT_s + fwd_s * 0.15 * Lwl_s);
                double const bF_s = hull->get_deadrise_at_x(xT_s + fwd_s * 0.60 * Lwl_s);
                s12_guard_ = (bA_s > 0.0 && bF_s > 0.0 &&
                              std::abs(bF_s - bA_s) > S12_WARP_SPREAD_FLOOR) ? 1 : 0;
            }
            if (s12_guard_ == 1) {
                s12_active = true;
                // Station fixed point: the station sits frac*lambda*b forward of the transom and
                // the load-coefficient lambda depends on (tau_eff, beta_eff). With the true
                // planform (default) the station comes from the tau-independent geometric mean
                // wetted length, so one pass is exact; without it, two passes converge the
                // weakly-coupled inversion.
                double lam_s = cfg.use_true_planform ? lambda_used : lambda;
                for (int it = 0; it < 2; it++) {
                    double const x_s = xT_s + fwd_s * cfg.s12_station_frac * std::max(lam_s, 0.2) * b;
                    double const bl_s = hull->get_deadrise_about(x_s, 0.25 * b);
                    double const beta_eff = (bl_s > 0.0)
                        ? std::max(1e-3, std::min(bl_s, radians(50.0))) : beta;
                    double const delta = hull->get_quarter_buttock_delta(x_s, 0.25 * b);
                    tau = std::max(radians(0.05), std::min(tau_raw + delta, radians(30.0)));
                    tau_deg = degrees(tau);
                    if (cfg.s12_mode == 1) { beta = beta_eff; }
                    // Re-invert the load-coefficient lambda at the effective attitude: it places
                    // the Wagner CoP and, without the true planform, spans the pressure planform.
                    double const C_L0_s = savitsky_clo(degrees(beta), C_Lbeta);
                    lambda = std::min(cfg.lambda_cap, std::max(0.2,
                        savitsky_lambda(std::pow(tau_deg, 1.1), Cv2, C_L0_s, 1e-7, 1e-6)));
                    if (cfg.use_true_planform) { break; }   // geometric station: one pass exact
                    lam_s = lambda;
                    beta_scale = lambda;
                    lambda_used = lambda;
                }
                sin_tau = std::sin(tau);
                cos_tau = std::cos(tau);
            }
        }

        if (out_lambda) { *out_lambda = lambda_used; }   // the exact mean planform λ the field used
        if (out_tau) { *out_tau = tau; }                 // the exact (clamped) trim the kernel used

        // Mean stagnation-line angle — used only for the wetted-planform geometry
        // (lambdaY below). The pressure is integrated against the LOCAL deadrise, so the
        // per-cell stagnation/wave-rise angles are recomputed inside the loop.
        double const tan_tau = std::tan(tau);
        double const tan_beta_mean = std::tan(beta);
        double const alpha_mean = std::atan2(M_PI * tan_tau, 2.0 * tan_beta_mean);
        double const tan_alpha = std::max(1e-3, std::tan(alpha_mean));

        // Reference point in the hull frame for moments: longitudinal distance of CG
        // from the transom, in beam units (positive forward).
        double const fwd = double(hull->get_fwd());
        double const x_cg = hull->get_cg(true).x;
        double const x_transom = (fwd > 0) ? hull->get_min(false).x : hull->get_max(false).x;
        double const LCG_b = (x_cg - x_transom) * fwd / b;

        // Local deadrise profile vs distance forward of the transom (beam units),
        // sampled once and interpolated in the integration loop. For monohedral hulls
        // it is flat (= mean); for warped hulls it rises forward, which is what shifts
        // the centre of pressure aft and lowers the running trim. Degenerate sections
        // fall back to the mean deadrise.
        int const MB = 11;
        std::vector<double> beta_prof(MB);
        for (int m = 0; m < MB; m++) {
            // s12 mode 1 (as-published equivalent prism): the constant effective deadrise
            // replaces the local profile — beta already holds beta_eff here.
            if (s12_active && cfg.s12_mode == 1) { beta_prof[m] = beta; continue; }
            double const d = (double(m) + 0.5) / double(MB) * std::max(beta_scale, 1e-3);
            double const bl = hull->get_deadrise_at_x(x_transom + fwd * d * b);
            beta_prof[m] = (bl > 0.0) ? std::min(bl, radians(50.0)) : beta;
        }
        auto beta_at = [&](double dist_fwd) -> double {
            double u = dist_fwd / std::max(beta_scale, 1e-3);
            u = std::min(1.0, std::max(0.0, u));
            double const t = u * double(MB) - 0.5;
            int i0 = std::min(MB - 1, std::max(0, int(std::floor(t))));
            int const i1 = std::min(MB - 1, i0 + 1);
            double f = std::min(1.0, std::max(0.0, t - double(i0)));
            return beta_prof[i0] * (1.0 - f) + beta_prof[i1] * f;
        };

        // Y grid is symmetric over [-0.5, 0.5]; integrate one side (Y in [0, 0.5])
        // and double, matching the localflow factor-of-2 convention.
        double const dY = 0.5 / double(cfg.ny);
        double const dX_unit = 1.0 / double(cfg.nx);  // X step is per-column, scaled by λ_Y

        double drag = 0.0, lift = 0.0, torque = 0.0;
        // Dynamic vs hydrostatic split of the planing lift+moment (Stage-B CoP diagnostic).
        double lift_dyn = 0.0, torque_dyn = 0.0, lift_stat = 0.0, torque_stat = 0.0;

        #pragma omp parallel for reduction(+:drag, lift, torque, lift_dyn, torque_dyn, lift_stat, torque_stat)
        for (int j = 0; j < cfg.ny; j++) {

            double const Y = (j + 0.5) * dY;  // midpoint, 0 < Y < 0.5
            double const Y_safe = std::min(Y, 0.499);

            // Wetted length to the transom at this transverse station. The wetted
            // planform is the triangle bounded by the straight stagnation (spray-root)
            // line and the transom, as the method describes: the keel (Y=0) is wetted
            // furthest forward and the chine (Y=0.5) least, with keel/chine asymmetry
            // λk - λc = 0.5/tanα = tanβ/(π·tanτ) (Savitsky's spray-root geometry) and the
            // mean wetted length preserved at λ (the apex straddles the mean at Y=0.25).
            // The forward boundary comes to a keel point; there is no fitted constant here.
            double const lambdaY = cfg.use_true_planform ? true_lamY[j]
                                                     : (lambda + (0.25 - Y) / tan_alpha);
            if (lambdaY <= 1e-3) {
                continue;  // chine runs dry where the triangle has tapered to zero
            }

            double const dX = lambdaY * dX_unit;
            double const dA = sq(b) * dX * dY;

            // Y-only stagnation factor.
            double const Y14 = std::pow(Y_safe, 1.4);
            double const tail = (0.5 - Y_safe) / std::max(0.51 - Y_safe, 1e-6);
            double const Py_stag = std::max(0.0, (1.02 - 0.25 * Y14) * tail);

            for (int i = 0; i < cfg.nx; i++) {

                double const X = (i + 0.5) * dX;

                // Local deadrise at this cell's longitudinal station ((λY - X) beam units
                // forward of the transom). Warped hulls carry high deadrise forward, so
                // the forward cells generate less lift; the centre of pressure — and hence
                // the running trim — come out lower than a single-deadrise model gives.
                double const dist_fwd = std::max(0.0, lambdaY - X);
                double const beta_l = beta_at(dist_fwd);
                double const tan_beta_l = std::max(1e-3, std::tan(beta_l));
                double const alpha_l = std::atan2(M_PI * tan_tau, 2.0 * tan_beta_l);
                double const alpha_W_l = std::atan2(tan_tau, tan_beta_l);
                double const tan_alpha_l = std::max(1e-3, std::tan(alpha_l));
                double const tan_alpha_W_l = std::max(1e-3, std::tan(alpha_W_l));
                double const sin_alpha2_l = sq(std::sin(alpha_l));

                // Transverse factor and longitudinal scaling at the local deadrise.
                double const Py = std::max(0.0, (1.02 - 0.05 * (degrees(beta_l) + 5.0) * Y14) * tail);
                double const Pmax_q_safe = std::max(Py_stag * sin_alpha2_l, 1e-6);
                double const C = Py * 0.006 * std::pow(tau_deg, 1.3);
                double const K = std::pow(std::max(C, 1e-12) / 2.588, 1.5) / std::pow(Pmax_q_safe, 1.5);

                // Longitudinal pressure distribution PL/q = C X^(1/3)/(X+K).
                double const Xc = std::pow(std::max(X, 1e-12), 1.0/3.0);
                double const PL_q = C * Xc / std::max(X + K, 1e-9);

                // Transom roll-off PT.
                double const rem = std::max(0.0, lambdaY - X);
                double const rem14 = std::pow(std::max(rem, 1e-12), 1.4);
                double const PT = rem14 / (rem14 + 0.05);

                double const Pdyn_q = PL_q * PT;

                // Hydrostatic contribution. Note: PSTATIC formula combines X and Y
                // in beam units; clip negative values that occur near the transom.
                double const Pstat_q_full = std::max(0.0,
                    2.0 * PT * Py * sin_tau * (X + Y * (1.0/tan_alpha_l - 1.0/tan_alpha_W_l)) / Cv2);
                // When use_hydrostatic is OFF the Pstat LIFT+MOMENT is handed to the hull buoyancy
                // (it double-counts the Archimedes CB), but its DRAG is KEPT: Archimedes buoyancy is
                // purely vertical, so the horizontal (drag) component of the bottom hydrostatic
                // pressure is not double-counted and must stay in the resistance.
                double const Pstat_q = cfg.use_hydrostatic ? Pstat_q_full : 0.0;

                double const P_q = Pdyn_q + Pstat_q;          // lift & moment pressure
                double const P = q * P_q;

                if (out_pressure_field) {
                    std::size_t idx = std::size_t(j) * std::size_t(cfg.nx) + std::size_t(i);
                    (*out_pressure_field)[idx] = P;
                }

                // Force components: pressure acts normal to the bottom (tilted by τ).
                // Vertical lift component: P·cos(τ); horizontal drag component: P·sin(τ).
                double const dF = P * dA;
                lift += dF * cos_tau;
                // Drag keeps the full hydrostatic bottom pressure (see note above).
                double const dF_drag = q * (Pdyn_q + Pstat_q_full) * dA;
                drag += dF_drag * sin_tau;

                // Moment about the lateral axis through CG.
                // Longitudinal lever (positive forward) from CG: x_cell_forward_of_transom - LCG.
                // x_cell_forward_of_transom in beam units is (λY - X) (Morabito §10.8).
                double const lever_x = ((lambdaY - X) - LCG_b) * b * cos_tau;
                // Vertical lever from CG: cell is at the hull bottom, height ~ X·sin(τ)
                // above the transom keel; for the cos-pressed lift-driven moment this is
                // small compared to the horizontal lever and we follow the localflow
                // convention of using planform lever only.
                torque += dF * cos_tau * lever_x;
                // Split the same lift/moment into dynamic (Pdyn) and hydrostatic (Pstat) parts.
                double const dF_dyn = q * Pdyn_q * dA, dF_stat = q * Pstat_q * dA;
                lift_dyn += dF_dyn * cos_tau; lift_stat += dF_stat * cos_tau;
                torque_dyn += dF_dyn * cos_tau * lever_x; torque_stat += dF_stat * cos_tau * lever_x;
            }
        }

        // Optionally re-place the DYNAMIC resultant at the Savitsky/Wagner CoP (frac*lambda fwd of
        // transom) instead of the integral's centroid, so the planing CoP marches forward with Cv.
        // The hydrostatic (Pstat) moment keeps its integral position; with Pstat off it is zero and
        // the hull buoyancy carries the hydrostatic restoring moment.
        if (cfg.dynamic_cop_fraction >= 0.0 || cfg.cop_savitsky) {
            double frac = cfg.dynamic_cop_fraction;
            if (cfg.cop_savitsky) {
                // Savitsky 1964 Eq31 speed-dependent CoP fraction (Cp/lambda) at the operating point,
                // replacing the constant 0.75 high-Cv asymptote. Parameter-free.
                double const ls = std::max(lambda, 1e-3);
                frac = std::min(1.0, std::max(0.3, 0.75 - 1.0 / (5.21 * sq(Cv / ls) + 2.39)));
            } else if (cfg.cop_cv_slope != 0.0) {
                // Warp gate: fore-aft deadrise spread sampled at FIXED Lwl-fractions (a stable
                // hull-shape property, not the running wetted length). The hinge (cop_warp_floor)
                // keeps prismatic AND mildly-warped hulls at the literature 3/4 CoP, so the ramp
                // activates only on strongly-warped forms (Begovic warp2/warp3).
                double const Lwl_g = std::max(1e-3, hull->get_length_wl());
                double const bA = hull->get_deadrise_at_x(x_transom + fwd * 0.15 * Lwl_g);
                double const bF = hull->get_deadrise_at_x(x_transom + fwd * 0.60 * Lwl_g);
                double warp = (bF > 0.0 && bA > 0.0) ? (bF - bA) : 0.0;
                warp = std::max(0.0, warp - cfg.cop_warp_floor);
                frac = std::min(1.0, std::max(0.3,
                    frac + cfg.cop_cv_slope * warp * (cfg.cop_cv_ref - Cv)));
            }
            double const lp = frac * lambda;                       // beam units fwd of transom
            double const lever_lp = (lp - LCG_b) * b * cos_tau;
            // Hook 2 (warped-transition lift closure): in the transition window on warped hulls, blend
            // the dynamic moment off the fixed Wagner CoP toward the enlarged planform's own integral
            // centroid (torque_dyn), driven by the SEPARATE warp_moment_amp so the bow-up moment is
            // dialled independently of the Hook-1 lift. wl_blend = 0 when moment OFF ⇒ byte-identical
            // to the Wagner placement (the Hook-1 lift then acts at the unchanged ¾-length CoP).
            double const torque_dyn_wagner = lift_dyn * lever_lp;
            double const wl_blend = std::min(1.0, cfg.warp_moment_amp * warp_lift_gate);
            torque = torque_stat + (1.0 - wl_blend) * torque_dyn_wagner + wl_blend * torque_dyn;
        }

        // Account for both sides (the loop covered half the hull, port or starboard).
        drag *= 2.0;
        lift *= 2.0;
        torque *= 2.0;

        if (out_split) {
            *out_split = glm::dvec4(lift_dyn * 2.0, torque_dyn * 2.0,
                                    lift_stat * 2.0, torque_stat * 2.0);
        }

        return {drag, lift, torque};
    }

    /// Mean stagnation/spray-root line angle α used to set the wetted-planform shape
    /// (per-strip wetted length λ_Y = λ − Y/tanα). tanα = π·tanτ/(2·tanβ): the π/2 is the
    /// Wagner spray-root pile-up, so the keel/chine wetted-length asymmetry is 2/π of the
    /// raw calm-water value. Uses the MEAN deadrise (the pressure itself uses local β).
    /// Exposed for validation of the planform geometry. Returns radians.
    double get_stagnation_angle(double speed) {
        if (speed <= 1e-2) {
            return 0.0;
        }
        double const beta = std::max(1e-3, std::min(hull->get_deadrise(), radians(50.0)));
        double const tau = std::max(radians(0.05), std::min(-hull->get_pitch(), radians(30.0)));
        return std::atan2(M_PI * std::tan(tau), 2.0 * std::tan(beta));
    }

    /// The Savitsky-inverted mean wetted length-to-beam ratio λ actually used by the
    /// bottom-pressure integral at the current attitude (the load-coefficient λ, NOT the
    /// geometric keel/chine value). Exposed for validation against measured wetted length
    /// (e.g. Fridsma ℓ_m/b). Mirrors the λ derivation in get_drag_lift_torque exactly.
    double get_lambda(double speed) {
        if (speed <= 1e-2) {
            return 0.0;
        }
        double const g = env->get_gravity();
        double const V = speed;
        double const b = std::max(hull->get_beam_chine(), 1e-3);
        double const Cv2 = std::max(sq(V / std::sqrt(g * b)), 1e-12);
        double const beta = std::max(1e-3, std::min(hull->get_deadrise(), radians(50.0)));
        double const tau = std::max(radians(0.05), std::min(-hull->get_pitch(), radians(30.0)));
        double const tau_deg = degrees(tau);
        double const vol = hull->get_reference_displacement();
        double const C_Lbeta = g * vol / (0.5 * sq(V) * sq(b));
        double const beta_deg_l = degrees(beta);
        double const C_L0 = savitsky_clo(beta_deg_l, C_Lbeta);
        double const lambda = std::min(cfg.lambda_cap, std::max(0.2,
            savitsky_lambda(std::pow(tau_deg, 1.1), Cv2, C_L0, 1e-7, 1e-6)));
        if (!cfg.use_true_planform) {
            return lambda;
        }
        double f_pile = cfg.pileup;
        if (cfg.auto_pileup) {
            double const fwd = double(hull->get_fwd());
            double const x_transom = (fwd > 0) ? hull->get_min(false).x : hull->get_max(false).x;
            double const ba = degrees(hull->get_deadrise_at_x(x_transom + 0.3 * b * fwd));
            double const b_aft = (ba > 0.0) ? ba : degrees(beta);
            f_pile = std::min(3.0, std::max(0.5, 0.1624 + 0.0589 * b_aft + 0.2221 * tau_deg));
        }
        TrueWettedProfile const tw = get_true_wetted_profile(hull, cfg.ny, b, f_pile);
        return tw.valid() ? std::max(0.2, tw.mean) : 0.0;
    }

private:

    /// Mean wetted length-to-beam ratio from keel and chine intersections with the
    /// current waterplane. λ = (Lk + Lc) / (2 b). Capped to Morabito's validity range.
    double compute_lambda(double b) {

        std::vector<double> xs, zs;

        // Keel buttock (y = 0)
        hull->slice_buttock(0.0, xs, zs, true);
        double Lk = polyline_x_range(xs);

        // Chine buttock (y at half-beam)
        xs.clear(); zs.clear();
        hull->slice_buttock(0.5 * b, xs, zs, true);
        double Lc = polyline_x_range(xs);

        double const mean = 0.5 * (Lk + Lc) / b;
        return std::max(0.2, std::min(mean, 5.5));
    }

    static double polyline_x_range(const std::vector<double>& xs) {
        if (xs.empty()) {
            return 0.0;
        }
        double xmin = xs.front(), xmax = xs.front();
        for (double x : xs) {
            if (x < xmin) xmin = x;
            if (x > xmax) xmax = x;
        }
        return std::abs(xmax - xmin);
    }
};


#endif // PLANING_H
