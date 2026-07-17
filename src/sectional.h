#ifndef SECTIONAL_H
#define SECTIONAL_H

#include <vector>
#include <cmath>
#include <algorithm>
#include "hull.h"

/// 2.5D sectional water-entry planing kernel (the warp closure). At the current pose and a
/// given speed, integrates a smooth added-mass momentum distribution over the wetted length to
/// give the planing LIFT and its CENTRE OF PRESSURE from ONE distribution -- the construction
/// that resolves the Morabito lift/moment decoupling (its lift is forward-loaded, its moment
/// Wagner-pinned). The aft-loaded added-mass force uses the LOCAL effective angle of attack
/// (wetted-bottom buttock slope), which == the keel trim on a prism (do-no-harm) and stays steep
/// on a warped stern, so the lift does not collapse with the keel trim.
///
/// This C++ kernel MIRRORS the validated Python prototype auto_validation/proto_2dt_solver.py
/// (smooth kappa(beta) Wagner added mass + monotone-q transom telescoping + displaced-volume LCB).
/// At the Begovic EFD attitude, with one Begovic-mono-pinned scale K (0.947 momentum-only; 0.482
/// with the default-on crossflow term, calibrate_K.py; the momentum-only pin generalises
/// out-of-sample to Series 62 at 3.9% lift -- proto_2dt_oos.py): lift MAPE mono 1.4 / warp1 1.8 /
/// warp2 7.9 / warp3 4.5 %, CoP-LCG RMS 0.025-0.036 L. warp3 closes BOTH lift and moment.
class Sectional {

    std::shared_ptr<Hull> hull;
    std::shared_ptr<Environment> env;
    // ---- Tunable levers (grouped; mutated only via the set_* forwarders below) ----
    struct Config {
        int nx = 480;   // x-stations for the wetted-band sampling. 480 (was 240) resolves the spray-root
                        // forward limit finely enough that it no longer steps between stations as the pose
                        // varies at high Fn / low trim -> removes the warp3 attitude saw-teeth (the moment(trim)
                        // fold). Only the warp-gated sectional hulls pay the cost (mono/warp1/prisms never call it).
        double smooth_frac = 0.10;   // moving-average window / wetted length
        // Begovic-mono-pinned added-mass scale (NOT re-fit on the predicted hulls). Pinned WITH the
        // crossflow term at its published coefficient (calibrate_K.py 1.0 -> 0.482); the momentum-only
        // pin (crossflow off) is 0.947. The lift splits ~half added-mass / half crossflow, consistent
        // with the Shuford/Savitsky linear-plus-crossflow decomposition.
        double K = 0.482;
        // Upper cap on the local effective AoA (wetted-bottom buttock slope). The aftmost station sits on
        // the near-vertical transom face, where slice_buttock reads an unphysical ~85deg slope; the old 20deg
        // cap still let that inflate the transom momentum flux -> a spurious aft force spike that placed the
        // CoP too far aft (the warp3 over-rotation / sinkage regression). A planing buttock slope is at most
        // a few x the keel trim, so a tighter cap removes the transom-edge artifact. Re-pin K when changing it.
        double tloc_cap = 20.0;      // degrees
        // Transom-edge slope fix (default OFF, byte-identical when off). The aftmost wetted station(s) sit on
        // the near-vertical transom face, where slice_buttock reads an unphysical ~85deg buttock slope that
        // inflates the transom momentum flux into a spurious AFT force spike -> the warp3 CoP lands ~0.08L too
        // far aft (the over-rotation / sinkage regression). When on, any station slope above tfix_thr is
        // replaced by the nearest forward valid (planing-surface) slope -- removing the transom-face spike
        // while PRESERVING the warp stern's legitimate buttock amplification. Pose-fixed do-no-harm on prisms.
        bool transom_fix = false;
        double tfix_thr = 12.0;      // degrees; slopes above this are the transom-face artifact
        // Faltinsen (2005) 2.5D transom correction. The flow ventilates at the transom: over the aft
        // ~0.5B of the wetted bottom the pressure is atmospheric (clean separation), so the added-mass
        // momentum rate there is NOT a sustained planing force. Zeroing it removes the spurious
        // near-transom load that the pure momentum-rate kernel puts at the worst lever (just fwd of the
        // transom), which biases the warped-hull centre of pressure aft -> under-trim. It FIXES warp3
        // trim+sinkage (winpc A/B: trim 1.40->3.22deg vs EFD 2.90, sinkage |err| 6.7->5.2mm) but, by
        // correcting the over-rise, sheds the spurious wetted-area friction that was cancelling the
        // model's pre-existing high-Fn warped-drag deficit, so it RAISES warp3 R MAPE (12.6->16.1%).
        // Default OFF to preserve the calibration-free resistance headline; enable via set_*_ventilation
        // / SECTIONAL_VENT. Gated upstream by the warp spread (prisms run Morabito; OOS/mono untouched).
        // See memory transom-ventilation-attitude-fix.
        bool transom_ventilation = false;
        double vent_beam_frac = 0.5;  // ventilated length / beam (Faltinsen 0.5B)
        // Optional nonlinear forward transition-lift term (default OFF). The best pose-fixed prototype
        // form is a spray-root stagnation load scaled by the local effective AoA and the forward-opening
        // wetted width gradient. It is kept separate and default-off so the current kernel stays byte-identical.
        bool forward_term = false;
        double forward_cs = 0.0;
        // Crossflow lift term (colleague round-3 #2, DEFAULT ON at the published value since the
        // diagcal campaign 2026-07-03). The pure momentum force
        // f_a = -U dq/dx is aft-loaded (loads where the wetted width grows), so it under-loads the
        // steep forward buttocks and floors the warped CoP -> warp3 under-trim (R2). The classical
        // remedy (Payne planing formulation; Shuford NACA TR-1355 crossflow lift) is a second
        // sectional lift term f_c = coef * 0.5 rho C_Dc(beta) * 2c(x) * (U tan theta(x))^2 cos beta.
        // Because it scales with tan^2(theta) it forward-loads exactly the steep sections the
        // momentum term ignores -> moves the sectional CoP forward -> bow-up moment that raises the
        // warp3 trim. Warp-gated by construction (the sectional closure runs only on warp2/warp3),
        // so mono/warp1/prisms are byte-identical. C_Dc = 4/3 is Shuford's (NACA TR-1355) measured
        // crossflow drag coefficient, found INDEPENDENT of deadrise -- so it is a constant, applied
        // as published (crossflow_coef=1.0 = the published value; the coef knob is retained only for
        // the sensitivity sweep; 0 restores the legacy momentum-only kernel, pair with K=0.947).
        // K is re-pinned on mono WITH this term ON, so no new free parameter is introduced.
        // Adoption provenance (fixes_ab/diagcal/DIAG.md): winpc Begovic 12.3 -> 10.8%, warp3
        // planing |dtrim| 1.61 -> 1.11 deg; the Phase-2 coef fine grid shows warp3 trim saturating
        // at ~1.0 while warp2's R-preference for higher coef is attitude error-compensation (proxy),
        // so the published value is also the evidence-backed one. Honest cost: warp2 sinkage
        // 5.6 -> 8.2 mm (price of the out-of-sample mono pin; see the Kw2 diagnostic in DIAG.md).
        double crossflow_coef = 1.0;
        // Added-mass kernel selector. false (DEFAULT) = the validated prototype kappa(beta) Wagner kernel
        // (proto_2dt_solver.py): width c = kappa(beta)*d/tan(beta), added mass rho*(pi/2)*c^2, monotone-q
        // transom telescoping. true = the Payne Cm(beta) similarity variant — DEPRECATED: it under-lifts
        // in the equilibrium solve (warp3 sinkage collapses negative) and must NOT be the default.
        bool similarity_model = false;
    } cfg;
    // STATIC hull properties captured at construction (the hull is floated then). get_beam_chine()
    // and get_length_wl() are POSE-DEPENDENT (they re-detect the chine / waterline at the running
    // attitude, ~5-10% low), so the chine cap and the x-sampling envelope must use the static values;
    // the waterline plane z_waterline is fixed, captured for d = wl - z_keel.
    double b_half0 = 0.0, Lwl0 = 0.0, wl0 = 0.0;

    // ---- small array helpers (mirror numpy where the Python prototype used it) ----
    static std::vector<double> movavg(std::vector<double> const& y, int win) {
        int const n = int(y.size());
        if (n < 3) { return y; }                          // reflect padding needs n>=3; nothing to smooth
        win = std::max(3, win | 1);                       // force odd
        win = std::min(win, 2 * n - 1);                   // clamp so the reflect pad (win/2) stays < n (OOB for smooth_frac>=2)
        int const pad = win / 2;
        std::vector<double> yp; yp.reserve(n + 2 * pad);
        for (int k = pad; k >= 1; k--) yp.push_back(y[k]);            // reflect: y[pad..1]
        for (int k = 0; k < n; k++) yp.push_back(y[k]);
        for (int k = 2; k <= pad + 1; k++) yp.push_back(y[n - k]);    // reflect: y[-2..-pad-1]
        std::vector<double> out(n, 0.0);
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int k = 0; k < win; k++) s += yp[i + k];
            out[i] = s / win;
        }
        return out;
    }

    // numpy.gradient on a UNIFORM grid (interior central, one-sided edges).
    static std::vector<double> gradient_uniform(std::vector<double> const& f, double dx) {
        int const n = int(f.size());
        std::vector<double> g(n, 0.0);
        if (n == 1) return g;
        g[0] = (f[1] - f[0]) / dx;
        g[n - 1] = (f[n - 1] - f[n - 2]) / dx;
        for (int i = 1; i < n - 1; i++) g[i] = (f[i + 1] - f[i - 1]) / (2.0 * dx);
        return g;
    }

    // numpy.gradient on a NON-UNIFORM grid (for buttock slopes).
    static std::vector<double> gradient_nonuniform(std::vector<double> const& f,
                                                   std::vector<double> const& x) {
        int const n = int(f.size());
        std::vector<double> g(n, 0.0);
        if (n < 2) return g;
        g[0] = (f[1] - f[0]) / (x[1] - x[0]);
        g[n - 1] = (f[n - 1] - f[n - 2]) / (x[n - 1] - x[n - 2]);
        for (int i = 1; i < n - 1; i++) {
            double const h1 = x[i] - x[i - 1], h2 = x[i + 1] - x[i];
            g[i] = (h1 * h1 * f[i + 1] + (h2 * h2 - h1 * h1) * f[i] - h2 * h2 * f[i - 1])
                   / (h1 * h2 * (h1 + h2));
        }
        return g;
    }

    static double trapz_uniform(std::vector<double> const& y, double dx) {
        int const n = int(y.size());
        if (n < 2) return 0.0;
        double s = 0.5 * (y[0] + y[n - 1]);
        for (int i = 1; i < n - 1; i++) s += y[i];
        return s * dx;
    }

    static double interp_lin(double x, std::vector<double> const& xp, std::vector<double> const& fp) {
        int const n = int(xp.size());
        if (n == 0 || x < xp.front() - 1e-9 || x > xp.back() + 1e-9) return std::nan("");
        if (x <= xp.front()) return fp.front();
        if (x >= xp.back()) return fp.back();
        int lo = 0, hi = n - 1;
        while (hi - lo > 1) { int m = (lo + hi) / 2; if (xp[m] <= x) lo = m; else hi = m; }
        double const t = (x - xp[lo]) / (xp[hi] - xp[lo]);
        return fp[lo] + t * (fp[hi] - fp[lo]);
    }

    static double kappa_deadrise(double beta_rad) {
        double const c = std::cos(beta_rad);
        return 1.0 + (M_PI / 2.0 - 1.0) * c * c;          // pi/2 at beta=0 -> 1 at 90deg
    }

    // Per-station wetted-mean |buttock slope| (the local effective AoA), == keel on a prism.
    struct ButtockField {
        std::vector<std::vector<double>> px, pslope;      // one profile per sampled buttock
    };
    ButtockField buttock_field(double b_half) {
        ButtockField bf;
        int const ny = 6;
        std::vector<double> bx, bz;
        for (int j = 0; j < ny; j++) {
            double const frac = 0.02 + (0.5 - 0.02) * double(j) / double(ny - 1);
            double const y = frac * b_half;
            hull->slice_buttock(y, bx, bz, true);
            if (bx.size() < 3) continue;
            // sort by x
            std::vector<int> idx(bx.size());
            for (int k = 0; k < int(bx.size()); k++) idx[k] = k;
            std::sort(idx.begin(), idx.end(), [&](int a, int b){ return bx[a] < bx[b]; });
            std::vector<double> sx(bx.size()), sz(bx.size());
            for (int k = 0; k < int(idx.size()); k++) { sx[k] = bx[idx[k]]; sz[k] = bz[idx[k]]; }
            std::vector<double> sl = gradient_nonuniform(sz, sx);
            std::vector<double> ang(sl.size());
            for (int k = 0; k < int(sl.size()); k++) ang[k] = std::atan(std::abs(sl[k]));
            bf.px.push_back(sx); bf.pslope.push_back(ang);
        }
        return bf;
    }
    double buttock_slope_at(ButtockField const& bf, double x) const {
        double sum = 0.0; int cnt = 0;
        for (std::size_t p = 0; p < bf.px.size(); p++) {
            double const v = interp_lin(x, bf.px[p], bf.pslope[p]);
            if (std::isfinite(v)) { sum += v; cnt++; }
        }
        return cnt ? sum / cnt : std::nan("");
    }

public:
    Sectional(std::shared_ptr<Hull> ship_hull, std::shared_ptr<Environment> environment)
        : hull(ship_hull), env(environment) { capture_static(); }

    /// Re-capture the static chine beam / waterline length / waterline plane. Call after the hull is
    /// floated if it was not yet at construction (the constructor captures automatically).
    void capture_static() {
        b_half0 = 0.5 * hull->get_beam_chine();
        Lwl0 = hull->get_length_wl();
        wl0 = hull->get_waterline();
    }

    void set_nx(int n) { cfg.nx = std::max(20, n); }
    void set_smooth_frac(double f) { cfg.smooth_frac = std::max(0.0, f); }
    void set_scale(double k) { cfg.K = k; }
    double get_scale() const { return cfg.K; }
    void set_tloc_cap(double deg) { cfg.tloc_cap = std::max(2.0, deg); }
    double get_tloc_cap() const { return cfg.tloc_cap; }
    void set_transom_fix(bool on) { cfg.transom_fix = on; }
    bool get_transom_fix() const { return cfg.transom_fix; }
    void set_forward_term(bool on, double cs) {
        cfg.forward_term = on;
        if (cs >= 0.0) cfg.forward_cs = cs;
    }
    bool get_forward_term() const { return cfg.forward_term; }
    double get_forward_cs() const { return cfg.forward_cs; }
    void set_similarity_model(bool on) { cfg.similarity_model = on; }
    bool get_similarity_model() const { return cfg.similarity_model; }
    /// Crossflow lift coefficient scale (colleague round-3 #2). 0 (default) = OFF, byte-identical.
    /// Folds the overall C_Dc scale; the deadrise shape 1.33 cos^2(beta) is built in.
    void set_crossflow(double coef) { cfg.crossflow_coef = std::max(0.0, coef); }
    double get_crossflow() const { return cfg.crossflow_coef; }
    void set_transom_ventilation(bool on) { cfg.transom_ventilation = on; }
    bool get_transom_ventilation() const { return cfg.transom_ventilation; }
    void set_vent_beam_frac(double f) { cfg.vent_beam_frac = std::max(0.0, f); }
    // Replace transom-face artifact slopes (> tfix_thr) by the nearest forward valid planing slope.
    void apply_transom_fix(std::vector<double>& tloc) const {
        if (!cfg.transom_fix) return;
        int const n = int(tloc.size());
        double const thr = radians(cfg.tfix_thr);
        for (int k = 0; k < n; k++) {
            if (tloc[k] <= thr) continue;
            int j = k + 1; while (j < n && tloc[j] > thr) j++;
            if (j < n) tloc[k] = tloc[j];
            else { int b = k - 1; while (b >= 0 && tloc[b] > thr) b--; if (b >= 0) tloc[k] = tloc[b]; }
        }
    }

    /// Raw sectional results at the current pose and speed U (NO scale applied to La/Lf):
    ///   x = {La, cop_a, lcb, L, x_t, Lf, cop_f}. Lift = K*La + Cs*Lf + buoyancy; net CoP from
    ///   (K*La,cop_a)+(Cs*Lf,cop_f)+(buoy,lcb).
    /// Returns nwet=0 (La=0, cop_a/lcb/cop_f=NaN) if fewer than 5 wetted stations; check nwet, not La.
    struct Result { double La, cop_a, lcb, L, x_t, Lf, cop_f, Lc, cop_c, Dc, Da; int nwet; };
    Result added_mass(double U) const {
        Result R{0.0, std::nan(""), std::nan(""), 0.0, 0.0, 0.0, std::nan(""), 0.0, std::nan(""), 0.0, 0.0, 0};
        double const rho = env->get_density();
        double const b_half = b_half0;                        // static chine half-beam (pose-stable)
        double const wl = wl0;                                // fixed waterline plane
        double const tau_keel = -hull->get_pitch();           // set_pitch(-trim): keel trim = -pitch
        double const Lwl = Lwl0;                              // static hull length (x-sampling envelope)

        // --- sample the wetted extent via keel submergence d(x) (robust where max(ys) is not) ---
        std::vector<double> xs_all(cfg.nx), d_all(cfg.nx), beta_all(cfg.nx);
        std::vector<double> ys, zs;
        for (int i = 0; i < cfg.nx; i++) {
            double const x = 0.0005 + (Lwl - 0.0005) * double(i) / double(cfg.nx - 1);
            xs_all[i] = x;
            hull->slice_section(x, ys, zs, true);
            double dd = 0.0;
            if (!zs.empty()) {
                double zmin = zs[0];
                for (double z : zs) zmin = std::min(zmin, z);
                dd = wl - zmin;
            }
            d_all[i] = dd;
            double const be = hull->get_deadrise_at_x(x);
            beta_all[i] = (be == be) ? be : 0.0;
        }
        // wetted band: the LONGEST CONTIGUOUS run of wetted stations. The raw [first,last] span would
        // include any interior dry gap (d<=thr) and inject spurious zeros into c_w/q; a planing hull's
        // wetted region is contiguous (verified on all Begovic hulls/poses), so the longest run IS it,
        // and this hardens the kernel for any hull where a gap could appear.
        int i0 = -1, i1 = -1, best_len = 0;
        for (int i = 0; i < cfg.nx; ) {
            if (d_all[i] <= 1e-4) { i++; continue; }
            int j = i; while (j < cfg.nx && d_all[j] > 1e-4) j++;     // [i, j) is one wetted run
            if (j - i > best_len) { best_len = j - i; i0 = i; i1 = j - 1; }
            i = j;
        }
        if (i0 < 0 || best_len < 5) return R;

        int const n = i1 - i0 + 1;
        double const dx = xs_all[1] - xs_all[0];
        std::vector<double> xs(n), d(n), beta(n), tloc(n);
        ButtockField const bf = buttock_field_const(b_half);
        for (int k = 0; k < n; k++) {
            xs[k] = xs_all[i0 + k]; d[k] = d_all[i0 + k]; beta[k] = beta_all[i0 + k];
            double const s = buttock_slope_at(bf, xs[k]);
            tloc[k] = std::isfinite(s) ? s : tau_keel;
        }
        apply_transom_fix(tloc);
        double const x_t = xs.front(), x_s = xs.back(), L = x_s - x_t;
        int const win = std::max(3, int(cfg.smooth_frac * n));

        // --- smooth geometry: Wagner wetted half-width, chine-capped at the constant chine half-beam ---
        std::vector<double> c_w(n), A_sub(n), tant(n);
        for (int k = 0; k < n; k++) {
            double const bc = std::min(std::max(beta[k], radians(2.0)), radians(80.0));
            double const tb = std::tan(bc);
            if (cfg.similarity_model) {
                double const Cc = M_PI / 2.0;
                double const c_dyn = d[k] * Cc / tb;
                c_w[k] = std::max(std::min(c_dyn, b_half), 0.0);
            } else {
                double const cg = d[k] / tb;
                double cw = std::min(kappa_deadrise(bc) * cg, b_half);
                c_w[k] = std::max(cw, 0.0);
            }
            double const d_chine = b_half * tb;
            A_sub[k] = std::max(0.0, (d[k] <= d_chine) ? d[k] * d[k] / tb
                                                       : 2.0 * b_half * d[k] - b_half * b_half * tb);
            double const tl = std::min(std::max(tloc[k], radians(0.1)), radians(cfg.tloc_cap));
            tant[k] = std::tan(tl);
        }
        c_w = movavg(c_w, win);
        tant = movavg(tant, win);

        // --- displaced-volume centroid (LCB) ---
        double vol = trapz_uniform(A_sub, dx);
        std::vector<double> A_sub_x(n);
        for (int k = 0; k < n; k++) A_sub_x[k] = A_sub[k] * xs[k];
        R.lcb = (vol > 1e-9) ? trapz_uniform(A_sub_x, dx) / vol : std::nan("");

        // --- added-mass momentum (aft-loaded; the warp-selective term) ---
        std::vector<double> q(n);
        for (int k = 0; k < n; k++) {
            double m_a;
            if (cfg.similarity_model) {
                double const bc = std::min(std::max(beta[k], radians(2.0)), radians(80.0));
                double const Cm = (M_PI / 2.0) * std::pow(1.0 - bc / (M_PI / 2.0), 2.0);
                m_a = rho * Cm * c_w[k] * c_w[k];
            } else {
                m_a = rho * (M_PI / 2.0) * c_w[k] * c_w[k];
            }
            q[k] = m_a * (U * tant[k]);
        }
        q = movavg(q, win);
        // section monotonically wets going aft (added mass only grows) until transom separation;
        // enforce q non-increasing forward so the entry rate f_a stays >= 0.
        for (int k = n - 2; k >= 0; k--) q[k] = std::max(q[k], q[k + 1]);
        std::vector<double> dq = gradient_uniform(q, dx);
        std::vector<double> f_a(n), f_a_x(n);
        // Transom-ventilation correction (see transom_ventilation note): zero the momentum-rate force
        // over the aft 0.5B and forbid suction ahead of it. vent-OFF reproduces the prior kernel exactly.
        double const vent_x = x_t + cfg.vent_beam_frac * 2.0 * b_half;
        for (int k = 0; k < n; k++) {
            double fz = -U * dq[k];
            if (cfg.transom_ventilation) fz = (xs[k] < vent_x) ? 0.0 : std::max(fz, 0.0);
            f_a[k] = fz; f_a_x[k] = fz * xs[k];
        }
        double const La = trapz_uniform(f_a, dx);
        R.La = La;
        R.cop_a = (std::abs(La) > 1e-9) ? trapz_uniform(f_a_x, dx) / La : std::nan("");
        // f9drag DIAGNOSTIC: water-entry (induced) drag of the momentum lift = the
        // horizontal projection of f_a on the local bottom inclined at buttock slope
        // tloc, Da = integral f_a*tan(tloc) dx. Warp-weighted (tloc >> keel trim aft).
        {
            std::vector<double> f_a_d(n);
            for (int k = 0; k < n; k++) f_a_d[k] = f_a[k] * tant[k];
            R.Da = trapz_uniform(f_a_d, dx);
        }
        if (cfg.forward_term && cfg.forward_cs > 0.0) {
            std::vector<double> dcw = gradient_uniform(c_w, dx);
            std::vector<double> f_f(n), f_f_x(n);
            for (int k = 0; k < n; k++) {
                double const grow = std::max(-dcw[k], 0.0);
                double const vloc = U * tant[k];
                f_f[k] = rho * vloc * vloc * grow;   // RAW Lf; the consumer (lift()) applies forward_cs
                f_f_x[k] = f_f[k] * xs[k];
            }
            R.Lf = trapz_uniform(f_f, dx);
            R.cop_f = (R.Lf > 1e-9) ? trapz_uniform(f_f_x, dx) / R.Lf : std::nan("");
        }
        // Crossflow lift (Payne/Shuford): forward-loaded second sectional term ~ c(x)*(U tanθ)^2.
        // Uses the SAME smoothed/capped c_w and tant the momentum term uses (tant = tan of the
        // capped buttock slope, so the transom-face artifact is excluded). Default OFF -> Lc=0.
        if (cfg.crossflow_coef > 0.0) {
            // Shuford (NACA TR-1355) crossflow drag coefficient: C_D,c = 4/3, determined from tests
            // and (Shuford's explicit finding) INDEPENDENT of deadrise ("as long as the angle of
            // deadrise was constant for the entire beam, C_D,c did not vary with the angle of
            // deadrise"). So it is a constant, NOT reduced by cos^2(beta). The crossflow drag opposes
            // the section's vertical entry, so the strip force is vertical (like the momentum f_a) --
            // no cos(beta) projection. Deadrise enters only through the wetted width c_w and slope.
            double const C_Dc = 4.0 / 3.0;
            std::vector<double> f_c(n), f_c_x(n), f_c_d(n);
            for (int k = 0; k < n; k++) {
                double const vloc = U * tant[k];              // local normal velocity (capped buttock slope)
                f_c[k] = cfg.crossflow_coef * 0.5 * rho * C_Dc * (2.0 * c_w[k]) * vloc * vloc;
                f_c_x[k] = f_c[k] * xs[k];
                // f9drag DIAGNOSTIC: the crossflow force is taken vertical (Lc); its dropped
                // horizontal (drag) projection on the local bottom inclined at buttock slope
                // tloc is f_c*tan(tloc). Warp-keyed (tant large on warped aft sections), ~U^2.
                f_c_d[k] = f_c[k] * tant[k];
            }
            R.Lc = trapz_uniform(f_c, dx);
            R.cop_c = (R.Lc > 1e-9) ? trapz_uniform(f_c_x, dx) / R.Lc : std::nan("");
            R.Dc = trapz_uniform(f_c_d, dx);
        }
        R.L = L; R.x_t = x_t; R.nwet = n;
        return R;
    }

    /// Debug: per wetted-station arrays at the current pose+speed (xs, d, beta_deg, tloc_deg, c_w, q, f_a).
    void debug_profiles(double U, std::vector<double>& o_xs, std::vector<double>& o_d,
                        std::vector<double>& o_beta, std::vector<double>& o_tloc,
                        std::vector<double>& o_cw, std::vector<double>& o_q,
                        std::vector<double>& o_fa) const {
        o_xs.clear(); o_d.clear(); o_beta.clear(); o_tloc.clear(); o_cw.clear(); o_q.clear(); o_fa.clear();
        double const rho = env->get_density();
        double const b_half = b_half0;
        double const wl = wl0;
        double const tau_keel = -hull->get_pitch();
        double const Lwl = Lwl0;
        std::vector<double> xs_all(cfg.nx), d_all(cfg.nx), beta_all(cfg.nx);
        std::vector<double> ys, zs;
        for (int i = 0; i < cfg.nx; i++) {
            double const x = 0.0005 + (Lwl - 0.0005) * double(i) / double(cfg.nx - 1);
            xs_all[i] = x;
            hull->slice_section(x, ys, zs, true);
            double dd = 0.0;
            if (!zs.empty()) { double zmin = zs[0]; for (double z : zs) zmin = std::min(zmin, z); dd = wl - zmin; }
            d_all[i] = dd;
            double const be = hull->get_deadrise_at_x(x);
            beta_all[i] = (be == be) ? be : 0.0;
        }
        int i0 = -1, i1 = -1, best_len = 0;       // longest contiguous wetted run (see added_mass)
        for (int i = 0; i < cfg.nx; ) {
            if (d_all[i] <= 1e-4) { i++; continue; }
            int j = i; while (j < cfg.nx && d_all[j] > 1e-4) j++;
            if (j - i > best_len) { best_len = j - i; i0 = i; i1 = j - 1; }
            i = j;
        }
        if (i0 < 0 || best_len < 5) return;
        int const n = i1 - i0 + 1;
        double const dx = xs_all[1] - xs_all[0];
        std::vector<double> xs(n), d(n), beta(n), tloc(n);
        ButtockField const bf = buttock_field_const(b_half);
        for (int k = 0; k < n; k++) {
            xs[k] = xs_all[i0 + k]; d[k] = d_all[i0 + k]; beta[k] = beta_all[i0 + k];
            double const s = buttock_slope_at(bf, xs[k]);
            tloc[k] = std::isfinite(s) ? s : tau_keel;
        }
        apply_transom_fix(tloc);
        int const win = std::max(3, int(cfg.smooth_frac * n));
        std::vector<double> c_w(n), tant(n);
        for (int k = 0; k < n; k++) {
            double const bc = std::min(std::max(beta[k], radians(2.0)), radians(80.0));
            double const tb = std::tan(bc);
            if (cfg.similarity_model) {
                double const Cc = M_PI / 2.0;
                c_w[k] = std::max(std::min(d[k] * Cc / tb, b_half), 0.0);
            } else {
                c_w[k] = std::max(std::min(kappa_deadrise(bc) * d[k] / tb, b_half), 0.0);
            }
            tant[k] = std::tan(std::min(std::max(tloc[k], radians(0.1)), radians(cfg.tloc_cap)));
        }
        c_w = movavg(c_w, win); tant = movavg(tant, win);
        std::vector<double> q(n);
        for (int k = 0; k < n; k++) {
            if (cfg.similarity_model) {
                double const bc = std::min(std::max(beta[k], radians(2.0)), radians(80.0));
                double const Cm = (M_PI / 2.0) * std::pow(1.0 - bc / (M_PI / 2.0), 2.0);
                q[k] = rho * Cm * c_w[k] * c_w[k] * (U * tant[k]);
            } else {
                q[k] = rho * (M_PI / 2.0) * c_w[k] * c_w[k] * (U * tant[k]);
            }
        }
        // NB: this diagnostic INTENTIONALLY diverges from the production added_mass() kernel, which
        // smooths q (q = movavg(q, win)) and then enforces monotone-q non-increasing forward. Both are
        // omitted here on purpose so o_q / o_fa expose the raw per-station bow-impact force. Do not treat
        // o_fa as the exact force added_mass() integrates.
        // q = movavg(q, win);
        // for (int k = n - 2; k >= 0; k--) q[k] = std::max(q[k], q[k + 1]);
        std::vector<double> dq = gradient_uniform(q, dx);
        // Transom-ventilation correction, identical to added_mass(): gated on the toggle, with the
        // ventilated zone measured as vent_beam_frac off the transom station x_t. (Previously this was
        // hard-coded ON at 0.5B with no x_t offset, so the diagnostic disagreed with the real kernel.)
        // vent-OFF (default) leaves the raw momentum-rate force, matching added_mass with vent off.
        double const vent_x = xs.front() + cfg.vent_beam_frac * 2.0 * b_half;
        for (int k = 0; k < n; k++) {
            o_xs.push_back(xs[k]); o_d.push_back(d[k]); o_beta.push_back(degrees(beta[k]));
            o_tloc.push_back(degrees(tloc[k])); o_cw.push_back(c_w[k]); o_q.push_back(q[k]);
            double f_z = -U * dq[k];
            if (cfg.transom_ventilation) f_z = (xs[k] < vent_x) ? 0.0 : std::max(f_z, 0.0);
            o_fa.push_back(f_z);
        }
    }

    /// Total planing lift [N] = K*La + Cs*Lf + buoyancy(=displaced weight), at the current pose and speed.
    double lift(double U) const {
        Result const r = added_mass(U);
        if (r.nwet < 5) return std::nan("");
        double const buoy = hull->get_displaced_mass() * env->get_gravity();
        return cfg.K * r.La + (cfg.forward_term ? cfg.forward_cs * r.Lf : 0.0)
               + (cfg.crossflow_coef > 0.0 && std::isfinite(r.Lc) ? r.Lc : 0.0) + buoy;
    }

private:
    // non-const wrapper needed because slice_buttock mutates the hull's scratch; keep added_mass const
    ButtockField buttock_field_const(double b_half) const {
        return const_cast<Sectional*>(this)->buttock_field(b_half);
    }
};

#endif
