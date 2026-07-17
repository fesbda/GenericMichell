#ifndef RESISTANCE_H
#define RESISTANCE_H

#include <cmath>
#include <stdexcept>
#include "thinship.h"
#include "localflow.h"
#include "holtrop.h"
#include "planing.h"
#include "sectional.h"

enum class ResistanceMethod {
    Auto,        // Fn_vol-blended dispatch: Unified ↔ Morabito/Savitsky
    Unified,     // Michell (wave) + Noblesse (local pressure) + Holtrop (viscous)
    Integral,    // Michell pressure-integral form
    Holtrop,     // Empirical Holtrop-Mennen
    Morabito,    // Morabito (2014) bottom pressures + Savitsky-style wave term
    Savitsky     // Savitsky (1964) closed-form total resistance
};

enum class ResistanceGridPrecision {
    Coarse,
    Fine,
    Finest
};

/// Sampled propeller-disc wake field, in non-dimensional V/V_ship units.
/// Polar grid: r[i] in [R_inner, R_outer], theta[j] in [0, 2*pi), measured
/// from the +z axis (top of disc) and increasing clockwise looking forward
/// (ITTC convention). The disc plane is the y-z plane through center.x.
/// Arrays are row-major (n_radial * n_azimuthal), theta varies fastest.
struct PropellerWakeField {
    glm::dvec3 center;
    double R_inner = 0.0;
    double R_outer = 0.0;
    int n_radial = 0;
    int n_azimuthal = 0;
    double ship_speed = 0.0;

    std::vector<double> r;            // size n_radial
    std::vector<double> theta;        // size n_azimuthal
    std::vector<double> axial;        // (1 - u_x*fwd/V) — freestream is 1
    std::vector<double> radial;       // outward positive
    std::vector<double> tangential;   // CW (positive when fluid swirls with right-hand
                                      //  rule about the +x_fwd axis seen from astern)

    // Area-weighted statistics over the annulus dA = r dr dtheta:
    double mean_axial = 1.0;          // V_A / V (= 1 - w_T)
    double wake_fraction = 0.0;       // Taylor wake fraction, 1 - mean_axial
    double axial_std = 0.0;
    double axial_min = 1.0;
    double axial_max = 1.0;
    double non_uniformity = 0.0;      // (axial_max - axial_min) / mean_axial
    double radial_rms = 0.0;
    double swirl_rms = 0.0;
    double volumetric_flow = 0.0;     // m^3/s, dimensional
    double disc_area = 0.0;           // m^2
};

/// Main class to predict the resistance of a vessel
class Resistance {

private:

    std::shared_ptr<Hull> original_hull;
    std::shared_ptr<Hull> hull; // a copy to freely transform
    std::shared_ptr<Environment> env;
    ThinShip michell;
    Holtrop holtrop;
    LocalFlow localflow;
    Morabito morabito;
    Savitsky savitsky;
    Sectional sectional;
    ResistanceGridPrecision mesh_precision = ResistanceGridPrecision::Finest;
    static constexpr bool do_fixes = true;
    // Frozen full-ship resistance for the drag-couple pitch moment (R*a), set by the
    // find_equilibrium outer fixed-point; 0 => the couple term in residual_core is inert.
    double drag_couple_R_ = 0.0;

    // 2.5D sectional water-entry attitude closure (DEFAULT ON, warp-gated). When ON,
    // find_equilibrium's planing LIFT and MOMENT come from the Sectional added-mass kernel (lift =
    // K*La, moment = lift*(cop_a - cg)) instead of the Morabito dynamic lift + Wagner-3/4 CoP. This
    // resolves the lift/moment decoupling that over-trims warped hulls. The Morabito DRAG path is
    // unchanged (resistance still uses it at the converged attitude); only the attitude solve moves.
    // WARP-GATED: the out-of-sample test showed the added-mass CoP is a WORSE prism CoP than the
    // calibrated Wagner-3/4 (it harms Fridsma/Taunton/Series62), so it is applied only where Wagner-3/4
    // fails — on warped hulls, gated by the fore-aft deadrise spread (same gate as cop_cv_ramp). The
    // floor excludes prismatic AND mildly-warped forms; 0 = ungated (the universal version).
    // ---- Tunable levers (grouped; mutated only via the set_* forwarders below) ----
    struct Config {
        bool use_sectional_attitude = true;
        bool sectional_delta_mode = false;
        bool attitude_wave_blend = false;
        // Immersion-scaled Michell wave lift on the Morabito planing path (wavelift campaign item 4,
        // F5, DEFAULT OFF). Batch-3's causal near-field keeps physical Michell wave lift in the
        // vertical balance through planing; on the mono/warp1 Morabito path (no sectional
        // replacement) it over-lifts -> the +10-14 mm planing over-rise. Scale wave_y/wave_z by the
        // immersed-volume fraction V(pose)/V(design): the wave-making, hence its lift, fades as the
        // hull rises -- tied to immersion, not the blanket (1-w_plan) fade that zeroed it and
        // collapsed Series 62 trim (physfix C6). ->1 in the displacement regime so that suite is
        // untouched by construction. Independent of attitude_wave_blend (that stays the C6 path).
        bool wave_lift_volume_scale = false;
        // Outer trim solve: when true, the equilibrium trim is the ITP (Interpolate-Truncate-
        // Project) root of the heave-balanced moment inside the steepest-restoring crossing
        // bracket, instead of the local linear-regression estimate. ITP is unconditionally
        // convergent (never worse than bisection) and superlinear, so it drives |moment| to a
        // true root (|M|/WL ~1e-6 vs the regression's ~1e-3) while keeping the same global scan
        // for branch selection. Default ON: a local Begović sweep (bare defaults) showed it is
        // do-no-harm-or-better — overall MAPE -0.31pp, warp3 -1.17pp (the multi-crossing fold hull),
        // mono/warp1/warp2 within +/-0.13pp. Set false to recover the linear-regression root.
        bool use_itp_root = true;
        double sectional_deficit_lo = 0.30;
        double sectional_deficit_hi = 0.80;
        double sectional_fn_lo = -1.0;
        double sectional_fn_hi = -1.0;
        double sectional_tloc_cap = 20.0;
        double sectional_warp_floor = 0.13;   // rad; sectional active only where (beta_fwd - beta_aft) > floor
        // F8 pressure-drag consistency (DEFAULT OFF). When the sectional kernel carries the
        // planing lift (replacement mode, warp-gated, K-scaled), the resistance still integrates
        // Morabito's bottom pressure computed from Morabito's own (smaller) lift. A planing
        // pressure field lifting L has the horizontal component L*tan(tau), so the drag the
        // assembly misses is (sec_lift - morabito_lift)*tan(tau_run), planing-weighted. Added
        // WITH the vp component so the residuary shortfall gate and the NSS total envelope
        // re-gate around it (regression supply is displaced, not stacked). Warp-gated by the
        // same spread floor as the sectional attitude path -> prisms/OOS byte-identical.
        bool sectional_pressure_drag = false;

        // Smooth Fn_vol-based blend between displacement-regime (Unified) and
        // planing-regime (Morabito/Savitsky) pressure forces. Cubic smoothstep
        // over [Fn_min, Fn_max]; exact 0 below Fn_min, exact 1 above Fn_max.
        // After the get_deadrise_at_x fix made the warp-resolved Morabito pressure correct, a
        // global sweep showed [0.6, 1.0] Pareto-dominated the old [0.7, 1.4] (full-set MAPE
        // 10.1->8.3%). A later blend-window sweep (2026-07, MATRIX=blend) then showed extending
        // the UPPER bound keeps the Michell displacement wave engaged through the semi-planing
        // band and lowers the Begovic warped-series MAPE monotonically (9.8->8.7% at [0.6,2.0]);
        // trim/sinkage stay ~unchanged, so the gain is resistance-side. Trade: it costs the
        // Series62 low-deadrise planing hull ~+1.3pt (global lever, not warp-gated). Adopted
        // [0.6, 2.0] as the balanced point.
        double planing_blend_fn_min = 0.6;
        double planing_blend_fn_max = 2.0;
        bool blend_use_fnb = false;  // false = blend on Fn_vol (default); true = on FnB (chine-beam)

        // Separate (decoupled) blend band used ONLY for the running-attitude solve
        // (find_equilibrium: the lift partition and the trim blend), NOT for the
        // resistance-component blends. The resistance band [planing_blend_fn_min,max] is
        // calibrated to the Begovic total and must stay; but driving the ATTITUDE off that
        // same band declares buoyancy ~0 at Fn_vol>=1 (the hump) and lifts the hull onto
        // plane far too early -- the SINKAGE_ISSUE.md premature-planing error (model rises
        // ~+90mm where the hull squats ~-18mm). A slower attitude band lets buoyancy hand
        // over gradually (squat through the hump, rise only at high speed) while leaving the
        // validated resistance partition untouched. Default now = the resistance band [0.6,2.0]
        // (2026-07 blend adoption). The old [0.6,3.0] attitude investigation (2026-06-25) had
        // reintroduced the SINKAGE_ISSUE premature-planing error, but at [0.6,2.0] on the current
        // attitude physics (true-planform lift + hydrostatic double-count removal + Wagner CoP) the
        // MATRIX=blend sweep confirms the Begovic trim/sinkage MAE stay ~baseline (|dsink|
        // 10.2->10.0mm), so the wider band is safe here.
        double planing_blend_attitude_fn_min = 0.6;
        double planing_blend_attitude_fn_max = 2.0;

        // F7 effective buoyancy (DEFAULT OFF). At planing attitude the hull does not carry
        // its full calm-plane displaced volume: the transom runs dry and the aft free surface
        // is depressed, so crediting g*displaced_mass over-counts the hydrostatic lift.
        // Savitsky's own buoyant term implies ~HALF the geometric pose buoyancy at full
        // planing (the classic dry-transom half-buoyancy result: 0.19 vs 0.37 W on the
        // Begovic mono), and the imposed-pose lift residual pins the removal at exactly
        // a = 0.500 with a chine-beam-Froude ramp (fixes_ab/f7buoy). Both the static lift
        // and the static torque scale by (1 - frac * w_plan * s(CvB)): uniform scaling of
        // the hydrostatic field, centroid kept at CB. w_plan keeps the displacement regime
        // byte-identical; the CvB smoothstep [cv_lo, cv_hi] phases the removal with transom
        // ventilation (dry by CvB ~ 3.5). NOTE: the sectional K was calibrated WITH the
        // over-count in the balance -- enabling this requires the K re-pin (~1.15 on the
        // warp2 anchor vs 0.482 legacy).
        double eff_buoy_frac = 0.0;    // fraction removed at full ramp; 0 = off, 0.5 = physics
        double eff_buoy_cv_lo = 0.1;   // CvB = U/sqrt(g*b_chine) smoothstep lower edge
        double eff_buoy_cv_hi = 3.5;   // upper edge (removal saturates here)
        // Moment treatment of the removal. true = scale static_torque too (uniform field
        // scaling: removal centroid at CB -- the pose CB sits 0.39-0.50 m from the transom,
        // ~0.25 m AFT of where the imposed-pose moment residuals put the removal, ~ the CG).
        // false = lift-only (static_torque untouched: removal acts at the CG, which matches
        // the measured x_rem = cg + Tres/Lres band-for-band within ~0.1 m). A/B lever.
        bool eff_buoy_torque = true;

        // kshape Phase 1 (DEFAULT OFF): CG-acting sectional lift relief at speed. The F7
        // equilibrium K-ladder (fixes_ab/kshape) shows the flat sectional K over-floats the
        // warped hulls at frv >= 3.3 (+9-10 mm, K_ds0 -> 0.35-0.7) while the trim wants K
        // HELD at ~1.2 (d dt/dK ~ -1.8 deg/K): the DOFs pull K opposite ways, so the excess
        // top-band lift must be shed moment-neutrally (at the CG), not via the K scale
        // (whose lift and moment are locked by the kernel CoP). Physically: 3D flow relief
        // the 2.5D kernel misses at high speed ratio. Relief = frac * s(CvB) * sec_lift,
        // lift only (plan_torque untouched), Replacement mode only -> same warp gate as the
        // sectional attitude path, prisms/OOS byte-identical.
        double sec_relief_frac = 0.0;   // fraction of sec_lift shed at full ramp; 0 = off
        double sec_relief_cv_lo = 2.0;  // CvB smoothstep lower edge (frv ~2.4: float onset)
        double sec_relief_cv_hi = 3.6;  // upper edge (top band, frv ~4.2)

        // kshape Phase 1 (DEFAULT OFF): weight the F8 pressure-drag charge by the dynamic
        // share of the load. At hump CvB the sectional supply substitutes BUOYANT load
        // (the equilibrium ladder shows the hump lift is right at flat K: |ds| <= 1.5 mm),
        // and hydrostatic lift carries no L*tan(tau) pressure-drag cost -- charging the
        // full (sec - morabito) difference over-prices the hump (warp3 2.0-2.6 eR +17).
        // Weight = the same CvB smoothstep partition F7 uses ([eff_buoy_cv_lo, cv_hi]).
        bool sec_pdrag_dynshare = false;

        // pdyn (2026-07-10 pressure-term campaign, DEFAULT OFF; fixes_ab/pdyn/TARGETS.md):
        // replace the planing viscous-pressure channel (the Morabito x-projection) by the
        // CFD-derived effective-angle law
        //     D_p = L_dyn * tan(tau_run + delta0 + k_warp * dbutt).
        // The seven pressure/shear-decomposed STAR sims (targets.csv) show the dynamic
        // bottom pressure acting 1.3-2.2 deg STEEPER than the keel trim on prismatic hulls
        // (near-constant over Fn_vol 2.3-5.7) and up to +3.45 deg on the most-warped
        // Begovic hull, while the Morabito integral projects at/below the keel trim.
        // The dynamic lift is PARAMETER-FREE: at equilibrium total lift = W, and the
        // hydrostatic part is the posed buoyancy B, so the dynamic (pressure) lift is
        // L_dyn = W - B_pose exactly. The RANS decomposition confirms this to the newton
        // (GPPH top: W-B = 996-140 = 856 N vs computed dynamic lift 856); it replaced a
        // fitted buoyancy-share 'a' that was both less accurate and less physical
        // (2026-07-12 winpc: a=0 improved GPPH 5.0->4.5, Taunton 8.3->7.7, S62 8.9->8.6).
        bool pdyn_on = false;
        double pdyn_delta0 = 1.6;   // effective-angle offset [deg] (CFD prisms: 1.3-2.2)
        double pdyn_kwarp = 0.0;    // deg of delta per deg of quarter-beam buttock delta
        // Whisker-spray CvB cap (pdyn campaign, DEFAULT OFF c=0): the whisker term's growth
        // outruns the CFD shear room at the very top of the speed range (GPPH 10.03 -> 12.24
        // m/s: whisker 18.9 -> 47.9 N while CFD shear grows as v^2 and leaves ~10 N of room;
        // solver friction alone is 185 vs CFD S 194.7 there). Scale the whisker drag by
        // (1 - c*smoothstep(CvB over [lo, hi])); at 10-11 m/s (CvB ~ 4) the whisker is REAL
        // shear per the CFD (removing it costs -7%), so the cap engages above that.
        double wsk_cap_c = 0.0;
        double wsk_cap_cv_lo = 4.0;
        double wsk_cap_cv_hi = 5.0;
        // Whisker-spray warp-only scoping (see get_whisker_spray_drag): apply the term only to
        // warped hulls (RANS: spray is ~3% on a clean prism, over-modelled by the correlation).
        // DEFAULT OFF (byte-identical); ON removes the need for the CvB cap above.
        bool whisker_warp_only = false;
        // Transom head-release contribution to the propeller-disc axial inflow (default OFF).
        // amp = 1 is the unscaled Bernoulli prediction; len_scale stretches Doctors' hollow
        // length. Affects get_propeller_wake_field only — no resistance path reads these.
        bool transom_hr_on = false;
        double transom_hr_amp = 1.0;
        double transom_hr_len = 1.0;
        // Handover window: below cv_lo the Morabito projection stands (the hump channel is
        // already adequate there); the law takes over through [cv_lo, cv_hi] in CvB. The
        // GPPH deficit lives at CvB >~ 4 (Fn_vol >~ 4.5); the F7 window ([0.1, 3.5]) engaged
        // at the hump and over-charged it (+10..17% at Fn_vol 2.6-3.1 in the first fit).
        double pdyn_cv_lo = 2.5;
        double pdyn_cv_hi = 4.5;

        // Pre-planing residuary "hump" correction (Lever 2). Thin-ship Michell
        // under-predicts the residuary peak through the semi-planing hump (Fn_vol ~1-1.8):
        // verified that ~11 N is missing at Fn_vol 1.3 even at the measured attitude, with
        // neither the Michell wave nor the Morabito planing pressure carrying it. This adds
        // a compact residuary bump in R/Delta vs Fn_vol, gated to planing-form hulls and the
        // hump band so it cannot touch the displacement validation or the Fn_vol>2 region.
        // STOPGAP: a transparent single-amplitude bump, NOT the full cited Mercier-Savitsky
        // (1973) pre-planing regression; to be replaced by that regression once sourced.
        // Begovic sweep (sweep_preplaning_hump.py): the best constant setting is
        // coef ~= 0.03 over [0.6, 1.5] (peak Fn_vol ~1.05), which moves the pooled planing
        // MAPE 11.3% -> ~10.4% (every hull improves or holds), but it CANNOT reach MaxSurf's
        // 9.0% and, crucially, a fixed coef fitted to the data is calibration -- so it is
        // OFF by default to keep the method calibration-free. The deficit is hull-specific
        // (peaks at Fn_vol ~1.0, broader for w1/w2, sharp for w3), which a geometry-scaled
        // Mercier-Savitsky regression would capture properly.
        double preplaning_hump_coef = 0.0;   // peak R/Delta of the bump; 0 = off (default)
        double preplaning_fn_lo = 0.6;
        double preplaning_fn_hi = 1.5;

        // Fraction of the planing trim taken from Savitsky's predictive relation versus the
        // warp-resolved bottom-pressure moment balance: at full planing the trim is
        // (1 - f) * tau_moment + f * tau_Savitsky. 0.0 = the pure Morabito moment-equilibrium
        // trim (self-consistent: attitude set by the same bottom-pressure field that gives the
        // forces); 1.0 = pure Savitsky predictive trim.
        // Physics-first default 0.0 (2026-06-25): the pure Morabito moment-equilibrium trim. Under the
        // true-planform lift + the hydrostatic-moment double-count removal + the Wagner 3/4-length CoP, the
        // moment-equilibrium trim is now CORRECT and generalises (it no longer over-predicts the warped
        // trim, which was the 2026-06-24 reason for the 1.0 restore). See ATTITUDE_PHYSICS_FINDINGS.md.
        // 1.0 = the legacy Savitsky-predictive-trim / risen-pose method (set via set_planing_legacy()).
        double savitsky_trim_weight = 0.0;
        // Warp-gated Savitsky trim weight (warped-hull sinkage closure). When warp_w > 0 the trim weight
        // is ramped 0 -> warp_w by the fore-aft deadrise spread (smoothstep over [warp_lo, warp_hi] rad),
        // so prismatic/hard-chine hulls (spread ~ 0) stay at 0 and only warped hulls get the partial blend.
        // 0 (default) = use the flat savitsky_trim_weight above. Ported from the Python harness gate.
        double savitsky_trim_warp_w = 0.0;
        double savitsky_trim_warp_lo = 0.05;   // warp-spread smoothstep lower edge (rad)
        double savitsky_trim_warp_hi = 0.15;   // warp-spread smoothstep upper edge (rad)
        // Savitsky (2012) effective-attitude mode (0 = OFF; forwarded to morabito+savitsky via
        // set_savitsky2012). Kept here to gate the predict_trim beta_override in the blend:
        // with s12 active predict_trim self-samples the effective deadrise and converts the
        // effective trim back to keel pitch, so the wetted-mean override must not double-adjust.
        int savitsky2012_mode = 0;

        // Planing-friction wetted length: true = running (equilibrium) trim λ (default,
        // shorter); false = Savitsky's predictive-trim λ (the residuary-calibrated wetted
        // length that gives standalone Savitsky its Fridsma accuracy). Experiment lever for
        // the residuary closure: running-trim λ under-feeds the friction on prismatic hulls.
        bool planing_friction_running_trim = true;

        // Whisker-spray drag (Savitsky, DeLorme & Datla 2007, "Inclusion of Whisker Spray
        // Drag in Performance Prediction Method for High-Speed Planing Hulls", Marine
        // Technology 44(1):35-56). The thin spray sheet forward of the stagnation line carries
        // a viscous drag the bottom-pressure integral misses; it is largest at HIGH deadrise /
        // LOW trim — exactly the regime where Auto under-predicts (audit F3, Fridsma beta30,
        // Begovic w3). Cited, calibration-free, so ON by default; toggle for ablation.
        bool whisker_spray_on = true;

        // Whisker-spray friction line (colleague round-2 B0, DEFAULT OFF). The Savitsky-DeLorme-Datla
        // closed form takes the spray sheet laminar (Blasius 1.328/sqrt(Re)) below Re_ww = 1.5e6 and
        // a transitional Schoenherr line (0.074/Re^0.2 - 4800/Re) above. Spray sheets are thin,
        // disturbed, high-shear flows; the ITTC HSMV 1-2% laminar guideline is model-scale-specific.
        // When on, force the fully-turbulent line Cf = 0.074/Re^0.2 at ALL Re_ww (drop the
        // -4800/Re transitional subtraction, which goes negative at model-scale Re_ww). This is a
        // BOUNDING diagnostic for the spray channel's high-Fn headroom, not a tuned constant.
        bool whisker_spray_turbulent = false;

        // Whisker-spray deadrise input (physfix C2, DEFAULT OFF). The Savitsky-DeLorme-Datla 2007
        // closed form takes a single deadrise; get_deadrise() is the aft-half STATIC mean, which on
        // a warped hull misreports the deadrise in the wetted spray region. When on, feed the
        // wetted-half-beam-weighted local mean over the actual wetted length (get_wetted_mean_deadrise)
        // instead. Guarded so constant-deadrise hulls (fore-aft spread ~0, i.e. every prism series)
        // take the exact original scalar path -> byte-identical. Resistance-only term (not in the
        // force balance) so it cannot move the attitude or break cold==warm. Parameter-free.
        bool whisker_local_beta = false;

        // Whisker-spray SPRAY-ROOT deadrise (colleague round-2 B1, DEFAULT OFF). The stagnation
        // (spray-root) line sits at the FORWARD end of the wetted length; on a warped hull it
        // crosses the high-deadrise forward sections (25-35 deg), not the low aft bottom the mean
        // reports. Spray drag rises steeply with deadrise at low trim, so the section-local
        // deadrise at the spray root is the physically-correct input (colleague's Eq.9 point,
        // same philosophy as the chine fix). When on, feed the forward-station deadrise (the
        // WARP_SAMPLE_FWD_FRAC sample the warp gate already reads) as beta. Guarded on the
        // fore-aft spread so constant-deadrise hulls (every prism series incl. Taunton) take the
        // exact scalar path -> byte-identical. Takes precedence over whisker_local_beta.
        // DEFAULT ON since the wspray campaign (2026-07-03): with the crossflow attitude fix in
        // place, the Fn-gated spray-root supplies the warp-keyed super-v^2 high-Fn deficit
        // (diagcal F1b/F2) -- winpc Begovic 10.8 -> 10.1% (warp2 -0.8 / warp3 -1.5, mono
        // byte-identical). Fires on the fore-aft deadrise spread, so it ALSO engages Series 62's
        // genuine forebody warp (13.6 deg aft -> 19.9 deg fwd) and improves it 8.6 -> 8.1%; the
        // true constant-deadrise prisms (Fridsma, Taunton) and the displacement suite stay
        // byte-identical (spread ~0). fixes_ab/wspray/FINDINGS_stageA.md.
        bool whisker_spray_root = true;

        // Fn-gate on the spray-root beta migration (wspray campaign item 2, 2026-07-03). The
        // spray root wets the forward (finer) sections only once the hull planes and the bow-wave
        // steepens; below the gate the aft-mean beta is kept (byte-identical to the default
        // whisker). Un-gated, the beta_F swap taxes the over-predicted mid-Fn band (warp3 2.6-3.3
        // +6 pt MAPE) while the real warp-keyed super-v^2 deficit it targets is at Fn_vol >= 3.3
        // (diagcal F1b/F2). Smoothstep-blend beta from aft-mean to beta_F over Fn_vol in [lo, hi];
        // hi <= lo disables the gate (full beta_F at all speeds = the raw campaign-2 behaviour).
        // Window [3.2, 3.8] = the full-planing transition; robust (gate sensitivity 10.1-10.3% over
        // [2.8,3.4]..[3.0,4.0], fixes_ab/wspray/), reported as a declared window like the blend bound.
        double whisker_spray_root_fn_lo = 3.2;
        double whisker_spray_root_fn_hi = 3.8;

        // Warp-keyed turbulent spray Cf (hiband campaign, 2026-07-03; DEFAULT OFF). On a
        // strongly-warped hull the spray sheet sweeps across sections of changing deadrise, so
        // the conical clean-sheet assumption behind the laminar (model-scale) Cf line breaks;
        // take the turbulent line there instead. Keyed to the SAME fore-aft spread floor as the
        // sectional/crossflow closure (sectional_warp_floor, 0.13 rad) -- the hulls whose flow
        // needs the crossflow closure are the ones with the disturbed spray root -- and phased
        // by the SAME Fn_vol window as the spray-root beta migration. By construction: mono,
        // Fridsma/Taunton (spread ~0), warp1 (0.094) and Series 62 (0.108) keep the laminar
        // line byte-identical; only warp2 (0.180) / warp3 (0.168) engage. At the imposed EFD
        // attitude this supplies the remaining warp-FLAT ~5-6 N high-frv deficit that survives
        // the spray-root adoption (warp2 -12.4 -> -4.5%, warp3 -10.3 -> -2.9% in the 4.0-4.8
        // band; Begovic pooled imposed 9.19 -> 8.16%). fixes_ab/hiband/.
        // DEFAULT ON since 2026-07-04 (user adoption): winpc A/B overall 10.1 -> 9.8%,
        // warp2 11.7->10.4 / 10.4->7.9, mono/warp1 + all OOS byte-identical, cold==warm
        // exact. The warp3/2011 equilibrium cost (8.8 -> 11.3) is bounded by the w3rise
        // findings: half the warp3 top-band deficit is the 2011-vs-2012 EFD inter-campaign
        // spread (F6b) and the residual is warp3's attitude-side pure-U^2 lift deficit
        // (w3rise F1-F4), not this drag channel.
        bool whisker_spray_turb_warp = true;

        // Local free-surface wave envelope (DEFAULT ON, warp-selective). Deployable now that the sectional
        // attitude closure removes the warp over-trim (whose excess induced drag previously compensated this
        // wave deficit -- the two had to deploy together). The thin-ship Michell wave is
        // centerplane-linearized and under-resolves the WARPED-hull wave-making at the
        // semi-planing hump and at high Fn (warp2/3 under-predict -22/-24% even at the
        // imposed EFD attitude, exactly where the Michell wave is small). The LocalFlow
        // Noblesse near-field pressure drag tracks that deficit and is warp-selective
        // (D/R_loc ~ 0.84-0.99 at the hump/high-Fn; ~0 on monohedral). Two estimates of the
        // same wave-making, so the Auto wave term becomes max(Michell, R_noblesse) -- the
        // ENVELOPE, not a sum (avoids double-count where Michell already resolves the wave at
        // its frv~2.3 peak), gated by the planing weight so the displacement regime is
        // byte-identical. Parameter-free. See diag_localfs_wave.py / localfs-wave-envelope.
        bool localfs_wave_envelope = true;

        // Above-water air resistance (windage). Calm-air drag of the above-water hull moving
        // through still air, R_AA = 1/2 * rho_air * C_AA * A_VT * V^2 (ITTC 7.5-02-03-01.4 /
        // Holtrop-Mennen windage form), with A_VT the above-water transverse projected (frontal)
        // area from the hull silhouette and C_AA the windage drag coefficient. Grows as V^2, so
        // it is largest at high speed; on the model-scale towing-tank hulls it is only ~1% of R_T
        // at the top speeds, but it becomes a several-percent term at full scale. Method-blind
        // (added to every method's total). DEFAULT OFF so the validated totals are byte-identical;
        // enable with set_air_resistance(true). C_AA default 0.8 (bluff above-water hull, the term
        // already lumps the above-water skin friction); set_air_drag_coef to tune, or
        // Hull::set_frontal_area to override A_VT (e.g. to add a superstructure). NOTE: this is the
        // vessel's own calm-air windage at uniform relative wind = ship speed; it does NOT model a
        // height-varying atmospheric wind profile or the air boundary layer over the sea surface
        // (negligible here -- in the ship frame the still air and still water both translate at the
        // ship speed, so there is no air/water shear to drive such a layer).
        bool air_resistance_on = false;
        double air_drag_coef = 0.8;   // C_AA

        // Mercier & Savitsky (1973, Davidson Lab R-1667 / DTIC AD-764958) pre-planing
        // RESIDUARY resistance, to fill the wave/residuary deficit that thin-ship Michell
        // misses through the hump and at high deadrise (audit F3). Equation (6): a 14-term
        // Doust regression for R_T/Delta (total, ref 100,000-lb ship, C_A=0) over Fn_vol=1-2;
        // we subtract the reference Schoenherr friction to get the scale-independent residuary
        // and use it (max'd with the Michell wave, ramped over the band) as the Auto wave term.
        // Calibration-free (Mercier's published regression + geometric form params); ON by
        // default after validation (2026-06): adds only the residuary SHORTFALL beyond the
        // Michell wave + Morabito induced our model already carries, gated to planing forms in
        // Fn_vol[0.9,3.0] so the displacement regime is untouched. Net effect small and never
        // negative across the three series; a free, cited improvement that does not reach into
        // the displacement gate. NB: read its form descriptors at the design pose (see
        // get_mercier_savitsky_residuary) or warm-start continuation breaks cold==warm.

        bool mercier_savitsky_on = true;
        bool radojcic_nss_on = true;
        // Radojcic NSS TOTAL-resistance envelope (default ON). When ON, the Auto planing total is
        // floored to the Radojcic regression's total RT (planing-weighted): R += w_plan * max(0, R_rado - R).
        // This carries the high-Fn drag the component model loses when the residuary hump-closures taper to
        // zero (see decomp: measured residuary keeps growing while MS/Radojcic-residuary -> 0). Radojcic is a
        // hard-chine planing regression, so it does-no-harm on the prism series (it floors them to their own
        // regression, which they already match). Slender hulls (lp/V^1/3 > 8) get no Radojcic total -> inert.
        // DEFAULT ON: closes the high-Fn under-prediction offset (Begovic pooled bias -5.0->-4.0%, monohedral
        // -8.9->-5.1%; Taunton signed -30.1->-28.9%) with do-no-harm on the prism series (Fridsma 11.7->11.7,
        // Series62 8.0->7.8). Gated to high Fn_vol (residuary-closure death) and deadrise>13deg (Radojcic
        // reliability). Set false to recover the pre-envelope behaviour.
        bool radojcic_total_envelope = true;
        // Deadrise gate on the Radojcic RESIDUARY-excess term (default ON, adopted 2026-07-05,
        // fixes_ab/nssbeta). The NSS regression has no deadrise input, and the envelope above
        // already restricts itself to beta > ~13 deg for exactly that reason ("over-predicts
        // very-low-deadrise planing prisms"); the sibling residuary-excess term never got the
        // same gate. On the corrected Fridsma grid the ungated residuary charges ~2 N mid-Cv on
        // beta=10 light prisms whose measured total is friction+induced alone (beta10 bias +18%,
        // worst cell +32%). When ON, scale the residuary excess by the SAME smoothstep([11,15]
        // deg) on the static area-mean deadrise: Fridsma beta10 degates; every other validated
        // hull is untouched by construction (area-mean deadrise: Begović 16.7-18.0, Series 62
        // 16.5 despite its 13-deg transom, Taunton 22.5, Fridsma beta20/30).
        bool radojcic_residuary_beta_gate = true;

        // Drag-couple pitch moment (R*a): the moment of the horizontal resistance about the CG,
        // omitted by the vertical-force-only moment balance. Bow-down, grows with Fn; parameter-free
        // rigid-body mechanics (no fitted constant). Default OFF (prototype). warp_floor 0 = ungated
        // (universal); >0 = applied only where the fore-aft deadrise spread exceeds the floor (rad).
        bool drag_couple_on = false;
        double drag_couple_warp_floor = 0.0;
        double drag_couple_warp_hi = 0.15;   // upper edge of the warp-spread smoothstep (rad)
    } cfg;
    // pdyn warp key (deg), cached once per instance from the DESIGN hull: the posed
    // readout is mesh-noisy (28 deg observed on the posed GPPH) and per-eval gates
    // flicker at scan extremes (see the s12_guard_ note in planing.h).
    double pdyn_dbutt_deg_ = -1e9;



    // Table VI: rows = 14 terms (mult: 1,X,U,W,XZ,XU,XW,ZU,ZW,W2,XW2,ZX2,UW2,WU2),
    // cols = Fn_vol 1.0,1.1,...,2.0.
    static constexpr double MS_A[14][11] = {
        {0.06473,0.10776,0.09483,0.03475,0.03013,0.03163,0.03194,0.04343,0.05036,0.05612,0.05967},
        {-0.48630,-0.88787,-0.63720,0,0,0,0,0,0,0,0},
        {-0.01030,-0.01634,-0.01540,-0.00978,-0.00664,0,0,0,0,0,0},
        {-0.06490,-0.13444,-0.13580,-0.05097,-0.05540,-0.10543,-0.08599,-0.13289,-0.15597,-0.18661,-0.19758},
        {0,0,-0.16046,-0.21880,-0.19359,-0.20540,-0.19442,-0.18062,-0.17813,-0.18288,-0.20152},
        {0.10628,0.18186,0.16803,0.10434,0.09612,0.06007,0.06191,0.05487,0.05099,0.04744,0.04645},
        {0.97310,1.83080,1.55972,0.43510,0.51820,0.58230,0.52349,0.78195,0.92859,1.18569,1.30026},
        {-0.00272,-0.00389,-0.00309,-0.00198,-0.00215,-0.00372,-0.00360,-0.00332,-0.00308,-0.00244,-0.00212},
        {0.01089,0.01467,0.03481,0.04113,0.03901,0.04794,0.04436,0.04187,0.04111,0.04124,0.04343},
        {0,0,0,0,0,0.08317,0.07366,0.12147,0.14928,0.18090,0.19769},
        {-1.40962,-2.46696,-2.15556,-0.92663,-0.95276,-0.70895,-0.72057,-0.95929,-1.12178,-1.38644,-1.55127},
        {0.29136,0.47305,1.02992,1.06392,0.97757,1.19737,1.18119,1.01562,0.93144,0.78414,0.78282},
        {0.02971,0.05877,0.05198,0.02209,0.02413,0,0,0,0,0,0},
        {-0.00150,-0.00356,-0.00303,-0.00105,-0.00140,0,0,0,0,0,0},
    };

public:
    void set_radojcic_nss(bool on) { cfg.radojcic_nss_on = on; }
    /// F8 sectional pressure-drag consistency toggle (see cfg.sectional_pressure_drag).
    void set_sectional_pressure_drag(bool on) { cfg.sectional_pressure_drag = on; }

    /// kshape: CG-acting sectional lift relief at speed (see cfg.sec_relief_frac).
    void set_sectional_lift_relief(double frac, double cv_lo = 2.0, double cv_hi = 3.6) {
        cfg.sec_relief_frac = frac;
        cfg.sec_relief_cv_lo = cv_lo;
        cfg.sec_relief_cv_hi = std::max(cv_lo + 1e-3, cv_hi);
    }

    /// kshape: weight the F8 pressure-drag charge by the dynamic load share (see
    /// cfg.sec_pdrag_dynshare).
    void set_sectional_pdrag_dynshare(bool on) { cfg.sec_pdrag_dynshare = on; }

    /// pdyn effective-angle pressure projection (see cfg.pdyn_on; fixes_ab/pdyn).
    /// delta0 in degrees; [cv_lo, cv_hi] = the CvB handover window from the Morabito
    /// projection to the law. The dynamic lift is parameter-free (W - B_pose).
    void set_pdyn_projection(bool on, double delta0_deg = 1.6,
                             double k_warp = 0.0, double cv_lo = 2.5, double cv_hi = 4.5) {
        cfg.pdyn_on = on;
        cfg.pdyn_delta0 = delta0_deg;
        cfg.pdyn_kwarp = k_warp;
        cfg.pdyn_cv_lo = cv_lo;
        cfg.pdyn_cv_hi = std::max(cv_lo + 1e-3, cv_hi);
        pdyn_dbutt_deg_ = -1e9;
    }

    /// pdyn: the planing pressure drag at the CURRENT pose when cfg.pdyn_on (the caller
    /// applies the displacement/planing blend weight): a smoothstep handover in CvB from
    /// the Morabito x-projection (hump, adequate there) to the CFD-derived law
    /// D_law = L_dyn * tan(tau_run + delta), with the parameter-free dynamic lift
    /// L_dyn = W - B_pose and delta = delta0 + k_warp * (quarter-beam buttock delta).
    double get_pdyn_pressure_drag(double speed) {
        Hull& H = const_cast<Hull&>(*hull);
        double const g = env->get_gravity();
        double const morabito_x = morabito.get_drag_lift_torque(speed).x;
        double const b = std::max(1e-6, H.get_beam_chine());
        double const cvb = speed / std::sqrt(g * b);
        double r = (cvb - cfg.pdyn_cv_lo) / (cfg.pdyn_cv_hi - cfg.pdyn_cv_lo);
        r = (r <= 0.0) ? 0.0 : (r >= 1.0 ? 1.0 : r * r * (3.0 - 2.0 * r));
        if (r <= 0.0) {
            return morabito_x;
        }
        double const W = g * const_cast<Hull&>(*original_hull).get_reference_mass();
        double const B = g * H.get_displaced_mass();       // buoyancy at this pose (hydrostatic lift)
        double const L_dyn = std::max(0.0, W - B);          // dynamic (pressure) lift, parameter-free
        double const tau_run = std::max(0.0, -degrees(H.get_pitch()));
        if (pdyn_dbutt_deg_ < -1e8) {
            double v = 0.0;
            if (cfg.pdyn_kwarp != 0.0) {
                Hull& O = const_cast<Hull&>(*original_hull);
                double const x_ap = O.get_ap();
                double const d = O.get_quarter_buttock_delta(
                    x_ap + 0.5 * (O.get_fp() - x_ap), 0.25 * b);
                if (d == d) v = std::max(0.0, std::min(5.0, degrees(d)));
            }
            pdyn_dbutt_deg_ = v;
        }
        double const delta = std::max(0.0, cfg.pdyn_delta0 + cfg.pdyn_kwarp * pdyn_dbutt_deg_);
        double const tau_eff = std::min(30.0, tau_run + delta);
        double const D_law = L_dyn * std::tan(radians(tau_eff));
        return (1.0 - r) * morabito_x + r * D_law;
    }
    void set_radojcic_total_envelope(bool on) { cfg.radojcic_total_envelope = on; }
    /// Deadrise gate on the NSS residuary excess (see cfg.radojcic_residuary_beta_gate).
    void set_radojcic_residuary_beta_gate(bool on) { cfg.radojcic_residuary_beta_gate = on; }

    // Raw Radojcic NSS regression core. Returns RT/Delta (total resistance / weight) for the
    // hard-chine planing-craft regression, or -1.0 if disabled / non-planing / out-of-scope.
    // out_RR_D = residuary/Delta (RT/Delta minus the reference friction line); out_W = weight
    // (vol*rho*g); out_grad = low-Fn phase-in taper; out_fallback = slender-hull (use M-S instead).
    double radojcic_rt_over_disp(double speed, double& out_W, double& out_grad,
                                 double& out_RR_D, bool& out_fallback) const {
        out_W = 0.0; out_grad = 0.0; out_RR_D = 0.0; out_fallback = false;
        if (!cfg.radojcic_nss_on || speed < 1e-2 || !hull) {
            return -1.0;
        }
        Hull& H = const_cast<Hull&>(*hull);
        if (!H.is_planing_form()) {
            return -1.0;
        }

        Hull& OH = const_cast<Hull&>(*original_hull);
        double const g = env->get_gravity();
        double const rho = env->get_density();
        double const vol = std::max(1e-9, OH.get_reference_displacement());
        double const vol13 = std::pow(vol, 1.0 / 3.0);
        double const Fnv = H.get_fn(speed, vol13);

        // Reference the DESIGN-float waterline geometry, not original_hull's current pose: the
        // caller (e.g. warm-start continuation) may leave original_hull at the previous running
        // attitude, which otherwise leaks the equilibrium-solver seed into this residuary and hence
        // the total (cold != warm at planing inception, where the heave equilibrium is multi-valued).
        // Same fix the Mercier-Savitsky residuary below applies; here via the captured reference.
        double const L = std::max(1e-3, OH.get_reference_length_wl());
        double const bpx = std::max(1e-3, OH.get_reference_beam_chine());
        double const lcg = OH.get_cg(false).x; // distance forward of transom

        double const lp_v13 = L / vol13;

        // Fallback to Mercier-Savitsky for highly slender hulls (e.g. Taunton) where Radojcic diverges
        if (lp_v13 > 8.0) {
            out_fallback = true;
            return -1.0;
        }

        double const lp_bpx = L / bpx;
        double const lcg_lp = lcg / L;
        
        double X[4] = {lp_v13, lp_bpx, lcg_lp, Fnv};
        double Pk[4] = {0.3214286, 0.2912621, 16.6666698, 0.2893891};
        double Rk[4] = {-1.5571429, -1.0597087, -5.4166677, -0.2770096};
        
        double A[9][4] = {
            {1.6841110, 3.0981180, 1.5282840, -16.4245400},
            {2.4092060, 0.3142931, 0.3908146, 0.2556609},
            {-1.9152880, -7.1345370, -3.6484370, -0.6078711},
            {-1.4699330, -10.6498300, -0.0225162, -0.5634111},
            {4.1680880, -4.5602380, 0.6937186, 2.7010170},
            {1.3626750, 0.1168425, 0.7784404, -10.0020100},
            {0.5479282, 2.5084270, -1.1066880, -3.2642590},
            {-0.9420888, -0.2510580, 0.0049670, -7.9977910},
            {-3.1893870, 1.2750320, -0.5223626, -3.9322540}
        };
        double a_j[9] = {3.1890700, -1.8597140, 1.1340150, 10.0365500, -2.0246550, -0.6886486, -0.3612035, 2.1865540, 4.7179180};
        
        double B[6][9] = {
            {-0.4082681, 0.0891881, 1.4520660, 0.0838715, 1.1443510, -4.0656230, 1.7600280, 3.3407880, 2.1149690},
            {0.6056641, 1.8613140, 0.3999820, 0.1692302, 0.0305501, 2.5074690, -0.4054057, -0.8660036, 1.2468480},
            {-0.1796350, -2.7916050, -2.9371210, 0.3929720, -1.9068460, -6.4689190, -0.7239615, -2.9223100, -2.3012630},
            {-5.3774030, -4.4911580, -2.2997850, -1.4914010, -3.7586970, -7.1752740, 0.5064505, 8.1555560, 5.5388670},
            {-0.2961807, -4.7618480, -2.7215400, -0.5896583, 2.0540800, 0.1658212, -1.8644930, -1.6496020, -8.5582700},
            {-6.2256920, 2.9541990, -2.5990310, 0.7695088, -3.3201360, -4.2225760, -1.7524090, 7.9941520, -3.7219060}
        };
        double b_i[6] = {-2.0033580, -3.5325410, 3.3583560, -4.1627820, 2.9999680, -0.1329643};
        
        double C1[6] = {-2.2633980, -7.7446910, 0.2705126, 4.3702590, 7.6477400, 0.0098969};
        double c1 = 3.2038660;
        
        double L1 = 6.4766839;
        double G1 = -0.2948834;
        
        auto sig = [](double x) { return 1.0 / (1.0 + std::exp(-x)); };
        
        double sum_i = 0.0;
        for (int i = 0; i < 6; i++) {
            double sum_j = 0.0;
            for (int j = 0; j < 9; j++) {
                double sum_k = 0.0;
                for (int k = 0; k < 4; k++) {
                    sum_k += A[j][k] * (Pk[k] * X[k] + Rk[k]);
                }
                sum_j += B[i][j] * sig(a_j[j] + sum_k);
            }
            sum_i += C1[i] * sig(b_i[i] + sum_j);
        }
        
        double RT_D = (sig(c1 + sum_i) - G1) / L1;
        
        // Subtract reference friction to obtain residuary resistance RR/Delta.
        // Design-pose B/T (like L/bpx above): H is at the running attitude, which would leak the
        // solver seed into the wetted-surface term (cold!=warm). Use the captured reference geometry.
        double const BT = OH.get_reference_beam_wl() / std::max(1e-4, OH.get_reference_draft());
        double const S_v23 = 2.262 * std::sqrt(L / vol13) * (1.0 + 0.046 * BT + 0.00287 * BT * BT);
        // Reference-craft Reynolds Re = V_ref*L_ref/nu for the 100,000-lb (salt water, 64 lb/ft^3)
        // reference craft at this Fn_vol: V_ref = Fnv*sqrt(g*vol13_ref), L_ref = (L/vol13)*vol13_ref,
        // vol13_ref = (100000/64)^(1/3), g = 32.2 ft/s^2, nu = 1.2817e-5 ft^2/s. Collapsing the
        // vol13_ref powers gives Re = Fnv*(L/vol13)*sqrt(32.2*100000/64)/nu. The prior form used
        // sqrt(L/vol13) and (100000/64)^(1/3) under the root -> ~27x too low -> CF_ref ~1.7x too high
        // -> residuary under-added (~1% of Delta at mid-Fnv). SOLVER_FIX_PLAN.md P3 (M&S 1973 / R-1667).
        double const Rn_ref = Fnv * (L / vol13)
            * std::sqrt(32.2 * 100000.0 / 64.0) / 1.2817e-5;
        double const CF_ref = 0.075 / sq(std::log10(std::max(1e3, Rn_ref)) - 2.0);  // ~Schoenherr
        
        double RR_D = RT_D - CF_ref * S_v23 * sq(Fnv) * 0.5;
        if (RR_D < 0.0) RR_D = 0.0;

        // Radojcic NSS applies over Fn_vol 1.1 to 4.2; taper off slightly at low end
        double g_rad = 1.0;
        if (Fnv < 0.9) g_rad = 0.0;
        else if (Fnv < 1.1) {
            double t = (Fnv - 0.9) / 0.2;
            g_rad = t * t * (3.0 - 2.0 * t);
        }

        out_W = vol * rho * g;
        out_grad = g_rad;
        out_RR_D = RR_D;
        return RT_D;
    }

    // Residuary contribution (RT/Delta minus reference friction), in Newtons. The shortfall-gated
    // form used by the Auto blend.
    double get_radojcic_nss_residuary(double speed) const {
        double W, grad, RR_D; bool fb;
        double const RT_D = radojcic_rt_over_disp(speed, W, grad, RR_D, fb);
        if (fb) return get_mercier_savitsky_residuary(speed);
        if (RT_D < 0.0) return 0.0;
        return grad * RR_D * W * (michell.is_demihull() + 1);
    }

    // Raw Radojcic NSS TOTAL resistance prediction, in Newtons (RT/Delta * weight). The regression's
    // native output; used as a high-Fn total-resistance envelope where the component model's residuary
    // closures have tapered to zero. Slender-hull fallback returns 0 (no Radojcic total available).
    double get_radojcic_nss_total(double speed) const {
        double W, grad, RR_D; bool fb;
        double const RT_D = radojcic_rt_over_disp(speed, W, grad, RR_D, fb);
        if (fb || RT_D < 0.0) return 0.0;
        return grad * RT_D * W * (michell.is_demihull() + 1);
    }


    Resistance(std::shared_ptr<Hull> ship_hull)
    : original_hull(ship_hull)
    , hull(std::make_shared<Hull>(*ship_hull)) // make a copy for this instance
    , michell(hull)
    , holtrop(hull)
    , localflow(hull)
    , morabito(hull, ship_hull->get_environment())
    , savitsky(hull, ship_hull->get_environment())
    , sectional(hull, ship_hull->get_environment())
    {
        env = hull->get_environment();
        // The trim=-pitch convention, buoyancy-torque sign, profilegrid transom scan, and the
        // localflow/thinship BL-slope frames all silently assume +X-forward (get_fwd()==1). Fail
        // loudly rather than produce silently wrong forces for an X-aft hull.
        if (hull->get_fwd() != 1) {
            throw std::runtime_error("Resistance requires a +X-forward hull (Hull fwd==1); "
                                     "the solver's sign/frame conventions assume it.");
        }
    }

    ~Resistance() {
//        std::cout << "Delete resistance instance" << std::endl;
    }

    /// Sync the internal hull copy's attitude (heave and pitch) from the original hull.
    /// Resistance holds a deep copy to allow find_equilibrium to mutate attitude
    /// without affecting the caller's hull, but this bypasses manual attitude
    /// updates on the original hull unless synced.
    void sync_attitude_from_original() {
        hull->set_heave(original_hull->get_heave());
        hull->set_pitch(original_hull->get_pitch());
    }

    std::shared_ptr<Hull> get_hull() const {
        return hull;
    }

    void set_params(double p1, double p2, double p3, double p4, double p5) {
        michell.set_params(p1, p2, p3, p4, p5);
        localflow.set_params(p1, p2, p3, p4, p5);
    }

    void set_deflection_table(const std::vector<double>& xs, const std::vector<double>& ys) {
        michell.set_deflection_table(xs, ys);
        localflow.set_deflection_table(xs, ys);
    }

    void set_deflection_field(const std::vector<double>& field) {
        michell.set_deflection_field(field);
        localflow.set_deflection_field(field);
    }

    void set_deflection_fields(
        const std::vector<double>& michell_field,
        const std::vector<double>& localflow_field
    ) {
        michell.set_deflection_field(michell_field);
        localflow.set_deflection_field(localflow_field);
    }

    void set_phase_model(double alpha, double gamma) {
        michell.set_phase_model(alpha, gamma);
    }

    ProfileGrid* get_grid() const {
        return michell.get_grid();
    }

    ProfileGrid* get_localflow_grid() const {
        return localflow.get_grid();
    }

    void set_grid_precision(ResistanceGridPrecision precision) {
        mesh_precision = precision;
    }

    bool prepare() {

        if (!is_ok()) {
            return false;
        }

        if (mesh_precision == ResistanceGridPrecision::Finest) {
            michell.set_grid(253, 53);
            localflow.set_grid(41, 17);
//            localflow.prepare(9);
        } else if (mesh_precision == ResistanceGridPrecision::Fine) {
            michell.set_grid(153, 33);
            localflow.set_grid(33, 13);
//            localflow.prepare(6);
        } else {
            michell.set_grid(67, 33);
            localflow.set_grid(31, 11);
//            localflow.prepare(3);
        }

        michell.fix_stuff(do_fixes, do_fixes, false);

        return true;
    }

    void set_transom_method(TransomShapeMethod hollow, TransomShapeMethod dryness) {
        michell.set_transom_method(hollow, dryness);
    }

    /// Raw thinship near-field (drag, lift, torque) triple — Batch 3 diagnostic access,
    /// single-hull convention (no demihull factor), matching the internal attitude use.
    glm::dvec3 get_michell_drag_lift_torque(double speed) {
        if (!is_ok()) { return {0, 0, 0}; }
        return michell.get_drag_lift_torque(speed);
    }

    /// Near-field / far-field Michell wave-drag ratio (Batch 3 acceptance diagnostic).
    double get_nearfar_drag_ratio(double speed) {
        if (!is_ok()) { return 0; }
        return michell.get_nearfar_drag_ratio(speed);
    }

    void set_demihull(double separation) {
        michell.set_demihull(separation);
    }

    bool is_demihull() {
        return michell.is_demihull();
    }

    double get_friction_coef(double speed, double length=0.0) {
        if (!is_ok()) { return 0; }
        if (length <= 0.0) {
            length = hull->get_length_ref();
        }
        return 0.075 / sq(std::log10(hull->get_rn(speed, length)) - 2);
    }

    double get_roughness_coef(double speed) {

        if (!is_ok()) {
            return 0;
        }

        if (hull->get_roughness() < 0.0) {
            return 0.0;
        }

        double dCf = 0.000125 + 0.044 * (std::pow(hull->get_roughness()/hull->get_length_ref(), 1.0/3.0) - 10 * std::pow(hull->get_rn(speed), -1.0/3.0));
        return std::max(dCf, 0.0);
    }

    /// Smooth Fn_vol-blend weight: 0 = pure Unified (displacement), 1 = pure planing.
    /// Cubic smoothstep over [Fn_min, Fn_max]; exactly 0/1 outside this range so
    /// planing equations never leak garbage into the displacement regime. The blend
    /// is also gated by hull geometry: pure displacement forms (round bilge, no
    /// chine at the waterline) return 0 regardless of Fn_vol, since the Savitsky /
    /// Morabito equations require a prismatic planing-style geometry.
    /// Predicted mean wetted length-to-beam ratio λ used by the Morabito bottom-pressure
    /// integral at the current attitude (validation hook; cf. Fridsma ℓ_m/b).
    double get_planing_lambda(double speed) {
        return morabito.get_lambda(speed);
    }

    /// Geometric mean wetted length-to-beam ratio (keel/chine waterline intersections).
    double get_planing_lambda_geometric() {
        return morabito.get_geometric_lambda();
    }

    /// Mean stagnation/spray-root line angle (rad) setting the wetted-planform shape.
    double get_planing_stagnation_angle(double speed) {
        return morabito.get_stagnation_angle(speed);
    }

    /// Cap the planing wetted length-to-beam ratio λ (both the bottom-pressure planform
    /// and the friction area) at a physical bound. ≥12 disables. Applies to Auto/Morabito.
    void set_planing_lambda_cap(double cap) {
        morabito.set_lambda_cap(cap);
        savitsky.set_lambda_cap(cap);
    }

    /// Warp-resolved-trim productionization (opt-in, OFF by default): integrate the bottom
    /// pressure over the TRUE per-strip wetted planform × a pile-up factor, instead of the
    /// Savitsky-λ triangle. See proto_phase*.py for the closure characterisation.
    void set_true_planform(bool on) { morabito.set_true_planform(on); savitsky.set_true_wetted(on); }
    /// Savitsky (2012) warped-hull effective-attitude corrections (default OFF). mode 1 =
    /// as-published equivalent prism (effective deadrise at the mean-wetted-length station +
    /// quarter-beam buttock trim substituted into the prismatic relations); mode 2 = tau-only
    /// hybrid (the kernel keeps its per-cell local-beta resolution). Applies to the Morabito
    /// pressure kernel, the running-trim planing friction, and the predictive-trim blend.
    /// Warp-guarded: constant-deadrise hulls take the production path byte-identically.
    void set_savitsky2012(int mode, double station_frac = 1.0) {
        morabito.set_savitsky2012(mode, station_frac);
        savitsky.set_savitsky2012(mode, station_frac);
        cfg.savitsky2012_mode = std::max(0, std::min(2, mode));
    }
    /// Decouple the FRICTION wetted length from the Morabito attitude planform: friction over
    /// the true calm-water profile (on, default) omits the spray-root pile-up that wets the
    /// bottom forward of the waterline; off uses Savitsky's residuary-calibrated (longer) λ,
    /// the empirical spray-inclusive wetted length. The attitude/CoP planform is unaffected,
    /// so the physical sinkage/trim are preserved -- the "separate wetted-length scales" lever.
    void set_friction_true_wetted(bool on) { savitsky.set_true_wetted(on); }
    /// Speed-ramped spray-root friction (default off): blend the true wetted length toward the
    /// spray-inclusive Savitsky λ as Fn_vol rises through [lo,hi]; closes the high-Fn friction
    /// deficit while leaving the low/mid-Fn (well-matched) range on the true length.
    void set_friction_cv_ramp(double lo, double hi) { savitsky.set_friction_cv_ramp(lo, hi); }
    /// pdyn Stage 3b: chines-dry friction-width correction (see Savitsky::calc).
    void set_friction_width(bool on) { savitsky.set_friction_width(on); }
    /// Planing-bottom correlation allowance Ca (default 7.0e-4, planing practice). Exposed only to
    /// quantify the MaxSurf benchmark asymmetry (that tool run at zero allowance) via a Ca=0 A/B.
    void set_planing_ca(double ca) { savitsky.set_planing_ca(ca); }
    /// Derived deadrise-scaled spray-root friction (default OFF, parameter-free). Blends the
    /// true-wetted length toward the Savitsky load λ by the lift-efficiency deficit (1-C_Lβ/C_L0);
    /// → 0 at β=0 (do-no-harm on flat hulls), grows with deadrise. Needs friction_true_wetted ON.
    void set_friction_deadrise_scaled(bool on) { savitsky.set_deadrise_scaled_friction(on); }
    /// Local free-surface wave envelope (default OFF): when ON, the Auto wave term becomes
    /// max(Michell, LocalFlow-Noblesse near-field drag), gated by the planing weight, to
    /// supply the warped-hull wave-making the centerplane-linearized Michell integral
    /// under-resolves at the hump and high Fn. Parameter-free. See localfs-wave-envelope.
    void set_localfs_wave_envelope(bool on) { cfg.localfs_wave_envelope = on; }
    /// Pure LocalFlow Noblesse near-field pressure drag (no transom term), full hull --
    /// the parameter-free local free-surface wave-resistance estimate that the envelope uses.
    double get_localflow_drag(double speed) {
        if (!is_ok() || speed < 1e-2) return 0;
        return localflow.get_drag_lift_torque_noblesse(speed).x * (michell.is_demihull() + 1);
    }
    /// LocalFlow Noblesse near-field VERTICAL lift at the current attitude (diagnostic). This is the
    /// lift the planing blend weights at (1 - planing_weight), i.e. switched OUT at full planing in
    /// favour of Morabito -- exposed to test whether it carries any of the warp lift deficit.
    double get_localflow_lift(double speed) {
        if (!is_ok() || speed < 1e-2) return 0;
        return localflow.get_drag_lift_torque_noblesse(speed).y * (michell.is_demihull() + 1);
    }
    void set_morabito_hydrostatic(bool on) { morabito.set_hydrostatic(on); }
    /// 2.5D sectional water-entry attitude closure (default ON; warp-gated). When ON, the
    /// running-attitude planing lift+moment come from the Sectional added-mass kernel (one mono-pinned
    /// scale K) instead of Morabito-dynamic+Wagner-CoP; the Morabito drag path is unchanged.
    void set_sectional_attitude(bool on, double K, double warp_floor) {
        cfg.use_sectional_attitude = on;
        cfg.sectional_delta_mode = false;
        if (K > 0.0) sectional.set_scale(K);
        if (warp_floor >= 0.0) cfg.sectional_warp_floor = warp_floor;   // <0 keeps the default; 0 = ungated
    }
    void set_sectional_delta(bool on, double K, double defl_lo, double defl_hi, double fn_lo, double fn_hi) {
        cfg.use_sectional_attitude = on;
        cfg.sectional_delta_mode = on;
        if (K > 0.0) sectional.set_scale(K);
        if (defl_lo >= 0.0) cfg.sectional_deficit_lo = defl_lo;
        if (defl_hi >  0.0) cfg.sectional_deficit_hi = std::max(cfg.sectional_deficit_lo + 1e-3, defl_hi);
        if (fn_lo >= 0.0) cfg.sectional_fn_lo = fn_lo;
        if (fn_hi >= 0.0) cfg.sectional_fn_hi = fn_hi;
    }
    /// Enable the ITP outer-trim root (see use_itp_root). Default ON (use_itp_root = true).
    void set_itp_root(bool on) { cfg.use_itp_root = on; }
    bool get_itp_root() const { return cfg.use_itp_root; }
    void set_sectional_tloc_cap(double deg) { sectional.set_tloc_cap(deg); }
    void set_sectional_nx(int n) { sectional.set_nx(n); }
    void set_sectional_smooth_frac(double f) { sectional.set_smooth_frac(f); }
    void set_sectional_transom_fix(bool on) { sectional.set_transom_fix(on); }
    void set_sectional_ventilation(bool on) { sectional.set_transom_ventilation(on); }
    void set_sectional_vent_beam_frac(double f) { sectional.set_vent_beam_frac(f); }
    void set_sectional_similarity_model(bool on) { sectional.set_similarity_model(on); }
    /// Crossflow sectional lift term (round-3 #2). coef 0 (default) = OFF, byte-identical.
    void set_sectional_crossflow(double coef) { sectional.set_crossflow(coef); }
    double get_sectional_crossflow() const { return sectional.get_crossflow(); }
    void set_sectional_forward_term(bool on, double cs) { sectional.set_forward_term(on, cs); }
    void set_attitude_wave_blend(bool on) { cfg.attitude_wave_blend = on; }
    bool get_attitude_wave_blend() const { return cfg.attitude_wave_blend; }
    /// Item 4 (F5): scale the Morabito-path Michell wave lift by the immersed-volume fraction
    /// V(pose)/V(design) instead of a blanket fade. Default OFF (byte-identical).
    void set_wave_lift_volume_scale(bool on) { cfg.wave_lift_volume_scale = on; }
    bool get_wave_lift_volume_scale() const { return cfg.wave_lift_volume_scale; }
    void debug_sectional_at(double heave, double pitch, double speed) {}
    
    // Add get_sectional_deficit for the python interface
    double get_sectional_deficit(double speed) {
        // Return the lift-deficit indicator at the current pose
        if (!cfg.use_sectional_attitude) return std::nan("");
        auto const planing_flow = morabito.get_drag_lift_torque(speed);
        auto const sr = sectional.added_mass(speed);
        if (sr.nwet < 5 || !std::isfinite(sr.La)) return std::nan("");
        double sec_lift = sectional.get_scale() * sr.La;
        if (sectional.get_forward_term()) {
            sec_lift += sectional.get_forward_cs() * sr.Lf;
        }
        double const denom = std::max(0.05 * hull->get_displaced_mass() * env->get_gravity(), std::abs(planing_flow.y));
        return (sec_lift - planing_flow.y) / denom;
    }
    void set_dynamic_cop_fraction(double frac) { morabito.set_dynamic_cop_fraction(frac); }
    void set_cop_savitsky(bool on) { morabito.set_cop_savitsky(on); }
    /// Cv-ramp the dynamic-CoP fraction, gated by the hull's fore-aft deadrise spread (warp).
    /// slope=0 (default) is OFF/byte-identical; ref is the neutral beam-Froude Cv. Warped hulls only.
    void set_cop_cv_ramp(double slope, double ref, double warp_floor) {
        morabito.set_cop_cv_ramp(slope, ref, warp_floor); }
    /// Warped-transition lift closure (Stage-B prototype, default-OFF). Raise the aft wetted area /
    /// pile-up through the Cv transition window on warped hulls. lift_amp drives the squat/phase fix
    /// (Hook 1: lift + induced drag); moment_amp drives the bow-up CoP shift (Hook 2), DECOUPLED so
    /// the squat fix need not over-trim. Both 0 (default) = OFF/byte-identical.
    void set_warp_lift_closure(double lift_amp, double moment_amp, double cv_lo, double cv_hi) {
        morabito.set_warp_lift_closure(lift_amp, moment_amp, cv_lo, cv_hi); }
    /// Use the chine-beam Froude number FnB (planing scale) for the displacement->planing blend
    /// instead of the volumetric Fn_vol. Default OFF (Fn_vol). Set the band via set_planing_blend.
    void set_blend_use_fnb(bool on) { cfg.blend_use_fnb = on; }

    /// Restore the pre-2026-06-25 calibrated risen-pose planing method (Savitsky predictive trim,
    /// full-weight Morabito on the Savitsky-lambda planform, hydrostatic included, integral CoP) for
    /// the paper's physical-vs-calibrated comparison. The physics-first config is the default.
    void set_planing_legacy() {
        cfg.savitsky_trim_weight = 1.0;
        set_true_planform(false);
        set_morabito_hydrostatic(true);
        set_dynamic_cop_fraction(-1.0);
    }
    void set_planing_pileup(double f) { morabito.set_pileup(f); }
    void set_planing_pileup_auto(bool on) { morabito.set_auto_pileup(on); }
    void set_planing_friction_factor(double f) { savitsky.set_friction_factor(f); }

    // Cubic-smoothstep Fn_vol blend weight on an explicit [fn_min, fn_max] band.
    double planing_weight_band(double speed, double fn_min, double fn_max) const {
        if (speed < 1e-2 || !hull) {
            return 0.0;
        }
        if (!const_cast<Hull&>(*hull).is_planing_form()) {
            return 0.0;
        }
        // Blend Froude number length scale: chine beam (FnB = V/sqrt(g·b), the planing/Savitsky
        // scale) when blend_use_fnb, else the design displacement^(1/3) (Fn_vol). Use the reference
        // (design) displacement, not the running one: as a planing hull rises get_displacement()
        // collapses, which would make the blend weight swing with attitude and trap the solve.
        double const len_scale = cfg.blend_use_fnb
            ? std::max(1e-6, const_cast<Hull&>(*hull).get_beam_chine())
            : std::pow(std::max(1e-12, const_cast<Hull&>(*hull).get_reference_displacement()), 1.0/3.0);
        double const Fn = const_cast<Hull&>(*hull).get_fn(speed, len_scale);
        double const span = std::max(1e-6, fn_max - fn_min);
        double t = (Fn - fn_min) / span;
        if (t <= 0.0) return 0.0;
        if (t >= 1.0) return 1.0;
        return t * t * (3.0 - 2.0 * t);
    }

    // Resistance-component blend weight (calibrated band).
    double planing_weight(double speed) const {
        return planing_weight_band(speed, cfg.planing_blend_fn_min, cfg.planing_blend_fn_max);
    }

    // Running-attitude blend weight (decoupled band; see planing_blend_attitude_fn_*).
    double planing_weight_attitude(double speed) const {
        return planing_weight_band(speed, cfg.planing_blend_attitude_fn_min, cfg.planing_blend_attitude_fn_max);
    }

    /// Tune the planing blend: Fn_vol below fn_min uses pure displacement, above fn_max pure planing.
    void set_planing_blend(double fn_min, double fn_max) {
        cfg.planing_blend_fn_min = fn_min;
        cfg.planing_blend_fn_max = std::max(fn_min + 1e-3, fn_max);
    }

    /// Tune the decoupled ATTITUDE blend (find_equilibrium lift+trim only); default = resistance band.
    void set_planing_blend_attitude(double fn_min, double fn_max) {
        cfg.planing_blend_attitude_fn_min = fn_min;
        cfg.planing_blend_attitude_fn_max = std::max(fn_min + 1e-3, fn_max);
    }

    /// Tune the pre-planing residuary hump (Lever 2). coef = peak R/Delta added at the
    /// hump centre (midpoint of [fn_lo, fn_hi]); coef = 0 disables it.
    void set_preplaning_hump(double coef, double fn_lo = 0.6, double fn_hi = 1.9) {
        cfg.preplaning_hump_coef = coef;
        cfg.preplaning_fn_lo = fn_lo;
        cfg.preplaning_fn_hi = std::max(fn_lo + 1e-3, fn_hi);
    }

    /// F7 effective buoyancy: remove frac * w_plan * s(CvB) of the pose hydrostatic lift+torque
    /// (dry-transom / free-surface-depression over-count; frac = 0.5 is the physics value, the
    /// classic planing half-buoyancy). CvB smoothstep [cv_lo, cv_hi] phases the removal with
    /// transom ventilation. frac <= 0 disables (default; byte-identical).
    void set_effective_buoyancy(double frac, double cv_lo = 0.1, double cv_hi = 3.5,
                                bool torque = true) {
        cfg.eff_buoy_frac = frac;
        cfg.eff_buoy_cv_lo = cv_lo;
        cfg.eff_buoy_cv_hi = std::max(cv_lo + 1e-3, cv_hi);
        cfg.eff_buoy_torque = torque;
    }
    /// Lever 1 experiment: weight on blending the equilibrium trim toward Savitsky's
    /// predict_trim. 1.0 = current; 0.0 = pure Morabito moment-equilibrium trim.
    void set_savitsky_trim_weight(double w) { cfg.savitsky_trim_weight = w; }
    /// Warp-gated Savitsky trim weight: ramp the trim weight 0 -> w by the fore-aft deadrise spread
    /// (smoothstep over [lo, hi] rad), so prismatic/hard-chine hulls (spread ~ 0, incl. the out-of-sample
    /// prism series) stay at weight 0 and only warped hulls get the partial Savitsky blend. Replaces the
    /// per-harness Python gate. w <= 0 disables the ramp (the flat set_savitsky_trim_weight applies).
    void set_savitsky_trim_weight_warp(double w, double lo, double hi) {
        cfg.savitsky_trim_warp_w = w;
        cfg.savitsky_trim_warp_lo = lo;
        cfg.savitsky_trim_warp_hi = hi;
    }
    /// Drag-couple pitch moment (R*a): add the moment of the horizontal resistance about the CG to
    /// the equilibrium moment balance (the vertical-force-only balance omits it). Parameter-free
    /// (a = CG height above the wetted-surface centroid; R = full ship resistance, frozen by an outer
    /// fixed-point). warp_floor 0 = ungated (universal physics); >0 = applied only on warped hulls
    /// (fore-aft deadrise spread > floor, smoothstep over [floor, warp_hi] rad). Default OFF.
    void set_drag_couple(bool on, double warp_floor = 0.0, double warp_hi = 0.15) {
        cfg.drag_couple_on = on;
        cfg.drag_couple_warp_floor = warp_floor;
        cfg.drag_couple_warp_hi = std::max(warp_floor + 1e-6, warp_hi);
    }

    /// Diagnostic: Savitsky predictive trim (deg) at a given effective deadrise
    /// (beta_override in radians; <0 = hull area-mean). Pins the design pose first, since
    /// predict_trim's LCG is only valid there, then RESTORES the caller's pose: this is a
    /// read-only diagnostic, and leaving the hull at the design pose poisons every
    /// pose-dependent getter evaluated after it (observed: Auto total 52.9 -> 202.8 N on a
    /// converged Taunton equilibrium). Used to invert for the deadrise a warped hull would
    /// need to reproduce its measured trim (warp-resolved-trim study).
    double get_savitsky_trim(double speed, double beta_override = -1.0) {
        double const heave0 = hull->get_heave(), pitch0 = hull->get_pitch();
        hull->set_heave(0.0); hull->set_pitch(0.0);
        double const trim = savitsky.predict_trim(speed, beta_override);
        hull->set_heave(heave0); hull->set_pitch(pitch0);
        return trim;
    }

    // ---- Shared residual core: the single source of truth for the lift/pitch-moment balance ----
    // find_equilibrium's inner residual (eval_at), the benchmark residual_at, and the CoP/force-
    // balance diagnostics (debug_force_balance, cop_test_at_trim, cop_components_at_trim) all route
    // through this one function, so the residual can never silently drift between them (audit F5).
    // The caller selects the displacement->planing blend band (the resistance band that the force
    // assembly uses vs. the decoupled attitude band that the solver uses) and the sectional-closure
    // mode; the blend weight is evaluated AFTER the pose is set, exactly mirroring each original call
    // site, so every caller is reproduced bit-for-bit.
    enum class BlendBand { Resistance, Attitude };
    enum class SectionalMode { None, Replacement, Delta };
    struct DeltaState { double locked_t = -1.0; double captured_t = 0.0; bool engaged = false; };
    struct PoseForces {   // field order matches debug_force_balance's 13-element vector
        double static_lift, static_torque, wave_lift, wave_torque,
               noblesse_lift, noblesse_torque, morabito_lift, morabito_torque,
               w_plan, lift_press, torque_press, total_lift, total_torque;
    };
    PoseForces residual_core(double heave, double pitch, double speed,
                             BlendBand band, bool use_sec,
                             SectionalMode mode, DeltaState* dstate) {
        hull->set_heave(heave);
        hull->set_pitch(pitch);
        prepare();
        double const g = env->get_gravity();
        double const mass = original_hull->get_reference_mass();
        double const W = g * mass;
        double static_lift = g * (hull->get_displaced_mass() - mass);
        double static_torque = g * hull->get_displaced_mass() * (hull->get_cb().x - hull->get_cg(true).x);
        auto const wave_flow = michell.get_drag_lift_torque(speed);
        auto const local_flow = localflow.get_drag_lift_torque_noblesse(speed);
        // ATTITUDE band (decoupled from the resistance blend) keeps buoyancy in the balance through
        // the hump so the hull squats instead of planing out early; the RESISTANCE band drives the
        // force assembly. Evaluated here (post-prepare) so it sees the current pose, as the originals did.
        double const w_plan = (band == BlendBand::Attitude) ? planing_weight_attitude(speed)
                                                            : planing_weight(speed);
        // F7 effective buoyancy: remove the dry-transom/free-surface-depression over-count from
        // the hydrostatic terms (see cfg.eff_buoy_frac). Uses the caller's blend band via w_plan,
        // so the displacement regime and non-planing hulls are byte-identical by construction.
        if (cfg.eff_buoy_frac > 0.0 && w_plan > 0.0) {
            double const cvb = speed / std::sqrt(g * std::max(1e-6, hull->get_beam_chine()));
            double const espan = std::max(1e-6, cfg.eff_buoy_cv_hi - cfg.eff_buoy_cv_lo);
            double et = (cvb - cfg.eff_buoy_cv_lo) / espan;
            et = (et <= 0.0) ? 0.0 : (et >= 1.0 ? 1.0 : et * et * (3.0 - 2.0 * et));
            double const escale = 1.0 - cfg.eff_buoy_frac * w_plan * et;
            static_lift -= (1.0 - escale) * g * hull->get_displaced_mass();
            if (cfg.eff_buoy_torque) static_torque *= escale;
        }
        double lift_press = local_flow.y;
        double torque_press = local_flow.z;
        glm::dvec3 planing_flow(0.0, 0.0, 0.0);
        bool sec_now = false;   // sectional kernel carried the planing forces this eval (replacement mode)
        if (w_plan > 0.0) {
            // Base planing forces: Morabito dynamic lift + Wagner-3/4 CoP (the calibrated prism law).
            planing_flow = morabito.get_drag_lift_torque(speed);
            double plan_lift = planing_flow.y;
            double plan_torque = planing_flow.z;
            if (mode != SectionalMode::None && cfg.use_sectional_attitude && use_sec) {
                auto const sr = sectional.added_mass(speed);
                // cop_a is NaN when |La|<=1e-9 (e.g. ventilation zeroes every station); guard it too, else
                // sec_torque = scale*La*(NaN - cg) poisons the moment residual. Falls back to base Morabito.
                if (sr.nwet >= 5 && std::isfinite(sr.La) && std::isfinite(sr.cop_a)) {
                    double const sec_add = sectional.get_forward_term() ? sectional.get_forward_cs() * sr.Lf : 0.0;
                    // Crossflow lift+moment (round-3 #2): forward-loaded second sectional term, raises the
                    // warped CoP. Lc already carries the coefficient; 0/NaN when OFF -> byte-identical.
                    double const sec_cross = (sectional.get_crossflow() > 0.0 && std::isfinite(sr.Lc)) ? sr.Lc : 0.0;
                    double const sec_lift = sectional.get_scale() * sr.La + sec_add + sec_cross;
                    double sec_torque = sectional.get_scale() * sr.La * (sr.cop_a - hull->get_cg(true).x);
                    if (sec_add > 0.0 && std::isfinite(sr.cop_f))
                        sec_torque += sec_add * (sr.cop_f - hull->get_cg(true).x);
                    if (sec_cross != 0.0 && std::isfinite(sr.cop_c))
                        sec_torque += sec_cross * (sr.cop_c - hull->get_cg(true).x);

                    if (mode == SectionalMode::Delta) {
                        double t = 0.0;
                        if (dstate && dstate->locked_t >= 0.0) {
                            t = dstate->locked_t;
                        } else {
                            double const denom = std::max(0.05 * W, std::abs(planing_flow.y));
                            double const defl = (sec_lift - planing_flow.y) / denom;
                            double const span = std::max(1e-6, cfg.sectional_deficit_hi - cfg.sectional_deficit_lo);
                            t = (defl - cfg.sectional_deficit_lo) / span;
                            t = (t <= 0.0) ? 0.0 : (t >= 1.0 ? 1.0 : t * t * (3.0 - 2.0 * t));  // cubic smoothstep
                            if (dstate) dstate->captured_t = t;
                        }
                        // Transition Froude-window gate: suppress the delta at the displacement hump where
                        // the warped hull under-trims (the deficit alone over-fires there); engage it only
                        // through the semi-planing-to-planing window where the over-trim actually lives.
                        if (t > 0.0 && cfg.sectional_fn_hi > cfg.sectional_fn_lo)
                            t *= planing_weight_band(speed, cfg.sectional_fn_lo, cfg.sectional_fn_hi);
                        if (t > 0.0 && std::isfinite(sec_torque)) {
                            plan_lift   = planing_flow.y + t * (sec_lift   - planing_flow.y);
                            plan_torque = planing_flow.z + t * (sec_torque - planing_flow.z);
                            if (dstate) dstate->engaged = true;
                        }
                    } else {   // SectionalMode::Replacement
                        plan_lift = sec_lift;
                        plan_torque = sec_torque;
                        // kshape: CG-acting lift relief (cfg.sec_relief_frac). Lift only --
                        // plan_torque untouched = the shed acts at the CG (moment-neutral),
                        // which is what the equilibrium K-ladder demands (K_ds0 falls to
                        // 0.35-0.7 at top while K_dt0 stays ~1.2; fixes_ab/kshape).
                        if (cfg.sec_relief_frac > 0.0) {
                            double const cvb = speed / std::sqrt(g * std::max(1e-6, hull->get_beam_chine()));
                            double const rspan = std::max(1e-6, cfg.sec_relief_cv_hi - cfg.sec_relief_cv_lo);
                            double rt = (cvb - cfg.sec_relief_cv_lo) / rspan;
                            rt = (rt <= 0.0) ? 0.0 : (rt >= 1.0 ? 1.0 : rt * rt * (3.0 - 2.0 * rt));
                            plan_lift -= cfg.sec_relief_frac * rt * sec_lift;
                        }
                        if (dstate) dstate->engaged = true;
                        sec_now = true;
                    }
                }
            }
            lift_press = (1.0 - w_plan) * local_flow.y + w_plan * plan_lift;
            torque_press = (1.0 - w_plan) * local_flow.z + w_plan * plan_torque;
        }
        // When the 2.5D sectional water-entry kernel carries the planing forces (replacement mode,
        // warped hulls at planing), it IS the planing lift+moment — the thin-ship Michell wave-lift
        // would double-count. Fade the Michell vertical lift/moment out by the planing weight so the
        // displacement regime is untouched. Inert unless sectional replacement engaged (sec_now).
        double wave_y = wave_flow.y, wave_z = wave_flow.z;
        // physfix C6: the same double-count exists on the MORABITO planing path (mono/warp1: no
        // sec_now). At full planing the Morabito dynamic lift is calibrated to carry the entire
        // load, so the un-faded Michell near-field wave lift is added on top (warp1 carries
        // 0.04-0.18 W of it after Batch 3's causal kernel -> the flagged 4.1->13.3 mm over-rise).
        // Fade it by the planing weight exactly as the sec_now branch does, gated by the
        // attitude_wave_blend flag (default OFF). Applied once (|| not double-fade); w_plan=0 in
        // the displacement regime -> byte-identical there and on already-faded sectional hulls.
        if (sec_now || cfg.attitude_wave_blend) {
            wave_y *= (1.0 - w_plan); wave_z *= (1.0 - w_plan);
        } else if (cfg.wave_lift_volume_scale && w_plan > 0.0) {
            // Item 4 (F5): immersion-scaled wave lift on the Morabito planing path. vfrac = 1 at
            // design float (byte-identical in the displacement regime) and falls to ~0.25 at full
            // planing, so the physical wave lift fades with the immersed volume instead of being
            // zeroed. sec_now hulls already faded above; this is the mono/warp1 path only.
            double const vfrac = std::min(1.0, hull->get_displaced_mass()
                                              / std::max(original_hull->get_reference_mass(), 1e-9));
            wave_y *= vfrac; wave_z *= vfrac;
        }
        double total_lift = static_lift + wave_y + lift_press;
        double total_torque = static_torque + wave_z + torque_press;
        // Drag-couple pitch moment. For a SELF-PROPELLED planing craft the thrust line sits low
        // (shaft/prop near the keel, below the drag resultant), so the net thrust-drag couple is
        // BOW-UP (increases trim) -- the accelerating-craft squat, not the braking-dive. (Pure
        // tow-through-CG would be bow-down R*a; that under-trims warp3 further -- A/B confirmed.) The
        // bow-up couple is +R*a in the bow-up-positive convention, with a = CG height above the
        // area-weighted wetted-surface centroid (lever between the drag line and the CG). R is frozen
        // by the find_equilibrium outer fixed-point (drag_couple_R_); 0 => inert (all other
        // residual_core callers stay couple-free). Centroid cached from get_displaced_mass above.
        if (drag_couple_R_ != 0.0) {
            double const a = hull->get_cg(true).z - hull->get_wetted_centroid().z;
            total_torque += drag_couple_R_ * a;
        }
        // Over-immersion past the sheer makes the hydrostatic integral non-finite; present it as large
        // positive lift so the heave balance steers back up to the physical waterline instead of railing
        // on a NaN (the Series 62 failure; auto_validation/diag_series62_heave_balance.py).
        if (!std::isfinite(total_lift)) { total_lift = OVERIMMERSION_LIFT_MULT * W; total_torque = 0.0; }
        return {static_lift, static_torque, wave_flow.y, wave_flow.z,
                local_flow.y, local_flow.z, planing_flow.y, planing_flow.z,
                w_plan, lift_press, torque_press, total_lift, total_torque};
    }

    /// Warp-resolved-trim study: impose a fixed running trim (deg) and heave-balance the
    /// vertical force (lift = 0), leaving the hull at that attitude so the caller can read
    /// the geometric wetted centroid (Hull::get_wetted_centroid_x) against the CG
    /// (Hull::get_cg_x). Returns the residual heave-balanced pitch MOMENT about the CG in
    /// W*Lwl units (>0 = bow-up = the pressure model wants more trim). Decides whether the
    /// geometric area centroid sits at the CG at the measured trim while the pressure moment
    /// is still bow-up (=> a geometric-CoP closure would shed the over-trim).
    double cop_test_at_trim(double speed, double trim_deg) {
        double const g = env->get_gravity();
        double const mass = original_hull->get_reference_mass();
        double const W = g * mass;
        hull->set_heave(0.0); hull->set_pitch(0.0); prepare();
        double const Lref = std::max(1e-3, hull->get_length_wl());
        double const pitch = -radians(trim_deg);
        double const tol_lift = 1e-3 * W;

        auto eval = [&](double heave, double& lift, double& torque) {
            PoseForces const pf = residual_core(heave, pitch, speed,
                                                BlendBand::Resistance, false, SectionalMode::None, nullptr);
            lift = pf.total_lift;
            torque = pf.total_torque;
        };

        double lift, torque, lo = -0.10 * Lref, hi = 0.10 * Lref;
        eval(lo, lift, torque); double const f_lo = lift;
        eval(hi, lift, torque); double const f_hi = lift;
        double h;
        if (f_lo <= 0.0) {
            h = lo;
        } else if (f_hi >= 0.0) {
            h = hi;
        } else {
            h = 0.0;
            for (int k = 0; k < 40; k++) {
                eval(h, lift, torque);
                if (std::abs(lift) < tol_lift) break;
                if (lift > 0.0) lo = h; else hi = h;
                h = 0.5 * (lo + hi);
            }
        }
        eval(h, lift, torque);             // leave the hull at the heave-balanced attitude
        return torque / (W * Lref);
    }

    /// Stage-A CoP diagnostic: like cop_test_at_trim, but returns the THREE pitch-moment
    /// components about CG at the heave-balanced imposed trim, each in W*Lref units:
    /// {buoyancy (static_torque), Michell wave (wave.z), planing (Noblesse/Morabito blend)}.
    /// Their sum equals cop_test_at_trim(). Lets us localise where the moment residual lives.
    glm::dvec3 cop_components_at_trim(double speed, double trim_deg) {
        double const g = env->get_gravity();
        double const mass = original_hull->get_reference_mass();
        double const W = g * mass;
        hull->set_heave(0.0); hull->set_pitch(0.0); prepare();
        double const Lref = std::max(1e-3, hull->get_length_wl());
        double const pitch = -radians(trim_deg);
        double const tol_lift = 1e-3 * W;

        auto eval = [&](double heave, double& lift,
                        double& t_static, double& t_wave, double& t_press) {
            PoseForces const pf = residual_core(heave, pitch, speed,
                                                BlendBand::Resistance, false, SectionalMode::None, nullptr);
            lift = pf.total_lift;
            t_static = pf.static_torque;
            t_wave = pf.wave_torque;
            t_press = pf.torque_press;
        };

        double lift, ts, tw, tp, lo = -0.10 * Lref, hi = 0.10 * Lref;
        eval(lo, lift, ts, tw, tp); double const f_lo = lift;
        eval(hi, lift, ts, tw, tp); double const f_hi = lift;
        double h;
        if (f_lo <= 0.0) { h = lo; }
        else if (f_hi >= 0.0) { h = hi; }
        else {
            h = 0.0;
            for (int k = 0; k < 40; k++) {
                eval(h, lift, ts, tw, tp);
                if (std::abs(lift) < tol_lift) break;
                if (lift > 0.0) lo = h; else hi = h;
                h = 0.5 * (lo + hi);
            }
        }
        eval(h, lift, ts, tw, tp);
        return glm::dvec3(ts, tw, tp) / (W * Lref);
    }

    /// Toggle the whisker-spray drag term (Savitsky-DeLorme-Datla 2007). On by default.
    void set_whisker_spray(bool on) { cfg.whisker_spray_on = on; }

    /// Apply the whisker-spray term only to warped hulls (see get_whisker_spray_drag).
    void set_whisker_warp_only(bool on) { cfg.whisker_warp_only = on; }

    /// Whisker-spray CvB cap (pdyn campaign; see cfg.wsk_cap_c). c=0 disables (default).
    void set_whisker_cvb_cap(double c, double cv_lo = 4.0, double cv_hi = 5.0) {
        cfg.wsk_cap_c = std::max(0.0, std::min(1.0, c));
        cfg.wsk_cap_cv_lo = cv_lo;
        cfg.wsk_cap_cv_hi = std::max(cv_lo + 1e-3, cv_hi);
    }
    /// Transom head-release term in the propeller-disc axial inflow (default OFF, so every
    /// other consumer is byte-identical). Behind a ventilated transom the free surface drops
    /// onto the transom edge, releasing rho*g*dz of hydrostatic head; the flow under it speeds
    /// up. The thin-ship kernel and the far-field wave spectrum carry no such near-field term,
    /// which is why the inviscid disc reads decelerated where RANS reads accelerated.
    /// amp = 1 is the unscaled Bernoulli prediction. See ThinShip::get_transom_head_release.
    void set_transom_head_release(bool on, double amp = 1.0, double len_scale = 1.0) {
        cfg.transom_hr_on = on;
        cfg.transom_hr_amp = std::max(0.0, amp);
        cfg.transom_hr_len = std::max(1e-3, len_scale);
    }

    /// Diagnostics for the transom head-release term; see ThinShip for the definitions.
    std::vector<double> get_transom_hollow_params(double speed) {
        return michell.get_transom_hollow_params(speed);
    }
    double get_transom_hollow_drop(double speed, double X, double Y, double len_scale = 1.0) {
        return michell.get_transom_hollow_drop(speed, X, Y, len_scale);
    }
    double get_transom_head_release(double speed, double X, double Y, double Z,
                                    double len_scale = 1.0) {
        return michell.get_transom_head_release(speed, X, Y, Z, len_scale);
    }

    /// Force the fully-turbulent spray friction line (colleague round-2 B0 bounding diagnostic).
    /// off (default) = the published laminar/transitional line, byte-identical.
    void set_whisker_spray_turbulent(bool on) { cfg.whisker_spray_turbulent = on; }
    bool get_whisker_spray_turbulent() const { return cfg.whisker_spray_turbulent; }
    /// Feed the spray-root (forward-section) deadrise to the whisker closed form (colleague B1).
    /// off (default) = the aft-half mean, byte-identical; warp-gated so prisms are unaffected.
    void set_whisker_spray_root(bool on) { cfg.whisker_spray_root = on; }
    bool get_whisker_spray_root() const { return cfg.whisker_spray_root; }
    /// Fn-gate (Fn_vol smoothstep window) over which the spray-root deadrise migrates from the
    /// aft-mean to the forward section. hi <= lo disables the gate (full beta_F at all speeds).
    void set_whisker_spray_root_fn_gate(double lo, double hi) {
        cfg.whisker_spray_root_fn_lo = lo; cfg.whisker_spray_root_fn_hi = hi; }
    void set_whisker_spray_turb_warp(bool on) { cfg.whisker_spray_turb_warp = on; }
    bool get_whisker_spray_turb_warp() const { return cfg.whisker_spray_turb_warp; }

    /// physfix C2: use the wetted-region local mean deadrise (not the aft-half static mean) as the
    /// whisker-spray deadrise. Parameter-free; byte-identical on constant-deadrise hulls.
    void set_whisker_local_beta(bool on) { cfg.whisker_local_beta = on; }
    bool get_whisker_local_beta() const { return cfg.whisker_local_beta; }

    /// Toggle the above-water air resistance (windage) term. Off by default; when on it is
    /// added to the total for every resistance method. See get_air_resistance.
    void set_air_resistance(bool on) { cfg.air_resistance_on = on; }
    bool get_air_resistance_on() const { return cfg.air_resistance_on; }

    /// Set the windage drag coefficient C_AA (default 0.8). Applied to the above-water
    /// transverse projected area A_VT in R_AA = 1/2 rho_air C_AA A_VT V^2.
    void set_air_drag_coef(double c) { cfg.air_drag_coef = std::abs(c); }
    double get_air_drag_coef() const { return cfg.air_drag_coef; }

    /// Above-water air resistance (windage) in N at the given speed:
    ///   R_AA = 1/2 * rho_air * C_AA * A_VT * V^2,
    /// with rho_air from the Environment, C_AA = air_drag_coef, and A_VT the above-water
    /// transverse projected (frontal) area from Hull::get_frontal_area (full-beam; a user
    /// override via Hull::set_frontal_area takes precedence). Computed from the design pose
    /// (original_hull) so it is a clean function of speed and does not couple into the
    /// running-attitude solve. Always returns the physical value (independent of the toggle)
    /// so the magnitude can be inspected; the toggle only governs whether get_total_resistance
    /// adds it. Returns 0 if there is no above-water area.
    double get_air_resistance(double speed) const {
        if (!is_ok() || speed < 1e-2) {
            return 0.0;
        }
        double const rho_air = env->get_air_density();
        double const A_vt = const_cast<Hull&>(*original_hull).get_frontal_area();
        if (rho_air <= 0.0 || A_vt <= 0.0) {
            return 0.0;
        }
        return 0.5 * rho_air * cfg.air_drag_coef * A_vt * speed * speed;
    }

    /// Planing-friction wetted-length mode: true (default) = running-trim λ; false =
    /// Savitsky's predictive-trim λ (residuary-calibrated). Experiment lever.
    void set_planing_friction_running_trim(bool on) { cfg.planing_friction_running_trim = on; }

    /// F8 sectional pressure-drag consistency (see cfg.sectional_pressure_drag). Returns the
    /// planing pressure drag the assembly misses when the sectional kernel carries the lift:
    /// w_plan * max(0, sec_lift - morabito_lift) * tan(tau_run). 0 unless the toggle is on,
    /// the hull passes the sectional warp gate, and the kernel has a valid wetted run.
    double get_sectional_pressure_drag(double speed) {
        if (!cfg.sectional_pressure_drag || !cfg.use_sectional_attitude || speed < 1e-2 || !hull) {
            return 0.0;
        }
        double const w_plan = planing_weight(speed);
        if (w_plan <= 0.0) {
            return 0.0;
        }
        Hull& H = const_cast<Hull&>(*hull);
        // Same warp gate as the sectional attitude path (residual_at / find_equilibrium):
        // prisms and the OOS series stay byte-identical.
        double const x_ap = H.get_ap(), x_fp = H.get_fp();
        double const bA = H.get_deadrise_at_x(x_ap + WARP_SAMPLE_AFT_FRAC * (x_fp - x_ap));
        double const bF = H.get_deadrise_at_x(x_ap + WARP_SAMPLE_FWD_FRAC * (x_fp - x_ap));
        if (!(bF > 0.0 && bA > 0.0) || (bF - bA) <= cfg.sectional_warp_floor) {
            return 0.0;
        }
        double const tau = -hull->get_pitch();   // set_pitch(-trim): bow-up trim = -pitch
        if (tau <= 0.0) {
            return 0.0;
        }
        auto const sr = sectional.added_mass(speed);
        if (sr.nwet < 5 || !std::isfinite(sr.La)) {
            return 0.0;
        }
        double const sec_add = sectional.get_forward_term() ? sectional.get_forward_cs() * sr.Lf : 0.0;
        double const sec_cross = (sectional.get_crossflow() > 0.0 && std::isfinite(sr.Lc)) ? sr.Lc : 0.0;
        double sec_lift = sectional.get_scale() * sr.La + sec_add + sec_cross;
        double const mor_lift = morabito.get_drag_lift_torque(speed).y;
        if (!std::isfinite(sec_lift) || !std::isfinite(mor_lift)) {
            return 0.0;
        }
        double const g = env->get_gravity();
        double const cvb = speed / std::sqrt(g * std::max(1e-6, H.get_beam_chine()));
        // kshape consistency: price the lift the balance actually carries -- apply the same
        // CG relief residual_core's Replacement branch applies (cfg.sec_relief_frac).
        if (cfg.sec_relief_frac > 0.0) {
            double const rspan = std::max(1e-6, cfg.sec_relief_cv_hi - cfg.sec_relief_cv_lo);
            double rt = (cvb - cfg.sec_relief_cv_lo) / rspan;
            rt = (rt <= 0.0) ? 0.0 : (rt >= 1.0 ? 1.0 : rt * rt * (3.0 - 2.0 * rt));
            sec_lift -= cfg.sec_relief_frac * rt * sec_lift;
        }
        double charge = w_plan * std::max(0.0, sec_lift - mor_lift) * std::tan(tau);
        // kshape: only the DYNAMIC share of the carried lift pays L*tan(tau) pressure drag;
        // at hump CvB the sectional supply substitutes buoyant load (cfg.sec_pdrag_dynshare).
        // Same smoothstep partition as the F7 effective-buoyancy ramp.
        if (cfg.sec_pdrag_dynshare && charge > 0.0) {
            double const dspan = std::max(1e-6, cfg.eff_buoy_cv_hi - cfg.eff_buoy_cv_lo);
            double dt = (cvb - cfg.eff_buoy_cv_lo) / dspan;
            dt = (dt <= 0.0) ? 0.0 : (dt >= 1.0 ? 1.0 : dt * dt * (3.0 - 2.0 * dt));
            charge *= dt;
        }
        return charge;
    }

    /// Whisker-spray drag (N), Savitsky-DeLorme-Datla 2007. Evaluated at the current
    /// (converged) running trim; gated to planing-form hulls and phased in by the planing
    /// blend weight. tan(alpha)=pi*tau/(2 tan beta); Dlam=cos(Theta)/(4 sin2a cos beta),
    /// Theta=2 alpha cos beta; Lww=b/(4 sin2a cos beta); Cf laminar/transitional on RNww;
    /// Rs = 1/2 rho V^2 b^2 Dlam Cf  (both sides; b = chine beam).
    double get_whisker_spray_drag(double speed) const {
        if (!cfg.whisker_spray_on || speed < 1e-2 || !hull) {
            return 0.0;
        }
        Hull& H = const_cast<Hull&>(*hull);
        if (!H.is_planing_form()) {
            return 0.0;
        }
        double const w = planing_weight(speed);
        if (w <= 0.0) {
            return 0.0;
        }
        // Warp-only scoping (default OFF = byte-identical). The Savitsky-DeLorme-Datla whisker
        // correlation OVER-predicts on clean constant-deadrise surfaces: the RANS decomposition
        // (Section 3) resolves the actual spray shear at only ~3% of the total, which the main
        // wetted-bottom friction already carries. The term's real value is the disturbed spray
        // over a WARPED bottom (whose sections cross varying deadrise), where it stands in for a
        // high-Fn warped-hull drag no closed form supplies. So apply it only to warped hulls
        // (same fore-aft deadrise-spread criterion the sectional closure uses); constant-deadrise
        // prisms take no whisker, which removes the over-growth the CvB cap was limiting.
        if (cfg.whisker_warp_only) {
            double const x_ap = H.get_ap(), x_fp = H.get_fp();
            double const bA = H.get_deadrise_at_x(x_ap + WARP_SAMPLE_AFT_FRAC * (x_fp - x_ap));
            double const bF = H.get_deadrise_at_x(x_ap + WARP_SAMPLE_FWD_FRAC * (x_fp - x_ap));
            if (!(bA > 0.0 && bF > 0.0 && std::abs(bF - bA) > cfg.sectional_warp_floor)) {
                return 0.0;
            }
        }
        // CvB cap (see cfg.wsk_cap_c; default 0 = byte-identical)
        double cap_scale = 1.0;
        if (cfg.wsk_cap_c > 0.0) {
            double const cvb_cap = speed / std::sqrt(env->get_gravity()
                                                     * std::max(1e-6, H.get_beam_chine()));
            double tc = (cvb_cap - cfg.wsk_cap_cv_lo)
                        / std::max(1e-6, cfg.wsk_cap_cv_hi - cfg.wsk_cap_cv_lo);
            tc = (tc <= 0.0) ? 0.0 : (tc >= 1.0 ? 1.0 : tc * tc * (3.0 - 2.0 * tc));
            cap_scale = 1.0 - cfg.wsk_cap_c * tc;
        }
        // Deadrise for the spray closed form. Default: aft-half static mean. physfix C2: on a warped
        // hull that mean misreports the deadrise in the wetted spray region -> use the
        // wetted-half-beam-weighted local mean over the actual wetted length. Guard on the fore-aft
        // spread (same stations as the sectional warp gate) so constant-deadrise hulls keep the exact
        // scalar path and stay byte-identical.
        double beta_in = H.get_deadrise();
        double g_turbwarp = 0.0;   // warp-keyed turbulent-Cf blend weight (see Config)
        if (cfg.whisker_spray_root || cfg.whisker_local_beta || cfg.whisker_spray_turb_warp) {
            double const x_ap = H.get_ap(), x_fp = H.get_fp();
            double const bA = H.get_deadrise_at_x(x_ap + WARP_SAMPLE_AFT_FRAC * (x_fp - x_ap));
            double const bF = H.get_deadrise_at_x(x_ap + WARP_SAMPLE_FWD_FRAC * (x_fp - x_ap));
            if (bA > 0.0 && bF > 0.0 && std::abs(bF - bA) > 1e-6) {   // warp gate: prisms byte-identical
                // Shared Fn_vol phase-in (see Config): below the window keep the aft-mean beta /
                // laminar line; blend to beta_F / turbulent over Fn_vol in [lo, hi]. Keeps both
                // warp closures off the over-predicted mid-Fn band. hi <= lo disables (gfn = 1).
                double gfn = 1.0;
                if (cfg.whisker_spray_root_fn_hi > cfg.whisker_spray_root_fn_lo) {
                    double const vol13w = std::pow(
                        std::max(hull->get_reference_displacement(), 1e-12), 1.0 / 3.0);
                    double const fnv = hull->get_fn(speed, vol13w);
                    double t = (fnv - cfg.whisker_spray_root_fn_lo)
                             / (cfg.whisker_spray_root_fn_hi - cfg.whisker_spray_root_fn_lo);
                    t = std::max(0.0, std::min(1.0, t));
                    gfn = t * t * (3.0 - 2.0 * t);   // smoothstep
                }
                if (cfg.whisker_spray_turb_warp && std::abs(bF - bA) > cfg.sectional_warp_floor) {
                    g_turbwarp = gfn;   // strongly-warped only (the sectional/crossflow set)
                }
                if (cfg.whisker_spray_root) {
                    // B1: forward (spray-root) section deadrise, phased in by gfn
                    beta_in = (1.0 - gfn) * H.get_deadrise() + gfn * bF;
                } else if (cfg.whisker_local_beta) {
                    double const wm = H.get_wetted_mean_deadrise();
                    if (wm > 0.0) beta_in = wm;
                }
            }
        }
        double const beta = std::max(radians(1.0), std::min(beta_in, radians(50.0)));
        double const tau = std::max(radians(0.25), -H.get_pitch());   // running trim
        double const b = std::max(H.get_beam_chine(), 1e-3);
        double const tan_alpha = M_PI * tau / (2.0 * std::tan(beta));
        double const alpha = std::atan(tan_alpha);
        double const sin2a = std::sin(2.0 * alpha);
        if (sin2a < 1e-4) {
            return 0.0;
        }
        double const cos_beta = std::cos(beta);
        double const Theta = 2.0 * alpha * cos_beta;            // spray dir rel. keel
        double const cosTheta = std::cos(Theta);
        if (cosTheta <= 0.0) {
            return 0.0;                                          // spray aft of chine, no drag
        }
        double const dlam = cosTheta / (4.0 * sin2a * cos_beta);
        double const Lww = b / (4.0 * sin2a * cos_beta);         // spray characteristic length
        double const nu = env->get_viscosity();
        double const RNww = std::max(1.0, speed * Lww / nu);
        double Cf;
        if (cfg.whisker_spray_turbulent) {
            Cf = 0.074 / std::pow(RNww, 0.2);                    // B0: forced fully-turbulent line (bounding)
        } else if (RNww < 1.5e6) {
            Cf = 1.328 / std::sqrt(RNww);                        // laminar/Blasius (model scale)
        } else {
            Cf = 0.074 / std::pow(RNww, 0.2) - 4800.0 / RNww;    // transitional/Schoenherr (full scale)
        }
        if (g_turbwarp > 0.0 && !cfg.whisker_spray_turbulent) {
            // Warp-keyed turbulent line (whisker_spray_turb_warp): the disturbed spray sheet on a
            // strongly-warped hull takes the turbulent Cf, phased in with the Fn_vol window.
            double const Cf_turb = 0.074 / std::pow(RNww, 0.2);
            Cf = (1.0 - g_turbwarp) * Cf + g_turbwarp * Cf_turb;
        }
        double const rho = env->get_density();
        double const Rs = 0.5 * rho * sq(speed) * sq(b) * dlam * std::max(0.0, Cf);
        return cap_scale * w * Rs;                               // phase in with planing
    }

    /// Toggle the Mercier-Savitsky pre-planing residuary term. Off by default.
    void set_mercier_savitsky(bool on) { cfg.mercier_savitsky_on = on; }

    /// Mercier-Savitsky (1973) pre-planing residuary resistance (N), gated to planing-form
    /// hulls and ramped over the band (0 below Fn_vol 0.9, full in [1,2], fading to 0 by
    /// 3.0). Evaluates Eq.(6) for R_T/Delta interpolated in Fn_vol, subtracts the reference
    /// (100,000-lb) Schoenherr friction to get the scale-independent residuary, x weight.
    double get_mercier_savitsky_residuary(double speed) const {
        if (!cfg.mercier_savitsky_on || speed < 1e-2 || !hull) {
            return 0.0;
        }
        Hull& H = const_cast<Hull&>(*hull);
        if (!H.is_planing_form()) {
            return 0.0;
        }
        // Form parameters are static design-pose descriptors -> read them from the ORIGINAL
        // hull (held at the design float), NOT the internal hull which is at the risen
        // running attitude when this is evaluated after equilibrium.
        Hull& OH = const_cast<Hull&>(*original_hull);
        double const g = env->get_gravity();
        double const rho = env->get_density();
        double const vol = std::max(1e-9, OH.get_reference_displacement());
        double const vol13 = std::pow(vol, 1.0 / 3.0);
        double const Fnv = H.get_fn(speed, vol13);
        // band ramp (cubic smoothstep edges): 0 below 0.9, 1 in [1,2], 0 above 3.0
        double g_ms;
        if (Fnv <= 0.9 || Fnv >= 3.0) {
            return 0.0;
        } else if (Fnv < 1.0) {
            double const t = (Fnv - 0.9) / 0.1; g_ms = t * t * (3.0 - 2.0 * t);
        } else if (Fnv <= 2.0) {
            g_ms = 1.0;
        } else {
            double const t = (3.0 - Fnv) / 1.0; g_ms = t * t * (3.0 - 2.0 * t);
        }
        // Read the geometric form descriptors at the DESIGN pose. The caller (warm-start
        // continuation) may leave original_hull at the PREVIOUS running attitude, so reading
        // these attitude-dependent quantities directly makes the residuary — and hence the
        // total resistance — depend on the start attitude (cold != warm at planing inception,
        // where the heave equilibrium is multi-valued). Reset to the design float, read, and
        // restore the caller's attitude so the term is seed-independent.
        double const oh_heave0 = OH.get_heave(), oh_pitch0 = OH.get_pitch();
        OH.set_heave(0.0); OH.set_pitch(0.0);
        double const L = std::max(1e-3, OH.get_length_wl());
        double const b = std::max(1e-3, OH.get_beam_chine());
        double const ie = std::min(35.0, std::max(2.0, OH.get_entrance_half_angle()));
        // Transom area: immersed section a hair forward of the AP (calc_transom_area reads
        // 0 on the coplanar transom face). A_X = true max section (the single-station
        // get_max_section_area can miss the transom-as-widest on planing hulls), so W<=1.
        double const fwd_dir = double(OH.get_fwd());
        double const At = std::max(0.0, OH.get_section_area_at(OH.get_ap() + fwd_dir * 0.01 * L));
        double const Ax = std::max(1e-9, std::max(OH.get_max_section_area(), At));
        OH.set_heave(oh_heave0); OH.set_pitch(oh_pitch0);
        // form parameters (design pose), clamped to the regression's range of applicability
        double const X = std::min(0.30, std::max(0.08, vol13 / L));
        double const U = std::sqrt(2.0 * ie);
        double const Z = std::min(1.2, std::max(0.30, vol / (b * b * b)));
        double const W = std::min(1.0, std::max(0.0, At / Ax));
        // R_T/Delta from Eq.(6), linearly interpolated between bracketing Fn_vol tenths
        double const fc = std::min(2.0, std::max(1.0, Fnv));
        int const i0 = std::min(9, std::max(0, int(std::floor((fc - 1.0) / 0.1))));
        double const frac = (fc - (1.0 + 0.1 * i0)) / 0.1;
        double const term[14] = {1.0, X, U, W, X*Z, X*U, X*W, Z*U, Z*W, W*W, X*W*W, Z*X*X, U*W*W, W*U*U};
        double rt0 = 0.0, rt1 = 0.0;
        for (int j = 0; j < 14; j++) { rt0 += MS_A[j][i0] * term[j]; rt1 += MS_A[j][i0 + 1] * term[j]; }
        double const RT_D = (1.0 - frac) * rt0 + frac * rt1;
        // subtract the reference (100,000-lb, C_A=0) friction to get residuary R_R/Delta.
        // Design-pose B/T (H is at the running attitude here): use the captured reference geometry,
        // consistent with the design-pose L/b/ie/At above (else B/T leaks the solver seed: cold!=warm).
        double const BT = OH.get_reference_beam_wl() / std::max(1e-4, OH.get_reference_draft());
        double const S_v23 = 2.262 * std::sqrt(L / vol13) * (1.0 + 0.046 * BT + 0.00287 * BT * BT);
        // Reference-craft Reynolds Re = V_ref*L_ref/nu for the 100,000-lb (salt water, 64 lb/ft^3)
        // reference craft at this Fn_vol: V_ref = Fnv*sqrt(g*vol13_ref), L_ref = (L/vol13)*vol13_ref,
        // vol13_ref = (100000/64)^(1/3), g = 32.2 ft/s^2, nu = 1.2817e-5 ft^2/s. Collapsing the
        // vol13_ref powers gives Re = Fnv*(L/vol13)*sqrt(32.2*100000/64)/nu. The prior form used
        // sqrt(L/vol13) and (100000/64)^(1/3) under the root -> ~27x too low -> CF_ref ~1.7x too high
        // -> residuary under-added (~1% of Delta at mid-Fnv). SOLVER_FIX_PLAN.md P3 (M&S 1973 / R-1667).
        double const Rn_ref = Fnv * (L / vol13)
            * std::sqrt(32.2 * 100000.0 / 64.0) / 1.2817e-5;
        double const CF_ref = 0.075 / sq(std::log10(std::max(1e3, Rn_ref)) - 2.0);  // ~Schoenherr
        double RR_D = RT_D - CF_ref * S_v23 * sq(Fnv) * 0.5;
        if (RR_D < 0.0) RR_D = 0.0;
        double const W_weight = vol * rho * g;
        return g_ms * RR_D * W_weight * (michell.is_demihull() + 1);
    }

    /// Pre-planing residuary hump resistance (N): a compact bump in R/Delta vs Fn_vol,
    /// active only for planing-form hulls inside the hump band. See member notes above.
    double get_preplaning_hump(double speed) const {
        if (speed < 1e-2 || !hull || cfg.preplaning_hump_coef <= 0.0) {
            return 0.0;
        }
        if (!const_cast<Hull&>(*hull).is_planing_form()) {
            return 0.0;
        }
        double const f = const_cast<Hull&>(*hull).get_fn(speed,
            std::pow(std::max(1e-12, const_cast<Hull&>(*hull).get_reference_displacement()), 1.0/3.0));
        if (f <= cfg.preplaning_fn_lo || f >= cfg.preplaning_fn_hi) {
            return 0.0;
        }
        double const span = cfg.preplaning_fn_hi - cfg.preplaning_fn_lo;
        double const bump = 4.0 * (f - cfg.preplaning_fn_lo) * (cfg.preplaning_fn_hi - f) / (span * span);
        double const W = const_cast<Hull&>(*hull).get_reference_displacement()
                         * env->get_density() * env->get_gravity();
        return cfg.preplaning_hump_coef * W * bump;
    }

    double get_friction_resistance(double speed, ResistanceMethod method) {

        if (!is_ok()) {
            return 0;
        }

        if (speed < 1e-2) {
            return 0;
        }

        // Savitsky has its own integral wetted-bottom friction; use it directly.
        if (method == ResistanceMethod::Savitsky) {
            return savitsky.calc(speed).friction * (michell.is_demihull() + 1);
        }

        auto Cf_roughness = get_roughness_coef(speed);
        auto Cf = (method == ResistanceMethod::Holtrop) ? get_friction_coef(speed) : michell.get_friction_coef(speed);
        double const Rf_disp = hull->get_force(speed, Cf + Cf_roughness) * (michell.is_demihull() + 1);

        if (method == ResistanceMethod::Auto) {
            double const w = planing_weight(speed);
            if (w <= 0.0) {
                return Rf_disp;
            }
            // Friction at the running (equilibrium) trim, consistent with the converged
            // attitude and the Morabito pressure — not the standalone predictive trim.
            // (Toggle: predictive-trim λ recovers Savitsky's residuary-calibrated friction.)
            double const Rf_plan = savitsky.calc(speed, cfg.planing_friction_running_trim).friction
                                   * (michell.is_demihull() + 1);
            if (w >= 1.0) {
                return Rf_plan;
            }
            return (1.0 - w) * Rf_disp + w * Rf_plan;
        }

        return Rf_disp;
    }

    double get_viscous_pressure_resistance(double speed, ResistanceMethod method) {

        if (!is_ok()) {
            return 0;
        }

        if (speed < 1e-2) {
            return 0;
        }

        auto Rp_local = 0.0;

        double const Rf = get_friction_resistance(speed, method);
        double const R_holtrop = holtrop.get_viscous_pressure_resistance(speed, Rf);
        double const R_transom = michell.get_transom_resistance(speed);

        if (method == ResistanceMethod::Morabito) {

            // Morabito bottom-pressure integral provides the dominant pressure drag;
            // the transom term is absorbed into the PT roll-off, so we omit R_transom.
            Rp_local = morabito.get_drag_lift_torque(speed).x;

        } else if (method == ResistanceMethod::Savitsky) {

            // Savitsky's pressure drag (W tan tau) is fully carried by get_wave_resistance
            // (calc().wave); the published method has no separate viscous-pressure term.
            // The old R_holtrop + R_transom lump was measured NOT small (savblend Stage -1,
            // 2026-07-03): Holtrop form drag from static hull coefficients is meaningless on
            // a planing prism and reached 56% of the Begovic mono standalone total (90.9 N
            // at v=8 vs warp3's 11.9 N) — the reason the standalone benchmark over-predicted
            // ~2x on mono while MaxSurf's Savitsky tracks EFD (fixes_ab/savblend/).
            Rp_local = 0.0;

        } else if (method == ResistanceMethod::Auto) {

            // Displacement-regime viscous pressure uses Holtrop's empirical form drag,
            // NOT the Noblesse local-pressure integral: the latter over-predicts form
            // drag 7-20x for full hulls (thin-ship slope^2 breakdown) unless a learned
            // deflection field corrects it. Michell wave + ITTC friction + Holtrop form
            // drag validates to ~7% on displacement hulls with no calibration; Morabito
            // takes over the pressure drag as the hull starts to plane. (The Unified
            // method below keeps Noblesse for the deflection-calibrated pipeline.)
            double const w = planing_weight(speed);
            if (w <= 0.0) {
                Rp_local = R_holtrop + R_transom;
            } else if (w >= 1.0) {
                Rp_local = cfg.pdyn_on ? get_pdyn_pressure_drag(speed)
                                       : morabito.get_drag_lift_torque(speed).x;
            } else {
                double const R_morabito = cfg.pdyn_on ? get_pdyn_pressure_drag(speed)
                                                      : morabito.get_drag_lift_torque(speed).x;
                // Transom term phases out as planing weight rises (Morabito handles the wake).
                Rp_local = (1.0 - w) * (R_holtrop + R_transom)
                           + w * R_morabito;
            }

        } else if (method == ResistanceMethod::Unified) {

            auto R_noblesse = localflow.get_drag_lift_torque_noblesse(speed).x;

            Rp_local = R_noblesse + R_transom;

        } else {

            Rp_local = R_holtrop + R_transom;

        }

        return Rp_local * (michell.is_demihull() + 1);
    }

    double get_wave_resistance(double speed, ResistanceMethod method) {

        if (!is_ok()) {
            return 0;
        }

        if (speed < 1e-2) {
            return 0;
        }

        switch (method) {
        case ResistanceMethod::Holtrop:
            return holtrop.get_wave_resistance(speed) * (michell.is_demihull() + 1);
        case ResistanceMethod::Integral:
            return michell.get_drag_lift_torque(speed).x;
        case ResistanceMethod::Savitsky:
            // Savitsky's closed-form lumps free-surface waves and lift-induced
            // pressure drag into a single mg·tan(τ) term, calibrated against
            // planing-craft measurements. Don't add anything else.
            return savitsky.calc(speed).wave * (michell.is_demihull() + 1);
        case ResistanceMethod::Morabito:
            // Free-surface wave-making (Michell) and Morabito's trim-induced bottom-
            // pressure drag are distinct mechanisms and not double-counted. Both remain
            // active across the semi-planing range, so keep the full Michell wave term.
            return michell.get_wave_resistance(speed);
        case ResistanceMethod::Auto: {
            // As Morabito, plus the pre-planing residuary hump (off by default). The
            // Mercier-Savitsky residuary is folded in at the get_total level as a
            // supplement (it overlaps the Morabito induced drag, so it cannot just be
            // added to the wave term — see get_total_resistance).
            double wave = michell.get_wave_resistance(speed) + get_preplaning_hump(speed);
            // Local free-surface wave envelope (default ON, planing/warp-gated): take the larger of the
            // thin-ship Michell wave and the LocalFlow near-field drag (two estimates of
            // the same wave-making) where the hull is semi-planing/planing. Closes the
            // warped-hull hump/high-Fn wave deficit Michell under-resolves; gated by the
            // planing weight so the displacement regime is byte-identical.
            if (cfg.localfs_wave_envelope && planing_weight(speed) > 1e-6) {
                wave = std::max(wave, get_localflow_drag(speed));
            }
            return wave;
        }
        default:
            return michell.get_wave_resistance(speed);
        }
    }

    double get_pressure_resistance(double speed, ResistanceMethod method) {

        if (!is_ok()) {
            return 0;
        }

        if (speed < 1e-2) {
            return 0;
        }

        return get_wave_resistance(speed, method) + get_viscous_pressure_resistance(speed, method);
    }

    /// Convenience entry point: solve the running equilibrium attitude first, then evaluate the
    /// total resistance at it. Mirrors the validation harnesses' (find_equilibrium; prepare;
    /// get_total_resistance) sequence in one call, so a caller cannot accidentally read resistance
    /// at the static design pose (get_total_resistance does NOT solve the attitude on its own).
    double get_total_resistance_equilibrium(double speed, ResistanceMethod method, std::size_t iters = 150) {
        find_equilibrium(speed, iters);
        prepare();
        return get_total_resistance(speed, method);
    }

    double get_total_resistance(double speed, ResistanceMethod method) {

        if (!is_ok()) {
            return 0;
        }

        if (speed < 1e-2) {
            return 0;
        }

        if (method == ResistanceMethod::Auto) {
            double const Rf = get_friction_resistance(speed, method);
            double const wave = get_wave_resistance(speed, method);              // Michell + hump
            double const vp = get_viscous_pressure_resistance(speed, method);    // Morabito induced (planing)
            // F8: sectional pressure-drag consistency (default OFF). Counted WITH vp below so the
            // residuary shortfall and the NSS total envelope re-gate around it (displaces, not stacks).
            double const sec_pd = get_sectional_pressure_drag(speed);
            double R = Rf + wave + vp + sec_pd + get_whisker_spray_drag(speed);
            // Mercier-Savitsky pre-planing RESIDUARY: add only the shortfall beyond the
            // wave + induced-pressure our model already carries (the MS residuary overlaps
            // the Morabito induced drag, so it must not be added on top — that is the F3
            // residuary deficit). Gated/ramped inside get_mercier_savitsky_residuary.
            // Residuary-closure precedence: Radojcic NSS (default) takes priority; if it is
            // disabled, fall back to Mercier-Savitsky; if both are off, no residuary is added.
            // (Previously this was a ternary that ignored mercier_savitsky_on whenever Radojcic
            //  was on, so set_mercier_savitsky(false) silently still ran a residuary term.)
            double ms = 0.0;
            if (cfg.radojcic_nss_on)          ms = get_radojcic_nss_residuary(speed);
            else if (cfg.mercier_savitsky_on) ms = get_mercier_savitsky_residuary(speed);
            if (ms > 0.0 && cfg.radojcic_residuary_beta_gate && cfg.radojcic_nss_on) {
                // The envelope's deadrise gate applied to the sibling residuary-excess term:
                // the NSS regression is deadrise-blind, so degate it on near-flat prisms
                // (same smoothstep([11,15] deg) on the static area-mean deadrise as below).
                double const beta_deg = const_cast<Hull&>(*original_hull).get_deadrise() * 180.0 / M_PI;
                double const tb = std::min(1.0, std::max(0.0, (beta_deg - 11.0) / 4.0));
                ms *= tb * tb * (3.0 - 2.0 * tb);
            }
            if (ms > 0.0) {
                R += std::max(0.0, ms - wave - vp - sec_pd);
            }
            // Radojcic NSS TOTAL-resistance envelope (default OFF): floor the planing total to the
            // regression's total RT. Carries the high-Fn drag the residuary hump-closures drop to zero.
            // planing-weighted so the displacement regime is untouched; max() so it only lifts
            // under-predictions (never lowers the warped over-trim, and inert where the model already
            // matches the regression, e.g. the prism series). Slender hulls return 0 -> inert.
            if (cfg.radojcic_total_envelope) {
                double const r_tot = get_radojcic_nss_total(speed);
                if (r_tot > 0.0) {
                    // High-Fn gate: the envelope only acts where the residuary hump-closures
                    // (MS / Radojcic-residuary) have tapered to zero and the component model loses
                    // the growing high-Fn drag (Fn_vol > ~2.5; see decomp). A cubic smoothstep over
                    // [2.0, 2.6] ramps it in there, so the mid-Fn hump (handled by the residuary, and
                    // where the prism series sit and already match) is left untouched -> do-no-harm.
                    double const vol13 = std::pow(std::max(1e-9,
                        const_cast<Hull&>(*original_hull).get_reference_displacement()), 1.0 / 3.0);
                    double const Fnv = const_cast<Hull&>(*hull).get_fn(speed, vol13);
                    double t = std::min(1.0, std::max(0.0, (Fnv - 2.0) / 0.6));
                    double const g_fn = t * t * (3.0 - 2.0 * t);
                    // Deadrise gate: the Radojcic NSS regression is deadrise-blind (no beta input), so its
                    // total over-predicts very-low-deadrise planing prisms (beta ~ 10 deg). Restrict the
                    // envelope to beta > ~13 deg (its reliable range) so it floors the moderate-deadrise
                    // hulls that need it (Begovic ~17 deg, Taunton 22.5 deg) without over-correcting the
                    // near-flat prisms. Smoothstep over [11, 15] deg on the static deadrise.
                    double const beta_deg = const_cast<Hull&>(*original_hull).get_deadrise() * 180.0 / M_PI;
                    double tb = std::min(1.0, std::max(0.0, (beta_deg - 11.0) / 4.0));
                    double const g_beta = tb * tb * (3.0 - 2.0 * tb);
                    R += g_fn * g_beta * planing_weight(speed) * std::max(0.0, r_tot - R);
                }
            }
            // Above-water windage (calm-air drag), off by default. Method-blind: added here
            // for Auto and below for the component methods.
            if (cfg.air_resistance_on) {
                R += get_air_resistance(speed);
            }
            return R;
        }
        double R = get_friction_resistance(speed, method) + get_pressure_resistance(speed, method);
        if (cfg.air_resistance_on) {
            R += get_air_resistance(speed);
        }
        return R;
    }

//    std::pair<double,double> get_lift_and_torque(double speed) {

//        auto const ret = michell.get_drag_lift_torque(speed);
//        std::cout << "L: " << ret.first << " N, Q: " << ret.second << std::endl;
//        return ret;

//    }

    BiLinearInterpolator get_wave_pressure_field(double speed) {

        std::vector<double> pressure_field;
        michell.get_drag_lift_torque(speed, &pressure_field);
        return BiLinearInterpolator(michell.get_grid()->get_xs(), michell.get_grid()->get_zs(), pressure_field);

    }

    /// Centerline wave elevation E(x, 0, 0) in metres, summing the wave (Michell/Havelock)
    /// and local-flow (Noblesse) contributions per paper Eqs. 5-6.
    std::vector<double> get_centerline_wave_elevation(double speed, const std::vector<double>& xs_eval) {

        std::vector<double> result(xs_eval.size(), 0.0);

        if (!is_ok() || speed < 1e-2) {
            return result;
        }

        auto E_wave = michell.get_wave_elevation_centerline(speed, xs_eval);
        auto E_local = localflow.get_local_elevation_centerline(speed, xs_eval);

        for (std::size_t i = 0; i < xs_eval.size(); i++) {
            result[i] = E_wave[i] + E_local[i];
        }
        return result;
    }

    /// Wave-only centerline elevation (no local-flow correction). Useful for comparison
    /// against thin-ship Michell-only references like Noblesse Fig. 16's wave component.
    std::vector<double> get_centerline_wave_elevation_michell(double speed, const std::vector<double>& xs_eval) {
        if (!is_ok() || speed < 1e-2) {
            return std::vector<double>(xs_eval.size(), 0.0);
        }
        return michell.get_wave_elevation_centerline(speed, xs_eval);
    }

    /// Local-only centerline elevation correction.
    std::vector<double> get_centerline_wave_elevation_local(double speed, const std::vector<double>& xs_eval) {
        if (!is_ok() || speed < 1e-2) {
            return std::vector<double>(xs_eval.size(), 0.0);
        }
        return localflow.get_local_elevation_centerline(speed, xs_eval);
    }

    /// 2-D far-field wave elevation Z(x, y) over the xs × ys grid (Tuck-Scullen-Lazauskas
    /// Fourier-Kochin form). Returns a flattened (ny × nx) array with x varying fastest.
    /// Z is inaccurate within ~1 ship length of the hull (no near-field correction).
    std::vector<double> get_wave_field(double speed,
                                       const std::vector<double>& xs,
                                       const std::vector<double>& ys,
                                       double viscosity = 5e-3,
                                       int n_theta = 401) {

        if (!is_ok() || speed < 1e-2) {
            return std::vector<double>(xs.size() * ys.size(), 0.0);
        }
        return michell.get_farfield_wave_field(speed, xs, ys, viscosity, n_theta);
    }

    /// Sample the disturbance velocity on a propeller disc and return non-dimensional
    /// axial/radial/tangential components plus area-weighted statistics. Combines the
    /// ThinShip Havelock kernel (Rankine + free-surface image + far-field wave term)
    /// with a near-field axial-deficit correction from the LocalFlow Noblesse Green
    /// function. include_localflow_axial=false skips the latter, useful to isolate the
    /// thin-ship potential prediction.
    PropellerWakeField get_propeller_wake_field(double speed,
                                                glm::dvec3 center,
                                                double inner_diameter,
                                                double outer_diameter,
                                                int n_radial = 12,
                                                int n_azimuthal = 36,
                                                bool include_localflow_axial = true) {

        PropellerWakeField wf;
        wf.center = center;
        wf.R_inner = 0.5 * inner_diameter;
        wf.R_outer = 0.5 * outer_diameter;
        wf.n_radial = n_radial;
        wf.n_azimuthal = n_azimuthal;
        wf.ship_speed = speed;

        if (!is_ok() || speed < 1e-2 || n_radial < 2 || n_azimuthal < 3
            || wf.R_outer <= wf.R_inner) {
            return wf;
        }

        wf.r.resize(n_radial);
        wf.theta.resize(n_azimuthal);
        std::size_t const N = static_cast<std::size_t>(n_radial) * n_azimuthal;
        wf.axial.assign(N, 1.0);
        wf.radial.assign(N, 0.0);
        wf.tangential.assign(N, 0.0);

        for (int i = 0; i < n_radial; i++) {
            wf.r[i] = wf.R_inner + (wf.R_outer - wf.R_inner) * (double(i) / double(n_radial - 1));
        }
        for (int j = 0; j < n_azimuthal; j++) {
            wf.theta[j] = (2.0 * M_PI) * (double(j) / double(n_azimuthal));
        }

        int const fwd = hull->get_fwd();

        // Pre-tabulate the spectrum once on the calling thread so the cache is warm
        // before the OpenMP region samples it.
        michell.get_velocity(speed, center.x, center.y, center.z);

        #pragma omp parallel for collapse(2)
        for (int i = 0; i < n_radial; i++) {
            for (int j = 0; j < n_azimuthal; j++) {

                double const r = wf.r[i];
                double const th = wf.theta[j];
                double const sT = std::sin(th);
                double const cT = std::cos(th);

                double const X = center.x;                  // disc plane is y-z at center.x
                double const Y = center.y + r * sT;         // theta = 0 is +z (top), increasing CW looking forward
                double const Z = center.z + r * cT;

                glm::dvec3 u = michell.get_velocity(speed, X, Y, Z);

                if (include_localflow_axial) {
                    glm::dvec3 const u_loc = localflow.get_velocity_local(speed, X, Y, Z);
                    u.x += u_loc.x;
                    u.y += u_loc.y;
                    u.z += u_loc.z;
                }

                // axial inflow into the disc, normalised by ship speed.
                // freestream in ship frame is (-speed*fwd, 0, 0); axial = -freestream_x = speed*fwd.
                // total = freestream + disturbance, so axial_inflow = speed*fwd - u.x*fwd
                //                                                   = (speed - u.x*fwd)
                // Normalised: V_A/V = 1 - u.x*fwd/speed.
                double Va = 1.0 - (u.x * fwd) / speed;

                if (cfg.transom_hr_on) {
                    Va += cfg.transom_hr_amp * michell.get_transom_head_release(
                              speed, X, Y, Z, cfg.transom_hr_len);
                }

                // In-plane components: outward radial unit vector r_hat = (0, sT, cT)
                //                       tangential unit (CW from astern about +x_fwd):
                //                       t_hat = fwd * (0, cT, -sT)
                double const Vr = (u.y * sT + u.z * cT) / speed;
                double const Vt = fwd * (u.y * cT - u.z * sT) / speed;

                std::size_t const idx = static_cast<std::size_t>(i) * n_azimuthal + j;
                wf.axial[idx] = Va;
                wf.radial[idx] = Vr;
                wf.tangential[idx] = Vt;
            }
        }

        // Stats (area-weighted on dA = r dr dtheta, trapezoidal in r, uniform in theta).
        double sum_a = 0.0, sum_a2 = 0.0;
        double sum_r2 = 0.0, sum_t2 = 0.0;
        double sum_w = 0.0;
        double amin = 1e9, amax = -1e9;
        double const dtheta = (2.0 * M_PI) / double(n_azimuthal);

        for (int i = 0; i < n_radial; i++) {
            double const r = wf.r[i];
            // trapezoidal weight in r
            double dr;
            if (n_radial == 1) {
                dr = 0.0;
            } else if (i == 0 || i == n_radial - 1) {
                dr = 0.5 * (wf.R_outer - wf.R_inner) / double(n_radial - 1);
            } else {
                dr = (wf.R_outer - wf.R_inner) / double(n_radial - 1);
            }
            double const w_r = r * dr;

            for (int j = 0; j < n_azimuthal; j++) {
                std::size_t const idx = static_cast<std::size_t>(i) * n_azimuthal + j;
                double const w = w_r * dtheta;
                double const a = wf.axial[idx];
                sum_a += a * w;
                sum_a2 += sq(a) * w;
                sum_r2 += sq(wf.radial[idx]) * w;
                sum_t2 += sq(wf.tangential[idx]) * w;
                sum_w += w;
                amin = std::min(amin, a);
                amax = std::max(amax, a);
            }
        }

        wf.disc_area = M_PI * (sq(wf.R_outer) - sq(wf.R_inner));
        if (sum_w > 0.0) {
            wf.mean_axial = sum_a / sum_w;
            double const var = std::max(0.0, sum_a2 / sum_w - sq(wf.mean_axial));
            wf.axial_std = std::sqrt(var);
            wf.radial_rms = std::sqrt(std::max(0.0, sum_r2 / sum_w));
            wf.swirl_rms = std::sqrt(std::max(0.0, sum_t2 / sum_w));
        }
        wf.wake_fraction = 1.0 - wf.mean_axial;
        wf.axial_min = amin;
        wf.axial_max = amax;
        if (wf.mean_axial > 1e-6) {
            wf.non_uniformity = (amax - amin) / wf.mean_axial;
        }
        wf.volumetric_flow = wf.mean_axial * speed * wf.disc_area;

        return wf;
    }

    BiLinearInterpolator get_local_pressure_field(double speed) {

        std::vector<double> pressure_field;
        localflow.get_drag_lift_torque_noblesse(speed, &pressure_field);
        return BiLinearInterpolator(localflow.get_grid()->get_xs(), localflow.get_grid()->get_zs(), pressure_field);

    }

    /// Morabito bottom-pressure field over the prismatic-hull (X, Y) plane.
    /// xs ∈ [0, λ], ys ∈ [-0.5, 0.5] in beam units; field is Pa.
    std::vector<double> get_morabito_pressure_field(double speed) {
        std::vector<double> pressure_field;
        morabito.get_drag_lift_torque(speed, &pressure_field);
        return pressure_field;
    }

    /// Diagnostic: the Morabito (drag, lift, torque-about-CG) at the current attitude AND the
    /// exact planform lambda the kernel used to produce them — [drag, lift, torque, lambda].
    /// Lets a port be validated at the integrated-force level (planing-kernel study).
    std::vector<double> get_morabito_force(double speed) {
        double lam = 0.0, tau = 0.0;
        glm::dvec3 f = morabito.get_drag_lift_torque(speed, nullptr, &lam, &tau);
        return {f.x, f.y, f.z, lam, tau};
    }

    /// Stage-B CoP diagnostic: split the Morabito planing lift+moment (about CG) at the
    /// CURRENT attitude into dynamic (Pdyn) and hydrostatic (Pstat) parts.
    /// Returns {lift_dyn, torque_dyn, lift_stat, torque_stat, torque_total} (N, N·m).
    std::vector<double> get_morabito_split(double speed) {
        glm::dvec4 split(0.0);
        glm::dvec3 f = morabito.get_drag_lift_torque(speed, nullptr, nullptr, nullptr, &split);
        return {split.x, split.y, split.z, split.w, f.z};
    }

    // ---- find_equilibrium tuning constants (fixed solver numerics, not user-facing levers) ----
    static constexpr double WARP_SAMPLE_AFT_FRAC = 0.15;   // aft deadrise sample, fraction of Lwl from AP
    static constexpr double WARP_SAMPLE_FWD_FRAC = 0.60;   // fwd deadrise sample, fraction of Lwl from AP
    static constexpr double SOLVER_EPS           = 1e-5;   // generic small epsilon (Newton denominator guard, Awp floor)
    static constexpr double TOL_FORCE_FRAC       = 1e-3;   // lift/moment convergence tolerance, fraction of W (×Lref for moment)
    static constexpr double HEAVE_LIM_FRAC       = 0.10;   // heave search bound, fraction of Lref
    static constexpr std::size_t EVAL_BUDGET_MIN      = 60; // floor on total residual evaluations
    static constexpr std::size_t EVAL_BUDGET_PER_ITER = 8;  // residual evaluations granted per requested iteration
    static constexpr int    HEAVE_NEWTON_ITERS   = 8;      // inner heave Newton/bisection iterations
    static constexpr double AWP_FLOOR_FRAC       = 0.05;   // waterplane-area floor in heave Newton, fraction of design Awp
    static constexpr double HEAVE_STEP_TOL_FRAC  = 1e-6;   // tiny-step convergence tolerance, fraction of Lref
    static constexpr double BALANCED_TOL_MULT    = 10.0;   // "balanced off the rails" threshold, multiple of tol_lift
    static constexpr double TRIM_WINDOW_LO_DEG   = -1.0;   // trim search window lower bound (deg)
    static constexpr double TRIM_WINDOW_HI_DEG   = 8.0;    // trim search window upper bound (deg)
    static constexpr int    TRIM_SCAN_N          = 37;     // trim scan samples (0.25-deg steps over the window)
    static constexpr double TRIM_REGRESSION_HALF_WIN_DEG = 1.5;  // local moment-regression half-window (deg)
    static constexpr double OVERIMMERSION_LIFT_MULT      = 10.0; // over-immersion sentinel lift, multiple of W
    static constexpr double ITP_TRIM_TOL_DEG     = 0.01;   // ITP trim-root tolerance (deg)
    static constexpr double ITP_K1               = 0.2;    // ITP truncation coefficient k1
    static constexpr double ITP_K2               = 2.0;    // ITP truncation exponent k2
    static constexpr int    ITP_N0               = 1;      // ITP slack iterations n0

    // Mutable working state threaded through the nested-1-D equilibrium solve (find_equilibrium ->
    // solve_trim_scan -> solve_heave_at_pitch -> eval_residual). Replaces the captured-by-reference
    // lambda state of the previous monolithic find_equilibrium, making the data flow explicit.
    struct EqContext {
        double speed = 0.0;
        double g = 0.0, rho = 0.0;          // gravity, water density
        double Lref = 0.0, Awp0 = 0.0;      // design-pose reference length / waterplane area
        double tol_lift = 0.0, heave_lim = 0.0;
        double tau_lo = 0.0, tau_hi = 0.0;  // trim search window (rad); trim = -pitch
        std::size_t evals = 0, eval_budget = 0;
        bool use_sec = false;               // warp gate: sectional water-entry closure active
        // Sectional-delta two-pass side-channel (inert unless sectional-delta mode is on):
        bool eval_delta_mode = false;
        double locked_t = -1.0;             // -1 => capture the blend fraction t dynamically
        double captured_t = 0.0;            // last dynamically-captured t
        bool sec_engaged = false;           // any eval engaged the sectional kernel (=> trim floor at 0)
        double h_guess = 0.0;               // warm-start heave threaded across the scan + ITP refinement
    };

    // Inner residual on the ATTITUDE band: sets the pose via residual_core, threading the eval counter
    // and the sectional-delta side-channel. (Was the eval_at lambda inside find_equilibrium.)
    void eval_residual(double heave, double pitch, EqContext& ctx, double& lift, double& torque) {
        SectionalMode const mode = ctx.eval_delta_mode
            ? SectionalMode::Delta
            : (cfg.sectional_delta_mode ? SectionalMode::None : SectionalMode::Replacement);
        DeltaState ds;
        ds.locked_t = ctx.locked_t;       // -1 => let the core capture t dynamically
        ds.captured_t = ctx.captured_t;   // preserved unless the core captures a new t
        PoseForces const pf = residual_core(heave, pitch, ctx.speed, BlendBand::Attitude, ctx.use_sec, mode, &ds);
        ++ctx.evals;   // one residual evaluation = one expensive wave-integral solve
        ctx.captured_t = ds.captured_t;
        ctx.sec_engaged = ctx.sec_engaged || ds.engaged;
        lift = pf.total_lift;
        torque = pf.total_torque;
    }

    // Inner solve: heave that zeroes vertical force at a fixed pitch. Warm-started from h_guess; returns
    // the balanced heave and the residual pitch moment there. Vertical force DECREASES with heave (raising
    // the hull reduces immersion), so the fully-immersed rail carries the most lift (f_lo > 0) and the
    // most-raised rail the least (f_hi < 0); safeguarded Newton (waterplane stiffness) + bisection.
    // (Was the solve_heave lambda.)
    double solve_heave_at_pitch(double pitch, EqContext& ctx, double h_guess, double& torque_out, bool& balanced) {
        double lift, torque;
        double lo = -ctx.heave_lim, hi = ctx.heave_lim;
        eval_residual(lo, pitch, ctx, lift, torque); double const f_lo = lift;
        eval_residual(hi, pitch, ctx, lift, torque); double const f_hi = lift;
        balanced = false;
        // If the force cannot be balanced inside the bound, settle on the rail (not a balanced point):
        // f_lo<=0 => too heavy even fully immersed; f_hi>=0 => still lifting fully raised (planing beyond bound).
        if (f_lo <= 0.0) { eval_residual(lo, pitch, ctx, lift, torque); torque_out = torque; return lo; }
        if (f_hi >= 0.0) { eval_residual(hi, pitch, ctx, lift, torque); torque_out = torque; return hi; }
        double h = glm::clamp(h_guess, lo, hi);
        for (int k = 0; k < HEAVE_NEWTON_ITERS && ctx.evals < ctx.eval_budget; k++) {
            eval_residual(h, pitch, ctx, lift, torque);
            if (std::abs(lift) < ctx.tol_lift) break;
            if (lift > 0.0) lo = h; else hi = h;      // excess lift => must rise (higher heave)
            double const Awp = std::max(hull->get_waterplane_area(), AWP_FLOOR_FRAC * ctx.Awp0);
            double hn = h + lift / (ctx.rho * ctx.g * Awp + SOLVER_EPS);   // Newton: f' = -rho g A_wp
            if (!(hn > lo && hn < hi)) hn = 0.5 * (lo + hi);   // bisection safeguard
            bool const tiny = std::abs(hn - h) < HEAVE_STEP_TOL_FRAC * ctx.Lref;
            h = hn;
            if (tiny) { eval_residual(h, pitch, ctx, lift, torque); break; }
        }
        torque_out = torque;
        balanced = std::abs(lift) < BALANCED_TOL_MULT * ctx.tol_lift;   // converged off the rails
        return h;
    }

    // Outer solve: trim whose heave-balanced pitch moment is zero. The heave-balanced moment(trim) curve
    // is restoring (slope < 0) but NOISY at ~1deg scale, so a single noisy crossing is seed-dependent.
    // Scan the window, then take the LOCAL linear-regression root about the STEEPEST-restoring up-trim
    // crossing (the dynamically stable branch, the same one at neighbouring speeds -> smooth attitude);
    // optionally refine to the true root with a bracketed ITP step. (Was the run_scan lambda.)
    double solve_trim_scan(EqContext& ctx) {
        int const Nscan = TRIM_SCAN_N;   // 0.25-deg steps over the window; fine enough that the local
                                         // moment-crossing regression has stable membership (less jitter)
        std::vector<double> scan_tau, scan_M;  // heave-balanced (trim, moment) samples
        double best_tau = 0.0, best_abs = std::numeric_limits<double>::max();
        for (int i = 0; i < Nscan && ctx.evals < ctx.eval_budget; i++) {
            double const tau = ctx.tau_lo + (ctx.tau_hi - ctx.tau_lo) * i / (Nscan - 1);
            double M; bool balanced;
            ctx.h_guess = solve_heave_at_pitch(-tau, ctx, ctx.h_guess, M, balanced);
            if (!balanced) continue;           // skip railed (non-balanced) points
            scan_tau.push_back(tau); scan_M.push_back(M);
            if (std::abs(M) < best_abs) { best_abs = std::abs(M); best_tau = tau; }
        }

        // Among every up-trim crossing (M:+ -> <=0), fit a local restoring regression and keep the one
        // with the STEEPEST (most negative) slope: the dynamically stable equilibrium, stable across speeds.
        double tau_star;
        std::size_t best_c = 0;                  // selected crossing's upper bracket index
        double const win = radians(TRIM_REGRESSION_HALF_WIN_DEG);   // local half-window about a crossing
        double best_slope = 0.0; bool found = false;
        for (std::size_t c = 1; c < scan_tau.size(); c++) {
            if (!(scan_M[c - 1] > 0.0 && scan_M[c] <= 0.0)) continue;
            double const tau_cross = 0.5 * (scan_tau[c - 1] + scan_tau[c]);
            double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0; int n = 0;
            for (std::size_t k = 0; k < scan_tau.size(); k++) {
                if (std::abs(scan_tau[k] - tau_cross) > win) continue;
                Sx += scan_tau[k]; Sy += scan_M[k];
                Sxx += scan_tau[k] * scan_tau[k]; Sxy += scan_tau[k] * scan_M[k]; ++n;
            }
            double const denom = n * Sxx - Sx * Sx;
            double const slope = (n >= 2 && std::abs(denom) > 1e-12) ? (n * Sxy - Sx * Sy) / denom : 0.0;
            double const root = (n >= 2 && slope < -1e-9)
                ? glm::clamp(-((Sy - slope * Sx) / n) / slope, ctx.tau_lo, ctx.tau_hi) : tau_cross;
            if (!found || slope < best_slope) { best_slope = slope; tau_star = root; best_c = c; found = true; }
        }
        if (!found) tau_star = best_tau;          // no sign change: trim with smallest |moment|

        // ITP refinement (use_itp_root, default ON): replace the regression estimate with the TRUE root
        // inside the selected crossing bracket, via the shared tools.h itp_root. g = -M is increasing through
        // the up-trim crossing; the closure threads the warm-start heave and the eval budget. The heave-
        // balanced moment is evaluated exactly as the scan does it. Toggle off => reverts to the regression.
        if (cfg.use_itp_root && found) {
            tau_star = itp_root(scan_tau[best_c - 1], scan_tau[best_c],
                                -scan_M[best_c - 1], -scan_M[best_c],
                                radians(ITP_TRIM_TOL_DEG), ITP_K1, ITP_K2, ITP_N0,
                [&](double x) {                       // g(x) = -M at trim x, threading h_guess + evals
                    double M; bool bal;
                    ctx.h_guess = solve_heave_at_pitch(-x, ctx, ctx.h_guess, M, bal);
                    return -M;
                },
                [&]() { return ctx.evals < ctx.eval_budget; });
        }

        // A planing hull does not run bow-down; floor the moment-balance trim at 0 once the sectional
        // closure is active (its heave-balanced moment goes noisy at the very-high-Fn / very-low-trim edge).
        if (ctx.use_sec || ctx.sec_engaged) tau_star = std::max(tau_star, 0.0);
        return tau_star;
    }

    /// Solve the running attitude (heave, pitch) that balances vertical force and
    /// pitch moment at the given speed, via a NESTED 1-D root find that decouples the
    /// two DOFs:
    ///   * inner — for a fixed trim, find the heave that zeroes the vertical force.
    ///     Vertical force rises monotonically with submergence, so the search interval
    ///     brackets the root; safeguarded Newton (waterplane stiffness) + bisection.
    ///   * outer — find the trim whose heave-balanced pitch moment is zero, via a scan +
    ///     linear-regression root over a physical trim window (the moment curve is
    ///     restoring but noisy at ~1° scale, so regression is used instead of chasing an
    ///     exact zero).
    /// Because the heave is always balanced before the moment is read and the trim is
    /// estimated by regression over a fixed scan, the converged attitude is the SAME
    /// solution regardless of the seed (cold or warm start) — unlike the previous coupled
    /// single-step scheme whose pitch secant was contaminated by the simultaneous heave
    /// move, so it oscillated and stuck against the ±7° trim clamp seed-dependently.
    glm::dvec3 find_equilibrium(double speed, std::size_t max_iterations) {
        double const g = env->get_gravity();
        double const rho = env->get_density();

        // Fixed reference (design) mass — independent of the original hull's current attitude, so the
        // caller can warm-start equilibrium from a planing pose.
        double const mass = original_hull->get_reference_mass();
        double const W = g * mass;

        // Reference geometry at the DESIGN pose (heave=0, pitch=0). These scale the heave bound and the
        // tolerances, so they MUST be seed-independent (reading them at a warm-start attitude reintroduced
        // seed-dependence). WARP-GATE the sectional closure: active only on warped hulls (fore-aft deadrise
        // spread above the floor), where Wagner-3/4 fails; on prisms/mild-warp it stays on Morabito.
        hull->set_heave(0.0); hull->set_pitch(0.0); prepare();
        // Fore-aft deadrise spread: the warp gate, shared by the sectional closure (use_sec) and the
        // Savitsky trim-weight ramp below. Sampled at fixed Lwl fractions => a stable hull-shape property.
        double const x_ap = hull->get_ap(), x_fp = hull->get_fp();
        double const bA = hull->get_deadrise_at_x(x_ap + WARP_SAMPLE_AFT_FRAC * (x_fp - x_ap));   // aft deadrise
        double const bF = hull->get_deadrise_at_x(x_ap + WARP_SAMPLE_FWD_FRAC * (x_fp - x_ap));   // fwd deadrise
        double const warp_spread = (bF > 0.0 && bA > 0.0) ? (bF - bA) : 0.0;
        bool use_sec = false;
        if (cfg.use_sectional_attitude) {
            sectional.capture_static();                            // static chine beam/length/waterline
            use_sec = warp_spread > cfg.sectional_warp_floor;
        }

        EqContext ctx;
        ctx.speed = speed;
        ctx.g = g; ctx.rho = rho;
        ctx.Lref = std::max(1e-3, hull->get_length_wl());
        ctx.Awp0 = std::max(SOLVER_EPS, hull->get_waterplane_area());
        ctx.tol_lift = TOL_FORCE_FRAC * W;
        ctx.heave_lim = HEAVE_LIM_FRAC * ctx.Lref;                 // physical heave bound (m)
        ctx.tau_lo = radians(TRIM_WINDOW_LO_DEG);                  // trim search window (rad); trim = -pitch
        ctx.tau_hi = radians(TRIM_WINDOW_HI_DEG);
        // Each residual evaluation is one expensive wave-integral solve, so cap the total cost.
        ctx.eval_budget = std::max<std::size_t>(EVAL_BUDGET_MIN, EVAL_BUDGET_PER_ITER * max_iterations);
        ctx.use_sec = use_sec;

        // Planing trim target for the Savitsky blend below: weighted by the attitude-band planing fraction.
        double const w_plan_eq = planing_weight_attitude(speed);

        bool const active_delta_mode = cfg.use_sectional_attitude && cfg.sectional_delta_mode && use_sec;
        ctx.eval_delta_mode = active_delta_mode;

        double tau_star;
        if (active_delta_mode) {
            // Option B: a first pass finds the baseline Morabito equilibrium; evaluate the deflection
            // there to lock the blend fraction, then run a blended second pass.
            ctx.eval_delta_mode = false;
            double const tau_base = solve_trim_scan(ctx);
            ctx.eval_delta_mode = true;
            ctx.locked_t = -1.0;                          // allow dynamic capture of t at the baseline
            double dummy_M; bool dummy_bal;
            solve_heave_at_pitch(-tau_base, ctx, 0.0, dummy_M, dummy_bal);
            ctx.locked_t = ctx.captured_t;                // lock t for the second pass
            tau_star = solve_trim_scan(ctx);
        } else {
            tau_star = solve_trim_scan(ctx);
        }

        // Blend the bottom-pressure moment-balance trim toward Savitsky's predictive trim as the regime
        // tends to planing (w_plan->1): at full planing the attitude is set by Savitsky's coupled lambda-tau
        // relation (correct falling trim). Heave is re-balanced at the blended trim afterwards.
        // Effective Savitsky trim weight: the flat cfg.savitsky_trim_weight, or — when the warp gate is
        // enabled (savitsky_trim_warp_w > 0) — ramped 0 -> warp_w by the fore-aft deadrise spread, so
        // prismatic hulls stay at 0 and only warped hulls get the partial blend. (Was a Python harness gate.)
        double trim_weight = cfg.savitsky_trim_weight;
        if (cfg.savitsky_trim_warp_w > 0.0) {
            double const span = std::max(cfg.savitsky_trim_warp_hi - cfg.savitsky_trim_warp_lo, 1e-6);
            double t = (warp_spread - cfg.savitsky_trim_warp_lo) / span;
            t = (t <= 0.0) ? 0.0 : (t >= 1.0 ? 1.0 : t * t * (3.0 - 2.0 * t));   // cubic smoothstep
            trim_weight = t * cfg.savitsky_trim_warp_w;
        }
        double const wt = w_plan_eq * trim_weight;
        if (wt > 0.0) {
            // Balance heave at the moment-balance trim, then read the effective planing deadrise there (the
            // wetted-width-weighted section-local deadrise over the physically-wetted bottom) — well below the
            // static area mean for a warped hull, which de-biases the otherwise over-predicted Savitsky trim.
            double m0; bool b0;
            solve_heave_at_pitch(-tau_star, ctx, ctx.h_guess, m0, b0);
            // s12 active: predict_trim self-samples the effective deadrise at the mean-wetted-
            // length station (and returns keel pitch); the wetted-mean override would double-
            // adjust the deadrise, so skip it.
            double const beta_eff = (cfg.savitsky2012_mode > 0)
                ? -1.0 : hull->get_wetted_mean_deadrise();
            // predict_trim mixes the design CG with body-frame geometry, so it MUST be evaluated at the design
            // pose (its LCG is corrupted at a bow-up attitude); restore it first, keeping only beta_eff.
            hull->set_heave(0.0); hull->set_pitch(0.0); prepare();
            double const tau_savitsky = (beta_eff > 0.0)
                ? radians(savitsky.predict_trim(speed, beta_eff))
                : radians(savitsky.predict_trim(speed));
            // predict_trim returns its bisection bracket end (~15 deg) when the CoP relation has no physical
            // root (e.g. a mis-set CG); blending toward that bogus value rails the trim. Only blend toward a
            // physically valid prediction; otherwise keep the scan-bounded moment-balance trim.
            if (tau_savitsky > 0.0 && tau_savitsky < radians(10.0)) {
                tau_star = (1.0 - wt) * tau_star + wt * tau_savitsky;
            }
        }
        // Final guard: the equilibrium trim cannot lie outside the physical search window.
        tau_star = glm::clamp(tau_star, ctx.tau_lo, ctx.tau_hi);

        // Drag-couple outer fixed-point: the couple-free trim above ignores the moment of the
        // horizontal resistance about the CG. Freeze the full ship R at that pose, then re-solve the
        // trim ONCE with the R*a moment active in residual_core. One pass suffices: the couple shifts
        // the trim ~0.5-1deg, across which the friction-dominated R changes <2%, so a second pass
        // would move the trim < ITP_TRIM_TOL_DEG. R enters only the MOMENT, never total_lift, so the
        // final heave-balance below is couple-independent; reset drag_couple_R_ first so a later bare
        // get_total_resistance at this pose is unaffected. get_total_resistance does NOT recurse into
        // find_equilibrium and mutates no pose.
        if (cfg.drag_couple_on && drag_couple_R_ == 0.0) {
            double gate = 1.0;
            if (cfg.drag_couple_warp_floor > 0.0) {   // optional warp gate (do-no-harm on prisms)
                double const span = std::max(cfg.drag_couple_warp_hi - cfg.drag_couple_warp_floor, 1e-6);
                double t = (warp_spread - cfg.drag_couple_warp_floor) / span;
                gate = (t <= 0.0) ? 0.0 : (t >= 1.0 ? 1.0 : t * t * (3.0 - 2.0 * t));   // cubic smoothstep
            }
            if (gate > 0.0) {
                double m0; bool b0;
                ctx.h_guess = solve_heave_at_pitch(-tau_star, ctx, ctx.h_guess, m0, b0);  // pin couple-free pose
                prepare();
                drag_couple_R_ = gate * get_total_resistance(speed, ResistanceMethod::Auto);
                ctx.evals = 0;                          // refresh the eval budget for the second trim solve
                tau_star = glm::clamp(solve_trim_scan(ctx), ctx.tau_lo, ctx.tau_hi);
                drag_couple_R_ = 0.0;                   // reset before the final heave => left pose couple-free
            }
        }

        // Leave the hull at the converged equilibrium attitude.
        double tq; bool bal;
        solve_heave_at_pitch(-tau_star, ctx, ctx.h_guess, tq, bal);

        return {hull->get_heave(), hull->get_pitch(), 0.0};
    }

    /// Diagnostic: force/moment breakdown at the CURRENT attitude and speed, with no
    /// equilibrium solve. Returns the same terms find_equilibrium balances, so a trim
    /// sweep reveals which moment contribution sets the equilibrium trim.
    /// [static_lift, static_torque, wave_lift, wave_torque, noblesse_lift, noblesse_torque,
    ///  morabito_lift, morabito_torque, w_plan, lift_press, torque_press, total_lift, total_torque]
    std::vector<double> debug_force_balance(double speed) {
        // Report at the hull's CURRENT pose (no equilibrium solve): re-apply it through the shared
        // core, which reproduces this diagnostic's resistance-band, sectional-free balance exactly.
        PoseForces const pf = residual_core(hull->get_heave(), hull->get_pitch(), speed,
                                            BlendBand::Resistance, false, SectionalMode::None, nullptr);
        return {pf.static_lift, pf.static_torque, pf.wave_lift, pf.wave_torque,
                pf.noblesse_lift, pf.noblesse_torque, pf.morabito_lift, pf.morabito_torque,
                pf.w_plan, pf.lift_press, pf.torque_press, pf.total_lift, pf.total_torque};
    }

    /// Faithful copy of find_equilibrium's eval_at residual (DEFAULT path: sectional-delta
    /// mode OFF) at an arbitrary (heave, pitch). Returns {lift, torque} — the exact two
    /// quantities find_equilibrium drives to zero — so an external (e.g. Python) root-finder
    /// can be benchmarked against the production nested-1-D solver on the identical residual.
    /// NOT used by find_equilibrium; added only for the convergence-method comparison.
    std::vector<double> residual_at(double heave, double pitch, double speed) {
        // Warp gate + static capture at the DESIGN pose (mirrors the find_equilibrium setup).
        hull->set_heave(0.0); hull->set_pitch(0.0); prepare();
        bool use_sec = false;
        if (cfg.use_sectional_attitude) {
            sectional.capture_static();
            double const x_ap = hull->get_ap(), x_fp = hull->get_fp();
            double const bA = hull->get_deadrise_at_x(x_ap + WARP_SAMPLE_AFT_FRAC * (x_fp - x_ap));
            double const bF = hull->get_deadrise_at_x(x_ap + WARP_SAMPLE_FWD_FRAC * (x_fp - x_ap));
            double const warp_spread = (bF > 0.0 && bA > 0.0) ? (bF - bA) : 0.0;
            use_sec = warp_spread > cfg.sectional_warp_floor;
        }
        // Replacement when sectional-delta mode is off (the default and the benchmarked path); when
        // delta mode is on the original applied no sectional here, so map that to None. Attitude band,
        // exactly as eval_at uses.
        SectionalMode const mode = cfg.sectional_delta_mode ? SectionalMode::None : SectionalMode::Replacement;
        PoseForces const pf = residual_core(heave, pitch, speed,
                                            BlendBand::Attitude, use_sec, mode, nullptr);
        return {pf.total_lift, pf.total_torque};
    }

    bool is_ok() const {
        return true;
    }

    static std::string get_version() {
#ifdef MICHELL_VERSION
        return std::string(MICHELL_VERSION);
#else
        return std::string("development");
#endif
    }

};

#endif // RESISTANCE_H
