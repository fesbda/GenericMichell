#ifndef HULL_H
#define HULL_H

#include "environment.h"
#include "mesh.h"
#include <limits>

/// Calculation of propreties and hydrostatics of an immersed object
class Hull {

private:

    static constexpr double empty_double = std::numeric_limits<double>::quiet_NaN();
    static constexpr glm::dvec3 empty_dvec3 = {empty_double, empty_double, empty_double};
    static bool is_empty(double const& x) { return std::isnan(x); }
    static bool is_empty(glm::dvec3 const& v) { return std::isnan(v.x) || std::isnan(v.y) || std::isnan(v.z); }

    std::shared_ptr<Mesh> mesh;
    std::shared_ptr<Environment> env;
    Transformation transform;

    // must be inited to perform all calculations
    int y_side, x_dir;

    // mass, buoyancy
    bool user_set_mass_centroid = false;
    glm::dvec3 mass_centroid;
    glm::dvec3 buoyancy_centroid;
    double displacement;
    double reference_displacement = empty_double;  // fixed design displacement for planing loading
    double reference_length_wl = empty_double;     // fixed design waterline length (Radojcic/MS envelopes)
    double reference_beam_chine = empty_double;     // fixed design chine beam (Radojcic/MS envelopes)
    double reference_beam_wl = empty_double;        // fixed design waterline beam (Radojcic/MS B/T ratio)
    double reference_draft = empty_double;          // fixed design draft (Radojcic/MS B/T ratio)
    glm::dvec3 bbox_min;
    glm::dvec3 bbox_max;

    // wetted stuff
    double wetted_area;
    glm::dvec3 wetted_centroid;
    glm::dvec3 bbox_min_wet;
    glm::dvec3 bbox_max_wet;

    // stuff at waterline
    double z_waterline = 0.0;
    double length_wl;
    double beam_wl;
    double x_beam_wl;
    double waterplane_area;

    // stuff at main frame
    double draft;
    double max_section_area;

    // perpendiculars
    double x_fp;
    double x_ap;

    // fore shoulder
    double x_shoulder;
    double arc_shoulder;

    // transom
    double beam_transom;
    double area_transom;

    // above-water transverse projected (frontal/windage) silhouette area
    double frontal_area;
    double user_frontal_area = empty_double;  // override (e.g. to add superstructure); NaN = auto

    // planing-form geometry (prismatic approximation, aft-half mean)
    bool user_set_deadrise = false;
    bool user_set_beam_chine = false;
    bool user_set_planing_form = false;
    bool deadrise_saturated = false;  // true when auto-estimate hit the upper-bound clamp
    bool planing_form_value = false;  // cached from geometry (or user override)
    double deadrise;     // radians
    double beam_chine;   // metres (full beam, not half)
    // Robust chine extraction (default OFF = legacy wetted-max-breadth, byte-identical).
    // When ON, the per-section chine is located by the full-section bottom-walk (the same
    // detector as get_deadrise_at_x), so beam_chine/deadrise stay the GEOMETRIC chine values
    // when the static waterline sits above the chine on a flared topside (GPPH-type hulls:
    // legacy over-reads beam 0.662 vs 0.6274 and deadrise 23.6 vs 17.2 deg) or when the
    // chine is dry at the static float (light loadings: legacy under-reads the beam, e.g.
    // Fridsma b30/c304 0.144 vs 0.229 -> standalone Savitsky -24% at speed).
    // is_planing_form keeps classifying by the legacy wetted-breadth beam (beam_chine_cls),
    // so the planing/displacement dispatch is unchanged in either mode.
    bool robust_chine = false;
    double beam_chine_cls;  // legacy wetted-max-breadth beam (classification only)

    // hull material properties (optional)
    double roughness = 0;
    double thickness = 0;

public:

    Hull(std::shared_ptr<Mesh> hull_mesh, std::shared_ptr<Environment> environment, int x_fwd_dir = 1)
        : mesh(hull_mesh)
        , env(environment)
        , x_dir(x_fwd_dir)
    {

        if (std::abs(mesh->get_min().y) > 1.0001 * std::abs(mesh->get_max().y)) {
            y_side = -1;
        } else {
            y_side = 1;
        }

        reset_calculated();

    }

    ~Hull() {
//        std::cout << "Delete hull instance" << std::endl;
    }

    // OPTIONS:

    /// Get the mesh that makes the hull
    std::shared_ptr<Mesh> get_mesh() {
        return mesh;
    }

    /// Get the environment of the hull
    std::shared_ptr<Environment> get_environment() {
        return env;
    }

    /// Get aft-to-fore direction along X
    int get_fwd() const {
        return x_dir;
    }

    /// Get used Y direction
    int get_side() const {
        return y_side;
    }

    // USER INPUT:

    /// Set the coordiantes of the mass centroid
    void set_cg(double x, double z) {
        mass_centroid.x = x;
        mass_centroid.y = 0;
        mass_centroid.z = z;
        transform.set_centre(mass_centroid);
        user_set_mass_centroid = true;
    }

    /// Get the coordinates of the mass centroid
    glm::dvec3 get_cg(bool transformed) {

        if (is_empty(mass_centroid)) {
            mass_centroid = get_cb();
            mass_centroid.z = get_waterline();
            user_set_mass_centroid = false;
        }

        return transformed ? transform.get(mass_centroid) : mass_centroid;
    }

    /// Set the vertical coordinate of the waterplane for the calculation
    void set_waterline(double z) {
        z_waterline = z;
        reset_calculated();
    }

    /// Get the vertical coordinate of the waterplane
    double get_waterline() {
        return z_waterline;
    }

    /// Set the heave/sinkage of the hull
    void set_heave(double z) {
        transform.set_translation({0, 0, z});
        reset_calculated();
    }

    /// Get the heave/sinkage of the hull
    double get_heave() const {
        return transform.get_translation().z;
    }

    /// Set the pitch/trim of the hull
    void set_pitch(double angle) {
        transform.set_rotation({0, angle, 0});
        reset_calculated();
    }

    /// Get the pitch/trim of the hull
    double get_pitch() const {
        return transform.get_rotation().y;
    }

    /// Set the hull surface roughness
    void set_roughness(double value) {
        roughness = value;
    }

    /// Get the hull surface roughness
    double get_roughness() const {
        return roughness;
    }

    Transformation get_transform() const {
        return transform;
    }

    glm::dvec3 get_vertex(std::size_t id) {
        return mesh->get_vertex(id, transform);
    }

//    glm::dvec3 get_normal(std::size_t vertex_id) {
//        return mesh->get_normal(vertex_id, transform);
//    }

    // SOLVERS:

    /// Solver to automatically set the sinkage by specifying the displacement
    bool find_heave(double target_displacement, double tolerance=1e-3) {

        if (target_displacement <= 0) {
            return false; // please check your input
        }

        if (std::abs(get_displacement() - target_displacement) / target_displacement < tolerance) {
            return true; // already at the target sinkage
        }

        if (target_displacement > get_displacement() && get_max(false).z < z_waterline) {
            return false; // whole volume already immersed, cannot get some more of it
        }

        if (target_displacement < get_displacement() && get_min(false).z > z_waterline) {
            return false; // whole volume already above waterline, cannot remove more of it
        }

        bool is_found = false;

        double const start_search = target_displacement > get_displacement() ? get_heave() - (get_max(false).z - z_waterline) : get_heave();
        double const end_search   = target_displacement > get_displacement() ? get_heave() : get_heave() + (z_waterline - get_min(false).z);

        double const found_z = bisection(
            [&,this](double z) { this->set_heave(z); return this->get_displacement(); },
            target_displacement,
            start_search,
            end_search,
            get_size(false).z * 5e-4,
            target_displacement*tolerance,
            &is_found);

        return is_found;
    }

    /// Solver to automatically set the waterline by specifying the displacement
    bool find_waterline(double target_displacement, double tolerance=1e-3) {

//        if (target_displacement <= 0) {
//            return false; // please check your input
//        }

//        if (std::abs(get_displacement() - target_displacement) / target_displacement < tolerance) {
//            return true; // already at the target waterline
//        }

//        bool is_found{false};

//        double const start_search = target_displacement > get_displacement() ? get_waterline() : get_min().z;
//        double const end_search   = target_displacement > get_displacement() ? get_max().z : get_waterline();

//        double const found_z = bisection(
//            [&,this](double z) { this->set_waterline(z); return this->get_displacement(); },
//            target_displacement,
//            start_search,
//            end_search,
//            get_size().z * 1e-3,
//            target_displacement*tolerance,
//            &is_found);

//        std::cout << "WL found = " << is_found << " for displ. = " << target_displacement << ", Z = " << found_z << std::endl;

//        set_waterline(found_z);
//        return is_found;
        return false;
    }


    // GETTERS:

    /// Get the minimum coordinate of the immersed hull
    glm::dvec3 get_min(bool wetted) {

        if (is_empty(bbox_min) || is_empty(bbox_min_wet)) {
            calc_bbox();
        }

        return wetted ? bbox_min_wet : bbox_min;
    }

    /// Get the maximum coordinate of the immersed hull
    glm::dvec3 get_max(bool wetted) {

        if (is_empty(bbox_max) || is_empty(bbox_max_wet)) {
            calc_bbox();
        }

        return wetted ? bbox_max_wet : bbox_max;
    }

    /// Get the bounding box size of the immersed hull
    glm::dvec3 get_size(bool wetted) {

        if (is_empty(bbox_max) || is_empty(bbox_min) || is_empty(bbox_min_wet) || is_empty(bbox_max_wet)) {
            calc_bbox();
        }

        return wetted ? bbox_max_wet - bbox_min_wet : bbox_max - bbox_min;
    }

    /// Get the area of waterplane cutting the hull
    double get_waterplane_area() {

        if (is_empty(waterplane_area)) {
            calc_length_beam_wl();
        }

        return waterplane_area;
    }

    /// Get the coefficient of fullness of waterplane
    double get_waterplane_coef() {

        if (get_length_wl() * get_beam_wl() < 1e-16) {
            return 0;
        }

        return get_waterplane_area() / (get_length_wl() * get_beam_wl());
    }

    /// Get the draught at the main frame location
    double get_draft() {

        if (is_empty(draft)) {
            calc_max_section_area();
        }

        return draft;
    }

    /// Get maximum draught
    double get_draft_max() {
        return std::max(0.0, get_waterline() - get_min(false).z);
    }

    /// Get the length at the waterline
    double get_length_wl() {

        if (is_empty(length_wl)) {
            calc_length_beam_wl();
        }

        return length_wl;
    }

    /// Get reference length
    double get_length_ref() {

        if (get_length_wl() == 0) {
            return get_mesh()->get_size().x;
        }

        return get_length_wl();
    }

    /// Get the maximum breadth of the waterplane section
    double get_beam_wl() {

        if (is_empty(beam_wl)) {
            calc_length_beam_wl();
        }

        return beam_wl;
    }

    /// Get the longitudinal location of the maximum breadth at the waterplane
    double get_x_beam_wl() {

        if (is_empty(x_beam_wl)) {
            calc_length_beam_wl();
        }

        return x_beam_wl;
    }

    /// Get the longitudinal location of the fore shoulder
    double get_x_shoulder() {

        if (is_empty(x_shoulder)) {
            calc_shoulder();
        }

        return x_shoulder;
    }

    /// Get the arc of the waterline section from the forward perpendicular to the shoulder
    double get_arc_shoulder() {

        if (is_empty(arc_shoulder)) {
            calc_arc_shoulder();
        }

        return arc_shoulder;
    }

    /// Get the longitudinal coordinate of the forward perpendicular, i.e. the start of the waterline
    double get_fp() {

        if (is_empty(x_fp)) {
            calc_length_beam_wl();
        }

        return x_fp;
    }

    /// Get the longitudinal coordinate of the aft perpendicular, i.e. the end of the waterline
    double get_ap() {

        if (is_empty(x_ap)) {
            calc_length_beam_wl();
        }

        return x_ap;
    }

    /// Get the breadth of the immersed transom
    double get_beam_transom() {

        if (is_empty(beam_transom)) {
            calc_length_beam_wl();
        }

        return beam_transom;
    }

    /// Get the immersed transom area
    double get_transom_area() {

        if (is_empty(area_transom)) {
            if (is_empty(x_ap)) {
                calc_length_beam_wl();
            }
            calc_transom_area();
        }

        return area_transom;
    }

    /// Does the hull has immersed transom
    bool has_transom() {
        return get_beam_transom() >= 1e-2;
    }

    /// Above-water transverse projected (frontal/windage) silhouette area, m^2, at the
    /// current attitude. The mesh is viewed bow-on: the silhouette half-width at each
    /// height above the waterline is the maximum |y| of the above-water surface there, and
    /// the area is 2*max|y| integrated over height -- so it is the full-beam projected area
    /// whether the mesh is a half- or full-hull. This is the A_VT used by the calm-air
    /// windage drag R_AA = 1/2 * rho_air * C_AA * A_VT * V^2. A user override
    /// (set_frontal_area) takes precedence; use it to add a superstructure/windage area.
    double get_frontal_area() {
        if (!is_empty(user_frontal_area)) {
            return user_frontal_area;
        }
        if (is_empty(frontal_area)) {
            calc_frontal_area();
        }
        return frontal_area;
    }

    /// Override the above-water frontal (windage) area, m^2. Pass a positive value to set it
    /// directly (e.g. measured A_VT including superstructure); pass <= 0 to revert to the
    /// auto value computed from the mesh silhouette.
    void set_frontal_area(double area) {
        user_frontal_area = (area > 0.0) ? area : empty_double;
    }

    /// Get the buoyancy centroid
    glm::dvec3 get_cb() {

        if (is_empty(buoyancy_centroid)) {
            calc_displacement_and_area();
        }

        if (std::isnan(buoyancy_centroid.x) || std::isnan(buoyancy_centroid.z)) {
            return mass_centroid;
        }

        return buoyancy_centroid;
    }

    /// Get the displaced volume
    double get_displacement() {

        if (is_empty(displacement)) {
            calc_displacement_and_area();
        }

        return displacement;
    }

    /// Get the displaced mass
    double get_displaced_mass() {
        return get_displacement() * env->get_density();
    }

    /// Reference (design) displaced volume — the fixed buoyant volume at the static
    /// floating condition. Planing methods must support the constant craft weight,
    /// NOT the instantaneous running displacement, which collapses as the hull rises
    /// onto plane. Capture this at the at-rest attitude before imposing trim/sinkage.
    void capture_reference_displacement() {
        reference_displacement = get_displacement();
        // Capture the design-float waterline geometry too, so the Radojcic/Mercier-Savitsky
        // envelopes reference the fixed design hull rather than the mutable running pose (which
        // otherwise leaks the equilibrium-solver seed into the envelope; see cold-vs-warm).
        reference_length_wl = get_length_wl();
        reference_beam_chine = get_beam_chine();
        reference_beam_wl = get_beam_wl();
        reference_draft = get_draft();
    }

    void set_reference_displacement(double vol) {
        reference_displacement = vol;
    }

    /// Fixed reference displaced volume; falls back to the running displacement
    /// when no reference has been captured (preserves legacy behaviour).
    double get_reference_displacement() {
        if (is_empty(reference_displacement)) {
            return get_displacement();
        }
        return reference_displacement;
    }

    /// Fixed design waterline length / chine beam captured at the reference float; fall back to
    /// the running values when no reference has been captured (preserves legacy behaviour).
    double get_reference_length_wl() {
        if (is_empty(reference_length_wl)) { return get_length_wl(); }
        return reference_length_wl;
    }
    double get_reference_beam_chine() {
        if (is_empty(reference_beam_chine)) { return get_beam_chine(); }
        return reference_beam_chine;
    }
    double get_reference_beam_wl() {
        if (is_empty(reference_beam_wl)) { return get_beam_wl(); }
        return reference_beam_wl;
    }
    double get_reference_draft() {
        if (is_empty(reference_draft)) { return get_draft(); }
        return reference_draft;
    }

    double get_reference_mass() {
        return get_reference_displacement() * env->get_density();
    }

    /// Get the wetted surface area
    double get_wetted_area() {

        if (is_empty(wetted_area)) {
            calc_displacement_and_area();
        }

        return wetted_area;
    }

    glm::dvec3 get_wetted_centroid() {

        if (is_empty(wetted_area)) {
            calc_displacement_and_area();
        }

        return wetted_centroid;
    }

    /// Longitudinal position of the geometric wetted-area centroid and of the CG, in the
    /// same (transformed) frame — for the warp-resolved-trim CoP study (bind-friendly scalars).
    double get_wetted_centroid_x() { return get_wetted_centroid().x; }
    double get_cg_x() { return get_cg(true).x; }

    /// Get the area of the main frame
    double get_max_section_area() {

        if (is_empty(max_section_area)) {
            calc_max_section_area();
        }

        return max_section_area;
    }

    /// Immersed transverse-section area at longitudinal station x (both sides), by the
    /// same per-face trapezoidal integration as the max-section/transom calcs. Robust at
    /// an arbitrary station (unlike calc_transom_area, which slices exactly on the coplanar
    /// transom face and so reads 0 there). Used for the Mercier-Savitsky A_T/A_X.
    double get_section_area_at(double x) {
        double area = 0.0;
        for (auto f : mesh->get_faces()) {
            glm::dvec3 tv[3] = {mesh->get_vertex(f.x, transform), mesh->get_vertex(f.y, transform), mesh->get_vertex(f.z, transform)};
            glm::dvec3 ints[3]; int c = 0;
            for (int i = 0; i < 3; i++) {
                c += intersect_edge_plane_x(tv[i], tv[(i+1)%3], x, ints[c]);
            }
            if (c >= 2 && (ints[0].y * y_side > 0 || ints[1].y * y_side > 0)) {
                if (ints[0].z < z_waterline && ints[1].z < z_waterline) {
                    area += std::abs((ints[0].y + ints[1].y) * (ints[0].z - ints[1].z));
                } else if (ints[0].z > z_waterline && ints[1].z < z_waterline) {
                    if (std::abs(ints[0].z - ints[1].z) < 1e-12) continue;
                    auto yw = ints[1].y + (z_waterline - ints[1].z) / (ints[0].z - ints[1].z) * (ints[0].y - ints[1].y);
                    area += std::abs((yw + ints[1].y) * (z_waterline - ints[1].z));
                } else if (ints[0].z < z_waterline && ints[1].z > z_waterline) {
                    if (std::abs(ints[1].z - ints[0].z) < 1e-12) continue;
                    auto yw = ints[0].y + (z_waterline - ints[0].z) / (ints[1].z - ints[0].z) * (ints[1].y - ints[0].y);
                    area += std::abs((yw + ints[0].y) * (z_waterline - ints[0].z));
                }
            }
        }
        return area;   // = 2 * sum(0.5*...) over both sides
    }

    /// Get the coefficient of fullness of the main frame
    double get_max_section_coef() {

        if (get_beam_wl() * get_draft() < 1e-16) {
            return 0.0;
        }

        return get_max_section_area() / (get_beam_wl() * get_draft());
    }

    /// Get the block coefficient, i.e. fullness of form compared to immersed box
    double get_block_coef() {

        if (get_beam_wl() * get_draft() * get_length_wl() < 1e-16) {
            return 0.0;
        }

        return get_displacement() / (get_beam_wl() * get_draft() * get_length_wl());
    }

    /// Get the prismatic coefficient, i.e. fullness of form compared to immersed prism
    double get_prismatic_coef() {

        if (get_max_section_area() * get_length_wl() < 1e-16) {
            return 0.0;
        }

        return get_displacement() / (get_max_section_area() * get_length_wl());
    }

    /// Get the form factor correction based on the Reynolds number change
    double get_form_factor_correction(double from_speed, double to_speed) {

        auto FCF_from = std::tanh(0.015064 * std::pow(std::log10(get_rn(from_speed)), 2.6752));
        auto FCF_to   = std::tanh(0.015064 * std::pow(std::log10(get_rn(to_speed)), 2.6752));
        return FCF_to / FCF_from;
    }

    /// Get Reynolds number for a speed and a specified length
    double get_rn(double speed, double length) {
        return speed * length / env->get_viscosity();
    }

    /// Get Reynolds number for a speed and length at the waterline
    double get_rn(double speed) {
        return get_rn(speed, get_length_ref());
    }

    /// Get Froude number for a speed and a specified length
    double get_fn(double speed, double length) {
        return speed / std::sqrt(env->get_gravity() * length);
    }

    /// Get Froude number for a speed and length at the waterline
    double get_fn(double speed) {
        return get_fn(speed, get_length_ref());
    }

    /// Get volume-based Froude number for a speed and the displacement
    double get_fn_vol(double speed) {
        return speed / std::sqrt(env->get_gravity() * std::pow(get_displacement(), 1.0/3.0));
    }

    /// Get speed from a Froude number
    double get_speed(double fn) {
        return fn * std::sqrt(get_length_ref() * env->get_gravity());
    }

    /// Get the non-dimensional coefficient from a force value
    double get_coef(double speed, double force) {
        return force / (0.5 * env->get_density() * get_wetted_area() * sq(speed));
    }

    /// Get the force value from the non-dimensional coefficient
    double get_force(double speed, double coef) {
        return coef * (0.5 * env->get_density() * get_wetted_area() * sq(speed));
    }

    /// Set the deadrise angle (radians) explicitly. Pass NaN to revert to mesh-based estimate.
    void set_deadrise(double angle) {
        if (std::isnan(angle)) {
            user_set_deadrise = false;
            deadrise = empty_double;
        } else {
            user_set_deadrise = true;
            deadrise = angle;
        }
    }

    /// Set the chine beam (metres) explicitly. Pass NaN to revert to mesh-based estimate.
    void set_beam_chine(double b) {
        if (std::isnan(b)) {
            user_set_beam_chine = false;
            beam_chine = empty_double;
        } else {
            user_set_beam_chine = true;
            beam_chine = b;
        }
    }

    /// Robust chine extraction (see calc_deadrise_and_chine). Default OFF preserves the
    /// legacy wetted-breadth extraction byte-identically; toggling invalidates the cached
    /// mesh-derived deadrise/beam so the next getter recomputes in the requested mode.
    void set_robust_chine(bool on) {
        if (robust_chine == on) {
            return;
        }
        robust_chine = on;
        if (!user_set_deadrise) deadrise = empty_double;
        if (!user_set_beam_chine) beam_chine = empty_double;
        beam_chine_cls = empty_double;
    }

    /// Override the planing-form classification (true: hull is treated as planing-suitable
    /// for the Auto blend; false: hull is treated as pure displacement). Pass NaN-equivalent
    /// (call clear_planing_form_override) to revert to mesh-based detection.
    void set_planing_form(bool value) {
        user_set_planing_form = true;
        planing_form_value = value;
    }

    /// Revert is_planing_form() to mesh-based auto-detection.
    void clear_planing_form_override() {
        user_set_planing_form = false;
    }

    /// True if the hull's geometry resembles a planing form (low-deadrise prismatic
    /// hull with a chine at the waterline). The Auto resistance dispatch only blends
    /// in Morabito/Savitsky contributions when this is true; pure displacement hulls
    /// (round bilge, no chine) avoid Savitsky's invalid-geometry output. Refined by
    /// two indicators: (1) auto-detected deadrise stays below the upper validity bound,
    /// (2) chine beam matches waterline beam (i.e. the wetted polyline's widest point
    /// is at the waterline, characteristic of chined hulls).
    bool is_planing_form() {
        if (user_set_planing_form) {
            return planing_form_value;
        }
        if (is_empty(deadrise) || is_empty(beam_chine)) {
            calc_deadrise_and_chine();
        }
        if (deadrise_saturated) {
            return false;
        }
        double const bwl = get_beam_wl();
        if (bwl <= 1e-6) {
            return false;
        }
        // Reject degenerate near-dry waterplanes: a flat-bottomed displacement ship
        // (e.g. KCS, KVLCC2) grazed right at the keel presents a zero-depth wetted
        // sliver whose beam trivially equals the waterline beam and whose deadrise
        // reads ~0, mimicking a chine-at-waterline planing section. A genuine planing
        // form has a wetted body of real depth, so require the draft to be a
        // meaningful fraction of the waterline beam before classifying.
        if (get_draft() < 0.03 * bwl) {
            return false;
        }
        // Classification uses the legacy wetted-breadth beam in either extraction mode
        // (beam_chine_cls == beam_chine when robust_chine is off): the "chine at the
        // waterline" indicator is DEFINED on the wetted section, so the robust geometric
        // chine (which can sit above or below the waterline) must not change dispatch.
        double const b_cls = user_set_beam_chine ? beam_chine : beam_chine_cls;
        return (b_cls / bwl) >= 0.99;
    }

    /// Mean deadrise angle (radians) over the aft half of the wetted hull, for planing methods.
    double get_deadrise() {

        if (is_empty(deadrise)) {
            calc_deadrise_and_chine();
        }

        return deadrise;
    }

    /// Mean chine beam (metres) over the aft half of the wetted hull, for planing methods.
    double get_beam_chine() {

        if (is_empty(beam_chine)) {
            calc_deadrise_and_chine();
        }

        return beam_chine;
    }

    /// Bottom walk shared by get_deadrise_at_x and the robust chine extraction
    /// (calc_deadrise_and_chine with robust_chine ON). Deadrise = slope of the BOTTOM,
    /// from keel to chine. The full section also contains the near-vertical topside above
    /// the chine; a plain max-|y| chine then grabbed the topside and reported a wild
    /// "deadrise" (w1 42.9 vs true ~19 deg; w2's chine knuckle is not even a mesh vertex,
    /// the section jumps bottom->topside). Robust fix: collapse to a half-section (unique
    /// |y|, lowest z), then walk outward from the keel along the bottom and stop at the
    /// first near-vertical step (slope > 1, i.e. > 45 deg, far steeper than any real
    /// deadrise) -- the chine is the last bottom point before that step.
    /// Returns false on a degenerate section. Exact code move from get_deadrise_at_x
    /// (2026-07-10), so the local-deadrise readout is bit-identical.
    static bool bottom_chine_walk(std::vector<double> const& ys, std::vector<double> const& zs,
                                  double* y_keel, double* z_keel,
                                  double* y_chine, double* z_chine) {
        if (ys.size() < 2) {
            return false;
        }
        std::vector<std::pair<double,double>> half;  // (|y|, z), outward
        {
            std::vector<std::pair<double,double>> pts;
            pts.reserve(ys.size());
            for (std::size_t k = 0; k < ys.size(); k++) pts.push_back({std::abs(ys[k]), zs[k]});
            std::sort(pts.begin(), pts.end());
            for (auto const& p : pts) {
                if (!half.empty() && p.first - half.back().first < 1e-4) {
                    half.back().second = std::min(half.back().second, p.second);  // bottom at this |y|
                } else {
                    half.push_back(p);
                }
            }
        }
        if (half.size() < 2) {
            return false;
        }
        *y_keel = half.front().first; *z_keel = half.front().second;
        *y_chine = *y_keel; *z_chine = *z_keel;
        for (std::size_t k = 1; k < half.size(); k++) {
            double const dyk = half[k].first - half[k-1].first;
            double const dzk = half[k].second - half[k-1].second;
            if (dyk < 1e-6) continue;
            if (dzk / dyk > 1.0) break;            // near-vertical topside reached
            *y_chine = half[k].first; *z_chine = half[k].second;
        }
        return true;
    }

    /// Bottom walk on the true lower ENVELOPE of the raw section segments (robust-chine
    /// extraction only; the legacy bottom_chine_walk above is kept exact for
    /// get_deadrise_at_x byte-stability). Point-based folds fail here in two ways
    /// (both observed): the deck's edge crossings interleave between bottom points at
    /// interior |y| and trip the >45 deg stop mid-bottom (the bottom is a long straight
    /// segment with no endpoint in that bin), and vertical topsides collapse to one
    /// arbitrary z in slice_section's y->z map. So: take SEGMENTS, fold to |y| (splitting
    /// any segment that crosses the centreline), evaluate each on the |y|-bins it spans,
    /// and keep the lowest z per bin — the true lower boundary. A near-vertical column
    /// (topside at the chine) deposits its LOWEST z, so the outboard-most wetted bin is
    /// the chine itself; deck bins beyond a flared chine rise at slope > 1 and stop the
    /// walk exactly there (GPPH: flare above the chine must not widen the beam).
    static bool bottom_chine_walk_envelope(std::vector<double> const& segs,
                                           double* y_keel, double* z_keel,
                                           double* y_chine, double* z_chine) {
        std::size_t const n = segs.size() / 4;
        if (n < 2) {
            return false;
        }
        double y_max = 0.0;
        for (std::size_t i = 0; i < n; i++) {
            y_max = std::max(y_max, std::max(std::abs(segs[4*i]), std::abs(segs[4*i+2])));
        }
        if (y_max < 1e-4) {
            return false;
        }
        int const NB = 256;
        double const w = y_max / NB;
        std::vector<double> env_y(std::size_t(NB), 0.0), env_z(std::size_t(NB), 0.0);
        std::vector<bool> used(std::size_t(NB), false);
        auto deposit = [&](double a0, double z0, double a1, double z1) {
            if (a0 > a1) { std::swap(a0, a1); std::swap(z0, z1); }
            int b0 = std::max(0, std::min(NB - 1, int(a0 / y_max * NB)));
            int b1 = std::max(0, std::min(NB - 1, int(a1 / y_max * NB)));
            for (int b = b0; b <= b1; b++) {
                double const yc = std::min(a1, std::max(a0, (b + 0.5) * w));
                double const z = (a1 - a0 > 1e-12)
                                 ? z0 + (yc - a0) / (a1 - a0) * (z1 - z0)
                                 : std::min(z0, z1);   // vertical column: its lowest point
                if (!used[std::size_t(b)] || z < env_z[std::size_t(b)]) {
                    env_y[std::size_t(b)] = yc; env_z[std::size_t(b)] = z;
                    used[std::size_t(b)] = true;
                }
            }
        };
        for (std::size_t i = 0; i < n; i++) {
            double const y0 = segs[4*i],     z0 = segs[4*i+1];
            double const y1 = segs[4*i+2],   z1 = segs[4*i+3];
            if (y0 * y1 < 0.0) {
                // crosses the centreline: split at y=0 so the |y| fold stays monotone
                double const t = y0 / (y0 - y1);
                double const zc = z0 + t * (z1 - z0);
                deposit(std::abs(y0), z0, 0.0, zc);
                deposit(0.0, zc, std::abs(y1), z1);
            } else {
                deposit(std::abs(y0), z0, std::abs(y1), z1);
            }
        }
        int first = -1;
        for (int b = 0; b < NB; b++) { if (used[std::size_t(b)]) { first = b; break; } }
        if (first < 0) {
            return false;
        }
        *y_keel = env_y[std::size_t(first)]; *z_keel = env_z[std::size_t(first)];
        *y_chine = *y_keel; *z_chine = *z_keel;
        int prev = first;
        for (int b = first + 1; b < NB; b++) {
            if (!used[std::size_t(b)]) continue;
            double const dyk = env_y[std::size_t(b)] - env_y[std::size_t(prev)];
            double const dzk = env_z[std::size_t(b)] - env_z[std::size_t(prev)];
            prev = b;
            if (dyk < 1e-6) continue;
            if (dzk / dyk > 1.0) break;            // near-vertical topside / flare reached
            *y_chine = env_y[std::size_t(b)]; *z_chine = env_z[std::size_t(b)];
        }
        return true;
    }

    /// Local deadrise angle (radians) of the section at longitudinal station x.
    /// Used for warped hulls, whose deadrise varies along the length: the planing
    /// pressure must be integrated against the *local* deadrise, not a single mean.
    /// Returns a negative value if the section is degenerate (caller should fall back).
    double get_deadrise_at_x(double x) {
        std::vector<double> ys, zs;
        slice_section(x, ys, zs, false);   // full section (attitude-stable geometry)
        double y_keel, z_keel, y_chine, z_chine;
        if (!bottom_chine_walk(ys, zs, &y_keel, &z_keel, &y_chine, &z_chine)) {
            return -1.0;
        }
        double const dy = y_chine - y_keel;
        if (dy <= 1e-4) {
            return -1.0;
        }
        return std::atan2(std::max(0.0, z_chine - z_keel), dy);
    }

    /// Wetted-width-weighted mean of the section-local deadrise over the currently-wetted
    /// bottom at the running attitude. This is the deadrise the hull is actually planing
    /// on: at planing trim the forward sections lift clear of the water, so only the lower
    /// aft sections are wetted and counted. For a constant-deadrise hull it equals the
    /// constant; for a warped hull it is well below the static area-mean deadrise, which
    /// is what the (mean-deadrise) Savitsky trim relation over-predicts from. Each station
    /// is weighted by its wetted half-beam (its share of the planform). Sampled, not tuned.
    /// Returns <0 if no wetted section is found (caller falls back to the area mean).
    double get_wetted_mean_deadrise() {
        double const fwd = double(x_dir);
        double const x_transom = (fwd > 0) ? get_min(false).x : get_max(false).x;
        double const x_bow     = (fwd > 0) ? get_max(false).x : get_min(false).x;
        double const span = std::abs(x_bow - x_transom);
        if (span < 1e-6) {
            return -1.0;
        }
        int const N = 40;
        double num = 0.0, den = 0.0;
        std::vector<double> ys, zs;
        for (int i = 0; i < N; i++) {
            double const d = (double(i) + 0.5) / double(N) * span;   // fwd of transom
            double const x = x_transom + fwd * d;
            slice_section(x, ys, zs, true);                          // wetted section
            double hb = 0.0;
            for (double y : ys) hb = std::max(hb, std::abs(y));      // wetted half-beam
            if (hb <= 1e-5) continue;                                // dry station
            double const bl = get_deadrise_at_x(x);                  // attitude-stable local β
            if (bl < 0.0) continue;
            double const w = hb;
            num += w * std::min(bl, radians(50.0));
            den += w;
        }
        return (den > 1e-9) ? (num / den) : -1.0;
    }

    /// Waterline half-angle of entrance (degrees): the angle the design waterline makes
    /// with the centreline at the bow. Estimated by a through-origin fit of the waterline
    /// half-beam vs distance aft of the stem over the forward ~15% of L_wl (half-beam is
    /// zero at the stem). Used by the Mercier-Savitsky (1973) pre-planing method.
    double get_entrance_half_angle() {
        double const L = get_length_wl();
        if (L < 1e-6) {
            return 10.0;
        }
        double const fwd = double(x_dir);
        double const x_stem = (fwd > 0) ? get_fp() : get_ap();   // bow end at the waterline
        int const N = 8;
        double sxx = 0.0, sxy = 0.0;
        std::vector<double> ys, zs;
        for (int i = 1; i <= N; i++) {
            double const d = (0.15 * L) * double(i) / double(N);  // distance aft of stem
            double const x = x_stem - fwd * d;
            slice_section(x, ys, zs, true);                       // wetted section at x
            double hb = 0.0;
            for (double y : ys) hb = std::max(hb, std::abs(y));   // waterline half-beam
            if (hb <= 0.0) continue;
            sxx += d * d; sxy += d * hb;
        }
        if (sxx < 1e-12) {
            return 10.0;
        }
        double const slope = sxy / sxx;                          // d(half-beam)/d(distance)
        double const ie = degrees(std::atan(std::max(0.0, slope)));
        return std::min(60.0, std::max(1.0, ie));
    }

    /// Get coordinates of a transverse section polyline (cut at constant x).
    void slice_section(double x, std::vector<double>& out_ys, std::vector<double>& out_zs, bool below_waterline) {

        std::map<double, double> section;  // key: y, value: z

        for (auto f : mesh->get_faces()) {

            glm::dvec3 tv[3] = {mesh->get_vertex(f.x, transform), mesh->get_vertex(f.y, transform), mesh->get_vertex(f.z, transform)};
            glm::dvec3 ints[3];
            int c = 0;
            for (int i = 0; i < 3; i++) {
                c += intersect_edge_plane_x(tv[i], tv[(i+1)%3], x, ints[c]);
            }
            if (c != 2) {
                continue;
            }

            // y_side filtering: include the user-active half only
            if (ints[0].y * y_side < 0 && ints[1].y * y_side < 0) {
                continue;
            }

            if (!below_waterline) {
                section[round(ints[0].y, 4)] = ints[0].z;
                section[round(ints[1].y, 4)] = ints[1].z;
                continue;
            }

            if (ints[0].z < z_waterline && ints[1].z < z_waterline) {
                section[round(ints[0].y, 4)] = ints[0].z;
                section[round(ints[1].y, 4)] = ints[1].z;
            } else if (ints[0].z > z_waterline && ints[1].z < z_waterline) {
                if (std::abs(ints[0].z - ints[1].z) < 1e-12) continue;
                double const y_at_wl = ints[1].y + (z_waterline - ints[1].z) / (ints[0].z - ints[1].z) * (ints[0].y - ints[1].y);
                section[round(y_at_wl, 4)] = z_waterline;
                section[round(ints[1].y, 4)] = ints[1].z;
            } else if (ints[0].z < z_waterline && ints[1].z > z_waterline) {
                if (std::abs(ints[1].z - ints[0].z) < 1e-12) continue;
                double const y_at_wl = ints[0].y + (z_waterline - ints[0].z) / (ints[1].z - ints[0].z) * (ints[1].y - ints[0].y);
                section[round(y_at_wl, 4)] = z_waterline;
                section[round(ints[0].y, 4)] = ints[0].z;
            }
        }

        out_ys.clear();
        out_zs.clear();

        for (auto const& s : section) {
            out_ys.push_back(s.first);
            out_zs.push_back(s.second);
        }
    }

    /// Raw section SEGMENTS (y0,z0,y1,z1 per crossing face, flattened) at station x.
    /// Preserves the segment topology that slice_section's y->z map collapses: the map
    /// keeps one z per rounded y, so vertical topsides and decks overwrite each other and
    /// the lower boundary cannot be reconstructed from its point list. Needed by the
    /// robust chine extraction's true lower envelope. Full section (no waterline clip).
    void slice_section_segments(double x, std::vector<double>& out) {
        out.clear();
        for (auto f : mesh->get_faces()) {
            glm::dvec3 tv[3] = {mesh->get_vertex(f.x, transform), mesh->get_vertex(f.y, transform), mesh->get_vertex(f.z, transform)};
            glm::dvec3 ints[3];
            int c = 0;
            for (int i = 0; i < 3; i++) {
                c += intersect_edge_plane_x(tv[i], tv[(i+1)%3], x, ints[c]);
            }
            if (c != 2) {
                continue;
            }
            if (ints[0].y * y_side < 0 && ints[1].y * y_side < 0) {
                continue;
            }
            out.push_back(ints[0].y); out.push_back(ints[0].z);
            out.push_back(ints[1].y); out.push_back(ints[1].z);
        }
    }

    /// Get coordinates of a buttock polyline
    /// Warp-resolved-planform study: per-transverse-station wetted longitudinal extent at the
    /// current (posed) attitude. For n stations at y = (k+0.5)/n * 0.5*beam_chine (keel→chine),
    /// returns a flattened [xmin0, xmax0, xmin1, xmax1, ...] of the wetted buttock (posed frame,
    /// metres); dry stations return (0,0). This is the GEOMETRIC wetted planform (the true curved
    /// boundary), to be compared against the kernel's straight-λ triangle.
    std::vector<double> get_wetted_buttock_profile(int n) {
        std::vector<double> out;
        out.reserve(std::size_t(2) * n);
        double const b = std::max(get_beam_chine(), 1e-6);
        std::vector<double> xs, zs;
        for (int k = 0; k < n; k++) {
            double const y = (double(k) + 0.5) / double(n) * 0.5 * b;
            xs.clear(); zs.clear();
            slice_buttock(y, xs, zs, true);
            if (xs.empty()) { out.push_back(0.0); out.push_back(0.0); continue; }
            double xmin = xs.front(), xmax = xs.front();
            for (double x : xs) { xmin = std::min(xmin, x); xmax = std::max(xmax, x); }
            out.push_back(xmin); out.push_back(xmax);
        }
        return out;
    }

    void slice_buttock(double y, std::vector<double>& out_xs, std::vector<double>& out_zs, bool below_waterline) {

        y *= y_side;

        std::map<double, double> buttock;

        for (auto f : mesh->get_faces()) {

            glm::dvec3 tv[3] = {mesh->get_vertex(f.x, transform), mesh->get_vertex(f.y, transform), mesh->get_vertex(f.z, transform)};
            glm::dvec3 ints[3];
            int c = 0;
            for (int i = 0; i < 3; i++) {
                c += intersect_edge_plane_y(tv[i], tv[(i+1)%3], y, ints[c]);
            }
            if (c != 2) {
                continue;
            }

            if (!below_waterline) {
                buttock[round(ints[0].x, 3)] = ints[0].z;
                buttock[round(ints[1].x, 3)] = ints[1].z;
                continue;
            }

            if (ints[0].z < z_waterline && ints[1].z < z_waterline) {
                buttock[round(ints[0].x, 3)] = ints[0].z;
                buttock[round(ints[1].x, 3)] = ints[1].z;
            } else if (ints[0].z > z_waterline && ints[1].z < z_waterline) {
                if (std::abs(ints[0].z - ints[1].z) < 1e-12) continue;
                double const x_waterline = ints[1].x + (z_waterline - ints[1].z) / (ints[0].z - ints[1].z) * (ints[0].x - ints[1].x);
                buttock[round(x_waterline, 3)] = z_waterline;
                buttock[round(ints[1].x, 3)] = ints[1].z;
            } else if (ints[0].z < z_waterline && ints[1].z > z_waterline) {
                if (std::abs(ints[1].z - ints[0].z) < 1e-12) continue;
                double const x_waterline = ints[0].x + (z_waterline - ints[0].z) / (ints[1].z - ints[0].z) * (ints[1].x - ints[0].x);
                buttock[round(x_waterline, 3)] = z_waterline;
                buttock[round(ints[0].x, 3)] = ints[0].z;
            }
        }

        out_xs.clear();
        out_zs.clear();

        for (auto const b : buttock) {
            out_xs.push_back(b.first);
            out_zs.push_back(b.second);
        }

    }

    /// Savitsky (2012) warped-hull effective-trim input: angle difference (radians) between
    /// the quarter-beam buttock (bottom at |y| = beam_chine/4) and the keel buttock (|y|=0),
    /// each least-squares fitted over 7 stations across [x_station - hw, x_station + hw].
    /// Built from DIRECT bottom-height samples (section half-collapse interpolated at the
    /// offset) — point z-samples carry only ~mm chord-sag noise on a coarse fan-triangulated
    /// mesh, whereas differentiating the fitted deadrise angle (or slicing buttock polylines,
    /// which deck crossings corrupt) amplifies mesh noise to whole degrees. Positive = the
    /// quarter-beam buttock runs further bow-up than the keel: a warped bottom carries its
    /// aft body at a higher incidence than the keel pitch reports. On a prismatic bottom the
    /// buttocks are parallel and the delta vanishes (to FP rounding). Angle differences are
    /// invariant under the pitch rotation, so posed-frame sampling is pose-pure (cold==warm
    /// safe). Returns 0 when the sections are degenerate (caller falls back to keel trim).
    double get_quarter_buttock_delta(double x_station, double half_window) {
        double const b = std::max(get_beam_chine(), 1e-6);
        double const f = double(get_fwd());
        double const h = std::max(half_window, 1e-3);
        int const K = 7;
        double sx = 0.0, sxx = 0.0, sq = 0.0, sxq = 0.0, sk = 0.0, sxk = 0.0;
        int n = 0;
        for (int i = 0; i < K; i++) {
            double const x = x_station + h * (2.0 * double(i) / double(K - 1) - 1.0);
            double const zq = bottom_z_at(x, 0.25 * b);
            double const zk = bottom_z_at(x, 0.0);
            if (std::isnan(zq) || std::isnan(zk)) {
                continue;
            }
            sx += x; sxx += x * x;
            sq += zq; sxq += x * zq;
            sk += zk; sxk += x * zk;
            n++;
        }
        if (n < 3) {
            return 0.0;
        }
        double const det = double(n) * sxx - sx * sx;
        if (det < 1e-12) {
            return 0.0;
        }
        double const slope_q = (double(n) * sxq - sx * sq) / det;
        double const slope_k = (double(n) * sxk - sx * sk) / det;
        return std::atan(f * slope_q) - std::atan(f * slope_k);
    }

    /// Mesh-robust local deadrise for the Savitsky-2012 effective values: mean of the valid
    /// get_deadrise_at_x samples over 7 stations across [x_station - hw, x_station + hw]
    /// (the per-station chine walk carries degree-level noise on coarse meshes; the window
    /// mean recovers the station value the method wants). < 0 when no station is valid.
    double get_deadrise_about(double x_station, double half_window) {
        double const h = std::max(half_window, 1e-3);
        int const K = 7;
        double sum = 0.0;
        int n = 0;
        for (int i = 0; i < K; i++) {
            double const x = x_station + h * (2.0 * double(i) / double(K - 1) - 1.0);
            double const bl = get_deadrise_at_x(x);
            if (bl > 0.0) { sum += bl; n++; }
        }
        return (n > 0) ? sum / double(n) : -1.0;
    }

private:

    /// z of the hull BOTTOM at station x and transverse offset |y| = y_off: the full section
    /// collapsed to a half-section (unique |y|, lowest z — the same collapse
    /// get_deadrise_at_x applies), linearly interpolated at y_off. y_off = 0 returns the
    /// local keel point. NaN when the section is degenerate or narrower than y_off.
    double bottom_z_at(double x, double y_off) {
        std::vector<double> ys, zs;
        slice_section(x, ys, zs, false);
        if (ys.size() < 2) {
            return empty_double;
        }
        std::vector<std::pair<double, double>> half;  // (|y|, z), outward
        {
            std::vector<std::pair<double, double>> pts;
            pts.reserve(ys.size());
            for (std::size_t k = 0; k < ys.size(); k++) pts.push_back({std::abs(ys[k]), zs[k]});
            std::sort(pts.begin(), pts.end());
            for (auto const& p : pts) {
                if (!half.empty() && p.first - half.back().first < 1e-4) {
                    half.back().second = std::min(half.back().second, p.second);
                } else {
                    half.push_back(p);
                }
            }
        }
        if (half.size() < 2) {
            return empty_double;
        }
        // Walk outward keeping only points that continue the bottom line (slope <= 1 vs the
        // last KEPT point, the same near-vertical threshold get_deadrise_at_x uses). This
        // SKIPS deck-interior fan-diagonal crossings — isolated points far above the bottom
        // that interleave the section at |y| values with no bottom counterpart — instead of
        // stopping at them; past the real chine everything rises steeper than 1 and is
        // skipped, so the kept polyline is the bottom only.
        std::vector<std::pair<double, double>> bot;
        bot.push_back(half.front());
        for (std::size_t k = 1; k < half.size(); k++) {
            double const dy = half[k].first - bot.back().first;
            if (dy < 1e-6) {
                continue;
            }
            if ((half[k].second - bot.back().second) / dy > 1.0) {
                continue;
            }
            bot.push_back(half[k]);
        }
        if (y_off <= bot.front().first) {
            return bot.front().second;
        }
        if (y_off > bot.back().first) {
            return empty_double;
        }
        for (std::size_t k = 1; k < bot.size(); k++) {
            if (bot[k].first >= y_off) {
                double const dy = bot[k].first - bot[k - 1].first;
                if (dy < 1e-12) { return bot[k].second; }
                double const t = (y_off - bot[k - 1].first) / dy;
                return bot[k - 1].second + t * (bot[k].second - bot[k - 1].second);
            }
        }
        return empty_double;
    }


    void reset_calculated() {

        // mass, buoyancy
        if (!user_set_mass_centroid) {
            mass_centroid = empty_dvec3;
        }
        buoyancy_centroid = empty_dvec3;
        displacement = empty_double;
        bbox_min = empty_dvec3;
        bbox_max = empty_dvec3;

        // wetted stuff
        wetted_area = empty_double;
        wetted_centroid = empty_dvec3;
        bbox_min_wet = empty_dvec3;
        bbox_max_wet = empty_dvec3;

        // stuff at waterline
        length_wl = empty_double;
        beam_wl = empty_double;
        x_beam_wl = empty_double;
        waterplane_area = empty_double;

        // stuff at main frame
        draft = empty_double;
        max_section_area = empty_double;

        // perpendiculars
        x_fp = empty_double;
        x_ap = empty_double;

        // fore shoulder
        x_shoulder = empty_double;
        arc_shoulder = empty_double;

        // transom
        beam_transom = empty_double;
        area_transom = empty_double;

        // windage
        frontal_area = empty_double;   // user_frontal_area is a user override, not reset

        // planing-form geometry
        if (!user_set_deadrise) {
            deadrise = empty_double;
        }
        if (!user_set_beam_chine) {
            beam_chine = empty_double;
        }
        beam_chine_cls = empty_double;

    }

    void calc_bbox() {

        constexpr double dmax = std::numeric_limits<double>::max();
        constexpr double dlow = std::numeric_limits<double>::lowest();
        bbox_min = bbox_min_wet = {dmax, dmax, dmax};
        bbox_max = bbox_max_wet = {dlow, dlow, dlow};

        std::size_t const n_vertices = mesh->get_vertices().size();

        for (std::size_t i = 0; i < n_vertices; i++) {

            auto const vtx = get_vertex(i);

            if (vtx.y * y_side < 0.0) {
                continue;
            }

            bbox_min = glm::min(bbox_min, vtx);
            bbox_max = glm::max(bbox_max, vtx);

            if (vtx.z > z_waterline) {
                continue;
            }

            bbox_min_wet = glm::min(bbox_min_wet, vtx);
            bbox_max_wet = glm::max(bbox_max_wet, vtx);
        }
    }

    void calc_displacement_and_area() {

        displacement = wetted_area = 0;
        buoyancy_centroid = {};
        wetted_centroid = {};
//        int n_tris = 0;

        auto const& normals = mesh->get_face_normals();

        for (std::size_t fi = 0; fi < mesh->get_faces().size(); fi++) {

            auto f = mesh->get_faces()[fi];
            glm::dvec3 tv[3] = {mesh->get_vertex(f.x, transform), mesh->get_vertex(f.y, transform), mesh->get_vertex(f.z, transform)};
            auto centre = (tv[0]+tv[1]+tv[2]) / 3.0;

            // triangle centre is on the wrong side
            if (centre.y * y_side < 0.0) {
                continue;
            }

            auto z_min = std::min({tv[0].z, tv[1].z, tv[2].z});
            auto z_max = std::max({tv[0].z, tv[1].z, tv[2].z});

            // triangle above waterline
            if (z_min > z_waterline) {
                continue;
            }

            // triangle intersects waterline -> clip to the below-water polygon and
            // fan-triangulate it. Walk the edges in order, emitting each below-water vertex
            // followed by the crossing point on that edge; this yields the ordered (convex)
            // below-water polygon (3 verts when one vertex is below, 4 when two are). Using
            // a per-below-vertex fan to the same two crossings (the previous scheme) does NOT
            // tile the trapezoid when two vertices are below — it overlaps high and misses a
            // wedge low, systematically UNDER-counting displacement and wetted area by a term
            // ~ (above-water rise)-dependent (negligible for fine meshes, large where a tall
            // triangle straddles the waterline, e.g. a single chine->deck topside facet).
            if (z_max > z_waterline && z_max > z_min+1e-6) {

                glm::dvec3 poly[4]; int np = 0;
                for (int i = 0; i < 3; i++) {
                    glm::dvec3 const& a = tv[i];
                    glm::dvec3 const& bb = tv[(i+1)%3];
                    if (a.z < z_waterline && np < 4) {
                        poly[np++] = a;
                    }
                    glm::dvec3 ip;
                    if (intersect_edge_plane_z(a, bb, z_waterline, ip) && np < 4) {
                        poly[np++] = ip;
                    }
                }

                for (int k = 1; k + 1 < np; k++) {   // fan poly[0]-poly[k]-poly[k+1]
                    glm::dvec3 tc = (poly[0] + poly[k] + poly[k+1]) / 3.0;
                    auto d_area = 2 * 0.5 * glm::length(glm::cross(poly[k]-poly[0], poly[k+1]-poly[0]));
                    auto d_vol = std::abs(tc.y * d_area * normals[fi].y);
                    wetted_area += d_area;
                    displacement += d_vol;
                    buoyancy_centroid.x += d_vol * tc.x;
                    buoyancy_centroid.z += d_vol * tc.z;
                    wetted_centroid.x += d_area * tc.x;
                    wetted_centroid.z += d_area * tc.z;
                }

            // triangle below waterline
            } else {

                auto d_area = 2 * mesh->get_face_areas()[fi];
                auto d_vol = std::abs(centre.y * d_area * normals[fi].y);
                wetted_area += d_area;
                displacement += d_vol;
                buoyancy_centroid.x += d_vol * centre.x;
                buoyancy_centroid.z += d_vol * centre.z;
                wetted_centroid.x += d_area * centre.x;
                wetted_centroid.z += d_area * centre.z;
            }

//            n_tris++;
        }

        if (displacement > 0.0) {
            buoyancy_centroid /= displacement;
        }

        if (wetted_area > 0.0) {
            wetted_centroid /= wetted_area;
        }

        // if user hasn't set CG, put it above CB at the waterline
        if (!user_set_mass_centroid) {
            mass_centroid = buoyancy_centroid;
            mass_centroid.z = z_waterline;
        }

        // include material thickness
        displacement += thickness * wetted_area;

        // todo fix for area cos of discrete triangles
//        auto discrete_fix = std::max(1.0, 1.001 - 9.888e-9*n_tris);
//        displacement *= cub(discrete_fix);
//        wetted_area *= sq(discrete_fix);
    }

    void calc_max_section_area() {

        max_section_area = 0;
        draft = 0;

        if (is_empty(x_beam_wl)) {
            calc_length_beam_wl();
        }

        for (auto f : mesh->get_faces()) {

            glm::dvec3 tv[3] = {mesh->get_vertex(f.x, transform), mesh->get_vertex(f.y, transform), mesh->get_vertex(f.z, transform)};
            glm::dvec3 ints[3]; int c = 0;
            for (int i = 0; i < 3; i++) {
                c += intersect_edge_plane_x(tv[i], tv[(i+1)%3], x_beam_wl, ints[c]);
            }

            // if found two intersected edges of a face, on the proper side of the ship
            if (c >= 2 && (ints[0].y * y_side > 0 || ints[1].y * y_side > 0)) {

                if (ints[0].z < z_waterline && ints[1].z < z_waterline) {
                    max_section_area += 2 * 0.5 * std::abs((ints[0].y + ints[1].y) * (ints[0].z - ints[1].z));
                } else if (ints[0].z > z_waterline && ints[1].z < z_waterline) {
                    if (std::abs(ints[0].z - ints[1].z) < 1e-12) continue;
                    auto y_waterline = ints[1].y + (z_waterline - ints[1].z) / (ints[0].z - ints[1].z) * (ints[0].y - ints[1].y);
                    max_section_area += 2 * 0.5 * std::abs((y_waterline + ints[1].y) * (z_waterline - ints[1].z));
                } else if (ints[0].z < z_waterline && ints[1].z > z_waterline) {
                    if (std::abs(ints[1].z - ints[0].z) < 1e-12) continue;
                    auto y_waterline = ints[0].y + (z_waterline - ints[0].z) / (ints[1].z - ints[0].z) * (ints[1].y - ints[0].y);
                    max_section_area += 2 * 0.5 * std::abs((y_waterline + ints[0].y) * (z_waterline - ints[0].z));
                }
            }

            // calc draft
            for (int i = 0; i < c; i++) {
                if (ints[i].z < z_waterline && z_waterline - ints[i].z > draft) { // no check for y side?
                    draft = z_waterline - ints[i].z;
                }
            }
        }

    }

    void calc_transom_area() {

        area_transom = 0.0;

        if (!has_transom()) {
            return;
        }

        for (auto f : mesh->get_faces()) {

            glm::dvec3 tv[3] = {mesh->get_vertex(f.x, transform), mesh->get_vertex(f.y, transform), mesh->get_vertex(f.z, transform)};
            glm::dvec3 ints[3];
            int c = 0;
            for (int i = 0; i < 3; i++) {
                c += intersect_edge_plane_x(tv[i], tv[(i+1)%3], x_ap, ints[c]);
            }

            if (c >= 2 && (ints[0].y * y_side > 0 || ints[1].y * y_side > 0)) {

                if (ints[0].z < z_waterline && ints[1].z < z_waterline) {
                    area_transom += 2 * 0.5 * std::abs((ints[0].y + ints[1].y) * (ints[0].z - ints[1].z));
                } else if (ints[0].z > z_waterline && ints[1].z < z_waterline) {
                    if (std::abs(ints[0].z - ints[1].z) < 1e-12) continue;
                    auto y_waterline = ints[1].y + (z_waterline - ints[1].z) / (ints[0].z - ints[1].z) * (ints[0].y - ints[1].y);
                    area_transom += 2 * 0.5 * std::abs((y_waterline + ints[1].y) * (z_waterline - ints[1].z));
                } else if (ints[0].z < z_waterline && ints[1].z > z_waterline) {
                    if (std::abs(ints[1].z - ints[0].z) < 1e-12) continue;
                    auto y_waterline = ints[0].y + (z_waterline - ints[0].z) / (ints[1].z - ints[0].z) * (ints[1].y - ints[0].y);
                    area_transom += 2 * 0.5 * std::abs((y_waterline + ints[0].y) * (z_waterline - ints[0].z));
                }
            }
        }
    }

    void calc_frontal_area() {

        frontal_area = 0.0;

        double const z_wl = z_waterline;
        double const z_top = get_max(false).z;
        if (!(z_top > z_wl)) {
            return;   // no freeboard -> no above-water area
        }

        // Bow-on silhouette: the half-width of the above-water hull at each height band is the
        // maximum |y| of the surface there, and the frontal area is 2*max|y| summed over the
        // bands (full beam from a half- or full-hull mesh). We sweep mesh EDGES through the
        // height bins -- not just vertices -- so a vertical topside spanning two distant rows
        // (e.g. the Fridsma/Taunton prisms, which carry only a chine-row and a deck-row vertex
        // with nothing between) fills every band it crosses rather than vanishing between bins.
        // Robust to mesh density without a union-of-triangles projection -- adequate for the
        // windage term.
        int const n_bins = 100;
        double const dz = (z_top - z_wl) / n_bins;
        std::vector<double> half_beam(n_bins, 0.0);

        auto sweep_edge = [&](glm::dvec3 const& va, glm::dvec3 const& vb) {
            double za = va.z, zb = vb.z, ya = va.y, yb = vb.y;
            if (za > zb) { std::swap(za, zb); std::swap(ya, yb); }   // now za <= zb
            double const lo = std::max(za, z_wl);
            double const hi = std::min(zb, z_top);
            if (hi <= lo) {
                return;   // edge does not reach above the waterline
            }
            double const span = zb - za;
            int klo = static_cast<int>((lo - z_wl) / dz);
            int khi = static_cast<int>((hi - z_wl) / dz);
            if (klo < 0) klo = 0;
            if (khi >= n_bins) khi = n_bins - 1;
            for (int k = klo; k <= khi; k++) {
                double const zc = z_wl + (k + 0.5) * dz;   // bin centre
                if (zc < lo || zc > hi) {
                    continue;
                }
                double const t = (span > 1e-12) ? (zc - za) / span : 0.0;
                double const y = std::abs(ya + t * (yb - ya));
                if (y > half_beam[k]) {
                    half_beam[k] = y;
                }
            }
        };

        for (auto f : mesh->get_faces()) {
            glm::dvec3 const tv[3] = {
                mesh->get_vertex(f.x, transform),
                mesh->get_vertex(f.y, transform),
                mesh->get_vertex(f.z, transform)
            };
            sweep_edge(tv[0], tv[1]);
            sweep_edge(tv[1], tv[2]);
            sweep_edge(tv[2], tv[0]);
        }

        for (int b = 0; b < n_bins; b++) {
            frontal_area += 2.0 * half_beam[b] * dz;
        }
    }

    void calc_length_beam_wl() {

        x_fp = -x_dir * 1e12;
        x_ap = -x_fp;
        beam_transom = 0.0;
        beam_wl = 0;
        waterplane_area = 0;
        length_wl = 0;

        auto const x_midship = (get_min(false).x + get_max(false).x) * 0.5;
        auto const x_eps = mesh->get_size().x * 0.25;
        auto const x_beam_min = x_midship - x_eps;
        auto const x_beam_max = x_midship + x_eps;

        for (auto f : mesh->get_faces()) {

            glm::dvec3 tv[3] = {mesh->get_vertex(f.x, transform), mesh->get_vertex(f.y, transform), mesh->get_vertex(f.z, transform)};
            glm::dvec3 ints[3];
            int c = 0;
            for (int i = 0; i < 3; i++) {
                c += intersect_edge_plane_z(tv[i], tv[(i+1)%3], z_waterline, ints[c]);
            }

            // if found two intersected edges of a face, on the proper side of the ship
            if (c >= 2 && (ints[0].y * y_side > 0 || ints[1].y * y_side > 0)) {

                if (x_dir > 0) {
                    x_fp = std::max(x_fp, std::max(ints[0].x, ints[1].x));
                    if (std::min(ints[0].x, ints[1].x) < x_ap) {
                        x_ap = std::min(ints[0].x, ints[1].x);
                        beam_transom = std::min(ints[0].y * y_side, ints[1].y * y_side);
                    }
                } else {
                    if (std::max(ints[0].x, ints[1].x) > x_ap) {
                        x_ap = std::max(ints[0].x, ints[1].x);
                        beam_transom = std::min(ints[0].y * y_side, ints[1].y * y_side);
                    }
                    x_fp = std::min(x_fp, std::min(ints[0].x, ints[1].x));
                }

                waterplane_area += 2 * 0.5 * std::abs((ints[0].y + ints[1].y) * (ints[0].x - ints[1].x));
            }

            for(int i = 0; i < c; i++) {
                if (ints[i].x > x_beam_min && ints[i].x < x_beam_max && ints[i].y * y_side > beam_wl) {
                    x_beam_wl = ints[i].x + 1e-5;
                    beam_wl = ints[i].y * y_side;
                }
            }

            /*
            c = 0;
            for (int i = 0; i < 3; i++) {
                c += intersect_edge_plane_z(tv[i], tv[(i+1)%3], z_waterline+0.002, ints[c]);
            }
            */


        }

        if (std::abs(x_fp - x_ap) > get_mesh()->get_size().x) {

            if (x_dir > 0) {
                x_fp = get_mesh()->get_max().x;
                x_ap = get_mesh()->get_min().x;
            } else {
                x_fp = get_mesh()->get_min().x;
                x_ap = get_mesh()->get_max().x;
            }

        } /*else*/ {

            beam_wl *= 2;
            length_wl = std::abs(x_fp - x_ap);

            beam_transom = 2 * std::max(0.0, std::abs(beam_transom));
            if (beam_transom < 1e-2 * beam_wl) {
                beam_transom = 0;
            }

            //std::cout << "beam of transom " << beam_transom << std::endl;
        }


    }

    void calc_shoulder() {

        if (is_empty(beam_wl)) {
            calc_length_beam_wl();
        }

        x_shoulder = x_fp;
        double dist_from_fp = get_length_wl();

        double const target_half_beam = (0.5 * beam_wl) * 0.99;

        for (auto f : mesh->get_faces()) {

            glm::dvec3 tv[3] = {mesh->get_vertex(f.x, transform), mesh->get_vertex(f.y, transform), mesh->get_vertex(f.z, transform)};
            glm::dvec3 ints[3];
            int c = 0;
            for (int i = 0; i < 3; i++) {
                c += intersect_edge_plane_z(tv[i], tv[(i+1)%3], z_waterline, ints[c]);
            }

            for(int i = 0; i < c; i++) {
                if (ints[i].y * y_side > target_half_beam && std::abs(x_fp - ints[i].x) < dist_from_fp) {
                    x_shoulder = ints[i].x;
                    dist_from_fp = std::abs(x_fp - ints[i].x);
                }
            }
        }

    }

    void calc_arc_shoulder() {

        if (is_empty(x_shoulder)) {
            calc_shoulder();
        }

        double arc = 0;
        std::map<double, double> xy_map;

        for (auto f : mesh->get_faces()) {

            glm::dvec3 tv[3] = {mesh->get_vertex(f.x, transform), mesh->get_vertex(f.y, transform), mesh->get_vertex(f.z, transform)};
            glm::dvec3 ints[3]; int c = 0;
            for (int i = 0; i < 3; i++) {
                c += intersect_edge_plane_z(tv[i], tv[(i+1)%3], z_waterline, ints[c]);
            }

            // if found two intersected edges of a face, on the proper side of the ship
            if (c >= 2 && (ints[0].y * y_side > 0 || ints[1].y * y_side > 0)) {
                // is front of x
                auto mid = (ints[0].x + ints[1].x) * 0.5;
                if (mid > std::min(x_shoulder, x_fp) && mid < std::max(x_shoulder, x_fp)) {
                    xy_map.insert(std::make_pair(ints[0].x, ints[0].y));
                    xy_map.insert(std::make_pair(ints[1].x, ints[1].y));
                }
            }
        }

        for (auto it = xy_map.begin(); it != xy_map.end(); ++it) {
            std::map<double, double>::const_iterator nx = it; ++nx;
            if (nx == xy_map.end()) {
                continue;
            }
            arc += std::sqrt(sq(nx->first - it->first) + sq(nx->second - it->second));
        }

        arc_shoulder = arc;
    }

    /// Compute mean deadrise angle and chine beam over the aft half of the wetted hull,
    /// for the prismatic-hull approximation used by Savitsky / Morabito planing methods.
    /// Honours user-supplied overrides; only fills the values the user has not set.
    void calc_deadrise_and_chine() {

        bool const need_beta = !user_set_deadrise;
        bool const need_b = !user_set_beam_chine;
        if (!need_beta && !need_b) {
            return;
        }

        // Aft half of the wetted hull: from midpoint of wetted X-range to the transom (aft end).
        auto const wmin = get_min(true);
        auto const wmax = get_max(true);
        double const x_mid = 0.5 * (wmin.x + wmax.x);
        double const x_aft = (x_dir > 0) ? wmin.x : wmax.x;
        if (std::abs(x_aft - x_mid) < 1e-6) {
            if (need_beta) deadrise = radians(10.0);  // benign fallback
            if (need_b) beam_chine = std::max(get_beam_wl(), 1e-3);
            beam_chine_cls = std::max(get_beam_wl(), 1e-3);
            return;
        }

        // Sample N sections, skipping the very ends to avoid degenerate slices.
        int const N = 11;
        double const margin = 0.05 * (x_mid - x_aft);  // 5% of the half-length
        double sum_beta = 0.0, sum_b = 0.0;
        int n_beta = 0, n_b = 0;
        double sum_b_cls = 0.0;                        // legacy wetted-breadth (classification)
        int n_b_cls = 0;

        std::vector<double> ys, zs;
        for (int i = 0; i < N; i++) {
            double const t = (i + 0.5) / double(N);
            double const x = (x_aft + margin) + t * ((x_mid - margin) - (x_aft + margin));
            ys.clear(); zs.clear();
            slice_section(x, ys, zs, true);
            if (ys.size() < 2) {
                continue;
            }

            // Half-beam at this section: largest |y| of the wetted polyline.
            double y_chine = 0.0;
            double z_chine = z_waterline;
            for (std::size_t k = 0; k < ys.size(); k++) {
                double const ay = std::abs(ys[k]);
                if (ay > y_chine) {
                    y_chine = ay;
                    z_chine = zs[k];
                }
            }
            if (y_chine < 1e-4) {
                continue;
            }
            sum_b_cls += 2.0 * y_chine;
            n_b_cls++;

            // Keel point: smallest |y| (keel-line). Take its z.
            double y_keel = std::numeric_limits<double>::max();
            double z_keel = z_chine;
            for (std::size_t k = 0; k < ys.size(); k++) {
                double const ay = std::abs(ys[k]);
                if (ay < y_keel) {
                    y_keel = ay;
                    z_keel = zs[k];
                }
            }

            // Robust mode: replace the wetted-breadth chine with the geometric chine from
            // the full-section bottom-walk (the get_deadrise_at_x detector). The wetted
            // slice above is WRONG in two float configurations: waterline above the chine
            // (grabs the flared topside: beam and deadrise over-read) and chine dry at
            // rest (returns the waterline breadth: beam under-read).
            if (robust_chine) {
                std::vector<double> segs;
                slice_section_segments(x, segs);
                double yk, zk, yc, zc;
                if (!bottom_chine_walk_envelope(segs, &yk, &zk, &yc, &zc)) {
                    continue;
                }
                y_keel = yk; z_keel = zk; y_chine = yc; z_chine = zc;
                if (y_chine < 1e-4) {
                    continue;
                }
            }

            double const dy = y_chine - y_keel;
            double const dz = z_chine - z_keel;
            if (dy > 1e-4) {
                sum_beta += std::atan2(std::max(0.0, dz), dy);
                n_beta++;
            }
            sum_b += 2.0 * y_chine;
            n_b++;
        }

        if (need_beta) {
            // Mean over sampled sections, clamped to Morabito's validity range [0, 40 deg].
            double const beta_raw = (n_beta > 0) ? (sum_beta / n_beta) : radians(10.0);
            double const beta_max = radians(40.0);
            deadrise_saturated = (beta_raw >= beta_max - 1e-6);
            deadrise = std::max(0.0, std::min(beta_raw, beta_max));
        }

        if (need_b) {
            double b_est = (n_b > 0) ? (sum_b / n_b) : std::max(get_beam_wl(), 1e-3);
            beam_chine = std::max(b_est, 1e-3);
        }
        beam_chine_cls = std::max((n_b_cls > 0) ? (sum_b_cls / n_b_cls)
                                                : std::max(get_beam_wl(), 1e-3), 1e-3);
    }

};

#endif // HULL_H
