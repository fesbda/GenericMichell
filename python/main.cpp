#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "../src/resistance.h"
#include "../src/sectional.h"
#include "../src/propeller.h"

namespace py = pybind11;

int main() {

    std::cout << "Hello from Michell Python!" << std::endl;
    return 0;
}

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
PYBIND11_MAKE_OPAQUE(std::vector<double>)


PYBIND11_MODULE(michell, resistance) {

    py::bind_vector<std::vector<double>>(resistance, "VectorDouble");

    py::class_<glm::dvec3>(resistance, "dvec3")
        .def_readwrite("x", &glm::dvec3::x)
        .def_readwrite("y", &glm::dvec3::y)
        .def_readwrite("z", &glm::dvec3::z)
    ;

    py::class_<Environment, std::shared_ptr<Environment>>(resistance, "Environment")
        .def(py::init<>())
        .def_property("gravity", &Environment::get_gravity, &Environment::set_gravity)
        .def_property("density", &Environment::get_density, &Environment::set_density)
        .def_property("air_density", &Environment::get_air_density, &Environment::set_air_density)
        .def_property("viscosity", &Environment::get_viscosity, &Environment::set_viscosity)
        .def_property("depth", &Environment::get_depth, &Environment::set_depth)
    ;

    py::class_<Mesh, std::shared_ptr<Mesh>>(resistance, "Mesh")
        .def(py::init<>())
        .def("add_vertex", &Mesh::add_vertex)
        .def("add_face", &Mesh::add_face)
        .def("update", &Mesh::update)
    ;

    py::class_<BiLinearInterpolator>(resistance, "BiLinearInterpolator")
        .def(py::init<>())
        .def_readonly("xi", &BiLinearInterpolator::xi)
        .def_readonly("yi", &BiLinearInterpolator::yi)
        .def_readonly("zi", &BiLinearInterpolator::fi)
        .def("get", &BiLinearInterpolator::operator())
    ;

    py::class_<PropellerWakeField>(resistance, "PropellerWakeField")
        .def(py::init<>())
        .def_readonly("center", &PropellerWakeField::center)
        .def_readonly("R_inner", &PropellerWakeField::R_inner)
        .def_readonly("R_outer", &PropellerWakeField::R_outer)
        .def_readonly("n_radial", &PropellerWakeField::n_radial)
        .def_readonly("n_azimuthal", &PropellerWakeField::n_azimuthal)
        .def_readonly("ship_speed", &PropellerWakeField::ship_speed)
        .def_readonly("r", &PropellerWakeField::r)
        .def_readonly("theta", &PropellerWakeField::theta)
        .def_readonly("axial", &PropellerWakeField::axial)
        .def_readonly("radial", &PropellerWakeField::radial)
        .def_readonly("tangential", &PropellerWakeField::tangential)
        .def_readonly("mean_axial", &PropellerWakeField::mean_axial)
        .def_readonly("wake_fraction", &PropellerWakeField::wake_fraction)
        .def_readonly("axial_std", &PropellerWakeField::axial_std)
        .def_readonly("axial_min", &PropellerWakeField::axial_min)
        .def_readonly("axial_max", &PropellerWakeField::axial_max)
        .def_readonly("non_uniformity", &PropellerWakeField::non_uniformity)
        .def_readonly("radial_rms", &PropellerWakeField::radial_rms)
        .def_readonly("swirl_rms", &PropellerWakeField::swirl_rms)
        .def_readonly("volumetric_flow", &PropellerWakeField::volumetric_flow)
        .def_readonly("disc_area", &PropellerWakeField::disc_area)
    ;

    py::class_<Hull, std::shared_ptr<Hull>>(resistance, "Hull")
        .def(py::init<std::shared_ptr<Mesh>, std::shared_ptr<Environment>, int>())
        .def("get_mesh", &Hull::get_mesh)
        .def("get_environment", &Hull::get_environment)
        .def("get_fwd", &Hull::get_fwd)
        .def("get_side", &Hull::get_side)
        .def("set_cg", &Hull::set_cg)
        .def("set_waterline", &Hull::set_waterline)
        .def("get_waterline", &Hull::get_waterline)
        .def("get_waterplane_area", &Hull::get_waterplane_area)
        .def("get_waterplane_coef", &Hull::get_waterplane_coef)
        .def("set_roughness", &Hull::set_roughness)
        .def("get_roughness", &Hull::get_roughness)
        .def("set_heave", &Hull::set_heave)
        .def("get_heave", &Hull::get_heave)
        .def("set_pitch", &Hull::set_pitch)
        .def("get_pitch", &Hull::get_pitch)
        .def("find_waterline", &Hull::find_waterline)
        .def("get_draft", &Hull::get_draft)
        .def("get_length_wl", &Hull::get_length_wl)
        .def("get_beam_wl", &Hull::get_beam_wl)
        .def("get_x_beam_wl", &Hull::get_x_beam_wl)
        .def("get_x_shoulder", &Hull::get_x_shoulder)
        .def("get_fp", &Hull::get_fp)
        .def("get_ap", &Hull::get_ap)
        .def("get_cb", &Hull::get_cb)
        .def("get_displacement", &Hull::get_displacement)
        .def("get_displaced_mass", &Hull::get_displaced_mass)
        .def("capture_reference_displacement", &Hull::capture_reference_displacement)
        .def("set_reference_displacement", &Hull::set_reference_displacement)
        .def("get_reference_displacement", &Hull::get_reference_displacement)
        .def("get_reference_mass", &Hull::get_reference_mass)
        .def("get_wetted_area", &Hull::get_wetted_area)
        .def("get_max_section_area", &Hull::get_max_section_area)
        .def("get_max_section_coef", &Hull::get_max_section_coef)
        .def("get_block_coef", &Hull::get_block_coef)
        .def("get_prismatic_coef", &Hull::get_prismatic_coef)
        .def("get_form_factor_correction", &Hull::get_form_factor_correction)
        .def("get_rn", (double (Hull::*)(double,double))(&Hull::get_rn))
        .def("get_rn_wl", (double (Hull::*)(double))(&Hull::get_rn))
        .def("get_fn", (double (Hull::*)(double,double))(&Hull::get_fn))
        .def("get_fn_wl", (double (Hull::*)(double))(&Hull::get_fn))
        .def("get_fn_vol", &Hull::get_fn_vol)
        .def("get_speed", &Hull::get_speed)
        .def("get_coef", &Hull::get_coef)
        .def("get_force", &Hull::get_force)
        .def("get_arc_shoulder", &Hull::get_arc_shoulder)
        .def("get_deadrise", &Hull::get_deadrise)
        .def("get_deadrise_at_x", &Hull::get_deadrise_at_x)
        .def("get_wetted_centroid_x", &Hull::get_wetted_centroid_x)
        .def("get_wetted_buttock_profile", &Hull::get_wetted_buttock_profile)
        .def("get_cg_x", &Hull::get_cg_x)
        .def("get_wetted_mean_deadrise", &Hull::get_wetted_mean_deadrise)
        .def("get_quarter_buttock_delta", &Hull::get_quarter_buttock_delta,
             py::arg("x_station"), py::arg("half_window"))
        .def("get_entrance_half_angle", &Hull::get_entrance_half_angle)
        .def("get_transom_area", &Hull::get_transom_area)
        .def("get_beam_transom", &Hull::get_beam_transom)
        .def("has_transom", &Hull::has_transom)
        .def("get_frontal_area", &Hull::get_frontal_area)
        .def("set_frontal_area", &Hull::set_frontal_area)
        .def("get_section_area_at", &Hull::get_section_area_at)
        .def("get_beam_chine", &Hull::get_beam_chine)
        .def("set_deadrise", &Hull::set_deadrise)
        .def("set_beam_chine", &Hull::set_beam_chine)
        .def("set_robust_chine", &Hull::set_robust_chine)
        .def("is_planing_form", &Hull::is_planing_form)
        .def("set_planing_form", &Hull::set_planing_form)
        .def("clear_planing_form_override", &Hull::clear_planing_form_override)
        .def("slice_section", [](Hull& h, double x, bool below_waterline) {
            std::vector<double> ys, zs;
            h.slice_section(x, ys, zs, below_waterline);
            return std::make_pair(ys, zs);
        }, py::arg("x"), py::arg("below_waterline") = true)
        .def("slice_buttock", [](Hull& h, double y, bool below_waterline) {
            std::vector<double> xs, zs;
            h.slice_buttock(y, xs, zs, below_waterline);
            return std::make_pair(xs, zs);
        }, py::arg("y"), py::arg("below_waterline") = true)
    ;

    py::enum_<ResistanceMethod>(resistance, "ResistanceMethod")
        .value("Auto", ResistanceMethod::Auto)
        .value("Holtrop", ResistanceMethod::Holtrop)
        .value("Unified", ResistanceMethod::Unified)
        .value("Integral", ResistanceMethod::Integral)
        .value("Morabito", ResistanceMethod::Morabito)
        .value("Savitsky", ResistanceMethod::Savitsky)
        .export_values()
    ;

    py::enum_<TransomShapeMethod>(resistance, "TransomMethod")
        .value("DoctorsDynamic", TransomShapeMethod::DoctorsDynamic)
        .value("DoctorsStatic", TransomShapeMethod::DoctorsStatic)
        .value("Robards", TransomShapeMethod::Robards)
        .value("Mixed", TransomShapeMethod::Mixed)
        .export_values()
    ;

    py::enum_<ResistanceGridPrecision>(resistance, "ResistanceGridPrecision")
        .value("Coarse", ResistanceGridPrecision::Coarse)
        .value("Fine", ResistanceGridPrecision::Fine)
        .value("Finest", ResistanceGridPrecision::Finest)
        .export_values()
    ;

    py::class_<Sectional, std::shared_ptr<Sectional>>(resistance, "Sectional")
        .def(py::init([](std::shared_ptr<Hull> h) {
            return std::make_shared<Sectional>(h, h->get_environment()); }))
        .def("set_nx", &Sectional::set_nx)
        .def("set_smooth_frac", &Sectional::set_smooth_frac)
        .def("set_scale", &Sectional::set_scale)
        .def("get_scale", &Sectional::get_scale)
        .def("set_tloc_cap", &Sectional::set_tloc_cap)
        .def("get_tloc_cap", &Sectional::get_tloc_cap)
        .def("set_transom_fix", &Sectional::set_transom_fix)
        .def("get_transom_fix", &Sectional::get_transom_fix)
        .def("set_transom_ventilation", &Sectional::set_transom_ventilation)
        .def("get_transom_ventilation", &Sectional::get_transom_ventilation)
        .def("set_vent_beam_frac", &Sectional::set_vent_beam_frac)
        .def("set_crossflow", &Sectional::set_crossflow)
        .def("get_crossflow", &Sectional::get_crossflow)
        .def("set_forward_term", &Sectional::set_forward_term)
        .def("get_forward_term", &Sectional::get_forward_term)
        .def("get_forward_cs", &Sectional::get_forward_cs)
        .def("lift", &Sectional::lift)
        .def("added_mass", [](Sectional& s, double U) {
            auto r = s.added_mass(U);
            return py::make_tuple(r.La, r.cop_a, r.lcb, r.L, r.x_t, r.nwet, r.Lc, r.cop_c, r.Dc, r.Da); })
        .def("debug_profiles", [](Sectional& s, double U) {
            std::vector<double> xs, d, beta, tloc, cw, q, fa;
            s.debug_profiles(U, xs, d, beta, tloc, cw, q, fa);
            return py::make_tuple(xs, d, beta, tloc, cw, q, fa); })
    ;

    py::class_<Resistance, std::shared_ptr<Resistance>>(resistance, "Resistance")
        .def_static("get_version", &Resistance::get_version)
        .def(py::init<std::shared_ptr<Hull>>())
        .def("get_hull", &Resistance::get_hull)
        .def("sync_attitude_from_original", &Resistance::sync_attitude_from_original)
        .def("is_ok", &Resistance::is_ok)
        .def("prepare", &Resistance::prepare)
        .def("set_demihull", &Resistance::set_demihull)
        .def("is_demihull", &Resistance::is_demihull)
        .def("get_friction_coef", &Resistance::get_friction_coef)
        .def("get_roughness_coef", &Resistance::get_roughness_coef)
        .def("get_friction_resistance", &Resistance::get_friction_resistance)
        .def("get_pressure_resistance", &Resistance::get_pressure_resistance)
        .def("get_viscous_pressure_resistance", &Resistance::get_viscous_pressure_resistance)
        .def("get_wave_resistance", &Resistance::get_wave_resistance)
        .def("get_total_resistance", &Resistance::get_total_resistance)
        .def("get_total_resistance_equilibrium", &Resistance::get_total_resistance_equilibrium,
             py::arg("speed"), py::arg("method"), py::arg("iters") = 150)
        .def("set_params", &Resistance::set_params)
        .def("set_deflection_table", [](Resistance& r, py::list xs, py::list ys) {
            std::vector<double> vxs, vys;
            for (auto h : xs) vxs.push_back(h.cast<double>());
            for (auto h : ys) vys.push_back(h.cast<double>());
            r.set_deflection_table(vxs, vys);
        })
        .def("set_deflection_field", [](Resistance& r, py::list field) {
            std::vector<double> v;
            for (auto h : field) v.push_back(h.cast<double>());
            r.set_deflection_field(v);
        })
        .def("set_deflection_fields", [](Resistance& r, py::list michell_field, py::list localflow_field) {
            std::vector<double> vm, vl;
            for (auto h : michell_field) vm.push_back(h.cast<double>());
            for (auto h : localflow_field) vl.push_back(h.cast<double>());
            r.set_deflection_fields(vm, vl);
        })
        .def("get_grid", &Resistance::get_grid, py::return_value_policy::reference)
        .def("get_localflow_grid", &Resistance::get_localflow_grid, py::return_value_policy::reference)
        .def("get_centerline_wave_elevation", [](Resistance& r, double speed, py::list xs) {
            std::vector<double> v;
            for (auto h : xs) v.push_back(h.cast<double>());
            return r.get_centerline_wave_elevation(speed, v);
        })
        .def("get_centerline_wave_elevation_michell", [](Resistance& r, double speed, py::list xs) {
            std::vector<double> v;
            for (auto h : xs) v.push_back(h.cast<double>());
            return r.get_centerline_wave_elevation_michell(speed, v);
        })
        .def("get_centerline_wave_elevation_local", [](Resistance& r, double speed, py::list xs) {
            std::vector<double> v;
            for (auto h : xs) v.push_back(h.cast<double>());
            return r.get_centerline_wave_elevation_local(speed, v);
        })
        .def("get_wave_field", [](Resistance& r, double speed, py::list xs, py::list ys, double viscosity, int n_theta) {
            std::vector<double> vxs, vys;
            for (auto h : xs) vxs.push_back(h.cast<double>());
            for (auto h : ys) vys.push_back(h.cast<double>());
            return r.get_wave_field(speed, vxs, vys, viscosity, n_theta);
        }, py::arg("speed"), py::arg("xs"), py::arg("ys"),
           py::arg("viscosity") = 5e-3, py::arg("n_theta") = 401)
        .def("set_phase_model", &Resistance::set_phase_model)
        .def("get_wave_pressure_field", &Resistance::get_wave_pressure_field)
        .def("get_local_pressure_field", &Resistance::get_local_pressure_field)
        .def("get_morabito_pressure_field", &Resistance::get_morabito_pressure_field)
        .def("get_morabito_force", &Resistance::get_morabito_force)
        .def("get_morabito_split", &Resistance::get_morabito_split)
        .def("planing_weight", &Resistance::planing_weight)
        .def("get_planing_lambda", &Resistance::get_planing_lambda)
        .def("get_planing_lambda_geometric", &Resistance::get_planing_lambda_geometric)
        .def("get_planing_stagnation_angle", &Resistance::get_planing_stagnation_angle)
        .def("set_planing_lambda_cap", &Resistance::set_planing_lambda_cap)
        .def("set_true_planform", &Resistance::set_true_planform)
        .def("set_savitsky2012", &Resistance::set_savitsky2012,
             py::arg("mode"), py::arg("station_frac") = 1.0)
        .def("set_morabito_hydrostatic", &Resistance::set_morabito_hydrostatic)
        .def("set_dynamic_cop_fraction", &Resistance::set_dynamic_cop_fraction)
        .def("set_cop_savitsky", &Resistance::set_cop_savitsky)
        .def("set_cop_cv_ramp", &Resistance::set_cop_cv_ramp,
             py::arg("slope"), py::arg("ref") = 2.0, py::arg("warp_floor") = 0.13)
        .def("set_warp_lift_closure", &Resistance::set_warp_lift_closure,
             py::arg("lift_amp"), py::arg("moment_amp") = 0.0,
             py::arg("cv_lo") = 1.0, py::arg("cv_hi") = 2.4)
        .def("set_blend_use_fnb", &Resistance::set_blend_use_fnb)
        .def("set_planing_legacy", &Resistance::set_planing_legacy)
        .def("set_planing_pileup", &Resistance::set_planing_pileup)
        .def("set_planing_pileup_auto", &Resistance::set_planing_pileup_auto)
        .def("set_planing_friction_factor", &Resistance::set_planing_friction_factor)
        .def("set_planing_blend", &Resistance::set_planing_blend)
        .def("set_planing_blend_attitude", &Resistance::set_planing_blend_attitude)
        .def("planing_weight_attitude", &Resistance::planing_weight_attitude)
        .def("set_preplaning_hump", &Resistance::set_preplaning_hump,
             pybind11::arg("coef"), pybind11::arg("fn_lo") = 0.6, pybind11::arg("fn_hi") = 1.9)
        .def("get_preplaning_hump", &Resistance::get_preplaning_hump)
        .def("set_sectional_pressure_drag", &Resistance::set_sectional_pressure_drag)
        .def("get_sectional_pressure_drag", &Resistance::get_sectional_pressure_drag)
        .def("set_sectional_lift_relief", &Resistance::set_sectional_lift_relief,
             pybind11::arg("frac"), pybind11::arg("cv_lo") = 2.0, pybind11::arg("cv_hi") = 3.6)
        .def("set_sectional_pdrag_dynshare", &Resistance::set_sectional_pdrag_dynshare)
        .def("set_pdyn_projection", &Resistance::set_pdyn_projection,
             py::arg("on"), py::arg("delta0_deg") = 1.6,
             py::arg("k_warp") = 0.0, py::arg("cv_lo") = 2.5, py::arg("cv_hi") = 4.5)
        .def("get_pdyn_pressure_drag", &Resistance::get_pdyn_pressure_drag)
        .def("set_effective_buoyancy", &Resistance::set_effective_buoyancy,
             pybind11::arg("frac"), pybind11::arg("cv_lo") = 0.1, pybind11::arg("cv_hi") = 3.5,
             pybind11::arg("torque") = true)
        .def("set_savitsky_trim_weight", &Resistance::set_savitsky_trim_weight)
        .def("set_savitsky_trim_weight_warp", &Resistance::set_savitsky_trim_weight_warp,
             py::arg("w"), py::arg("lo") = 0.05, py::arg("hi") = 0.15)
        .def("set_drag_couple", &Resistance::set_drag_couple,
             py::arg("on"), py::arg("warp_floor") = 0.0, py::arg("warp_hi") = 0.15)
        .def("get_savitsky_trim", &Resistance::get_savitsky_trim,
             py::arg("speed"), py::arg("beta_override") = -1.0)
        .def("cop_test_at_trim", &Resistance::cop_test_at_trim,
             py::arg("speed"), py::arg("trim_deg"))
        .def("cop_components_at_trim", &Resistance::cop_components_at_trim,
             py::arg("speed"), py::arg("trim_deg"))
        .def("set_whisker_spray", &Resistance::set_whisker_spray)
        .def("set_whisker_warp_only", &Resistance::set_whisker_warp_only)
        .def("set_whisker_cvb_cap", &Resistance::set_whisker_cvb_cap,
             py::arg("c"), py::arg("cv_lo") = 4.0, py::arg("cv_hi") = 5.0)
        .def("set_friction_width", &Resistance::set_friction_width)
        .def("set_transom_head_release", &Resistance::set_transom_head_release,
             py::arg("on"), py::arg("amp") = 1.0, py::arg("len_scale") = 1.0)
        .def("get_transom_hollow_params", &Resistance::get_transom_hollow_params,
             py::arg("speed"))
        .def("get_transom_hollow_drop", &Resistance::get_transom_hollow_drop,
             py::arg("speed"), py::arg("X"), py::arg("Y"), py::arg("len_scale") = 1.0)
        .def("get_transom_head_release", &Resistance::get_transom_head_release,
             py::arg("speed"), py::arg("X"), py::arg("Y"), py::arg("Z"),
             py::arg("len_scale") = 1.0)
        .def("set_whisker_spray_turbulent", &Resistance::set_whisker_spray_turbulent)
        .def("get_whisker_spray_turbulent", &Resistance::get_whisker_spray_turbulent)
        .def("set_whisker_spray_root", &Resistance::set_whisker_spray_root)
        .def("get_whisker_spray_root", &Resistance::get_whisker_spray_root)
        .def("set_whisker_spray_root_fn_gate", &Resistance::set_whisker_spray_root_fn_gate,
             py::arg("lo"), py::arg("hi"))
        .def("set_whisker_spray_turb_warp", &Resistance::set_whisker_spray_turb_warp)
        .def("get_whisker_spray_turb_warp", &Resistance::get_whisker_spray_turb_warp)
        .def("set_wave_lift_volume_scale", &Resistance::set_wave_lift_volume_scale)
        .def("get_wave_lift_volume_scale", &Resistance::get_wave_lift_volume_scale)
        .def("set_sectional_crossflow", &Resistance::set_sectional_crossflow)
        .def("get_sectional_crossflow", &Resistance::get_sectional_crossflow)
        .def("set_whisker_local_beta", &Resistance::set_whisker_local_beta, py::arg("on"))
        .def("get_whisker_local_beta", &Resistance::get_whisker_local_beta)
        .def("set_planing_friction_running_trim", &Resistance::set_planing_friction_running_trim)
        .def("set_friction_true_wetted", &Resistance::set_friction_true_wetted)
        .def("set_friction_deadrise_scaled", &Resistance::set_friction_deadrise_scaled)
        .def("set_friction_cv_ramp", &Resistance::set_friction_cv_ramp,
             py::arg("lo"), py::arg("hi"))
        .def("set_planing_ca", &Resistance::set_planing_ca)
        .def("set_localfs_wave_envelope", &Resistance::set_localfs_wave_envelope)
        .def("get_localflow_drag", &Resistance::get_localflow_drag)
        .def("get_localflow_lift", &Resistance::get_localflow_lift)
        .def("set_sectional_attitude", &Resistance::set_sectional_attitude,
             py::arg("on"), py::arg("K"), py::arg("warp_floor") = -1.0)
        .def("set_sectional_delta", &Resistance::set_sectional_delta,
             py::arg("on"), py::arg("K") = -1.0, py::arg("defl_lo") = -1.0, py::arg("defl_hi") = -1.0,
             py::arg("fn_lo") = -1.0, py::arg("fn_hi") = -1.0)
        .def("get_sectional_deficit", &Resistance::get_sectional_deficit, py::arg("speed"))
        .def("set_sectional_tloc_cap", &Resistance::set_sectional_tloc_cap, py::arg("deg"))
        .def("set_sectional_transom_fix", &Resistance::set_sectional_transom_fix, py::arg("on"))
        .def("set_sectional_ventilation", &Resistance::set_sectional_ventilation, py::arg("on"))
        .def("set_sectional_vent_beam_frac", &Resistance::set_sectional_vent_beam_frac, py::arg("f"))
        .def("set_sectional_similarity_model", &Resistance::set_sectional_similarity_model, py::arg("on"))
        .def("set_sectional_forward_term", &Resistance::set_sectional_forward_term,
             py::arg("on"), py::arg("cs")=0.0)
        .def("debug_sectional_at", &Resistance::debug_sectional_at,
             py::arg("heave"), py::arg("pitch"), py::arg("speed"))
        .def("set_attitude_wave_blend", &Resistance::set_attitude_wave_blend, py::arg("on"))
        .def("get_attitude_wave_blend", &Resistance::get_attitude_wave_blend)
        .def("get_whisker_spray_drag", &Resistance::get_whisker_spray_drag)
        .def("set_air_resistance", &Resistance::set_air_resistance)
        .def("get_air_resistance_on", &Resistance::get_air_resistance_on)
        .def("set_air_drag_coef", &Resistance::set_air_drag_coef)
        .def("get_air_drag_coef", &Resistance::get_air_drag_coef)
        .def("get_air_resistance", &Resistance::get_air_resistance)
        .def("set_mercier_savitsky", &Resistance::set_mercier_savitsky)
        .def("set_radojcic_nss", &Resistance::set_radojcic_nss)
        .def("set_radojcic_total_envelope", &Resistance::set_radojcic_total_envelope)
        .def("set_radojcic_residuary_beta_gate", &Resistance::set_radojcic_residuary_beta_gate)
        .def("set_sectional_nx", &Resistance::set_sectional_nx)
        .def("set_sectional_smooth_frac", &Resistance::set_sectional_smooth_frac)
        .def("get_radojcic_nss_residuary", &Resistance::get_radojcic_nss_residuary)
        .def("get_radojcic_nss_total", &Resistance::get_radojcic_nss_total)
        .def("get_mercier_savitsky_residuary", &Resistance::get_mercier_savitsky_residuary)
        .def("get_propeller_wake_field",
             [](Resistance& r, double speed, py::tuple center,
                double inner_diameter, double outer_diameter,
                int n_radial, int n_azimuthal, bool include_localflow_axial) {
                 glm::dvec3 c(center[0].cast<double>(),
                              center[1].cast<double>(),
                              center[2].cast<double>());
                 return r.get_propeller_wake_field(speed, c, inner_diameter, outer_diameter,
                                                   n_radial, n_azimuthal,
                                                   include_localflow_axial);
             },
             py::arg("speed"), py::arg("center"),
             py::arg("inner_diameter"), py::arg("outer_diameter"),
             py::arg("n_radial") = 12, py::arg("n_azimuthal") = 36,
             py::arg("include_localflow_axial") = true)
        .def("find_equilibrium", &Resistance::find_equilibrium)
        .def("set_itp_root", &Resistance::set_itp_root, py::arg("on"))
        .def("get_itp_root", &Resistance::get_itp_root)
        .def("debug_force_balance", &Resistance::debug_force_balance)
        .def("residual_at", &Resistance::residual_at,
             py::arg("heave"), py::arg("pitch"), py::arg("speed"))
        .def("set_grid_precision", &Resistance::set_grid_precision)
        .def("set_transom_method", &Resistance::set_transom_method)
        .def("get_michell_drag_lift_torque", &Resistance::get_michell_drag_lift_torque,
             py::arg("speed"))
        .def("get_nearfar_drag_ratio", &Resistance::get_nearfar_drag_ratio, py::arg("speed"))
    ;

    py::class_<ProfileGrid>(resistance, "ProfileGrid")
        .def(py::init<std::shared_ptr<Hull>>())
        .def("size_x", &ProfileGrid::size_x)
        .def("size_z", &ProfileGrid::size_z)
        .def("get_dx", &ProfileGrid::get_dx)
        .def("get_dz", &ProfileGrid::get_dz)
        .def("setup", &ProfileGrid::setup)
        .def("get_xs", &ProfileGrid::get_xs)
        .def("get_zs", &ProfileGrid::get_zs)
        .def("get_beams", &ProfileGrid::get_beams)
        .def("get_slopes", &ProfileGrid::get_slopes)
        .def("get_slopes0", &ProfileGrid::get_slopes0)
    ;
    py::enum_<PropellerType>(resistance, "PropellerType")
        .value("BSeries", PropellerType::BSeries)
        .value("Gawn", PropellerType::Gawn)
        .value("KCA", PropellerType::KCA)
        .value("SurfacePiercing", PropellerType::SurfacePiercing)
        .export_values()
    ;

    py::class_<Propeller, std::shared_ptr<Propeller>>(resistance, "Propeller")
        .def_static("get_version", &Propeller::get_version)
        .def(py::init<std::shared_ptr<Environment>>())
        .def_property("type", &Propeller::get_type, &Propeller::set_type)
        .def_property("rpm", &Propeller::get_rpm, &Propeller::set_rpm)
        .def_property("pitch", &Propeller::get_pitch, &Propeller::set_pitch)
        .def_property("blades", &Propeller::get_blades, &Propeller::set_blades)
        .def_property("diameter", &Propeller::get_diameter, &Propeller::set_diameter)
        .def_property("area_ratio", &Propeller::get_area_ratio, &Propeller::set_area_ratio)
        .def_property("shaft_draft", &Propeller::get_shaft_draft, &Propeller::set_shaft_draft)
        .def_property("inclination", &Propeller::get_inclination, &Propeller::set_inclination)
        .def_property("advance_speed", &Propeller::get_advance_speed, &Propeller::set_advance_speed)

        .def("set_rpm_by_j", &Propeller::set_rpm_by_j, py::arg("J"))
        .def("get_rps", &Propeller::get_rps)
        .def("get_advance_coef", &Propeller::get_advance_coef)
        .def("get_thrust_coef", &Propeller::get_thrust_coef)
        .def("get_torque_coef", &Propeller::get_torque_coef)
        .def("get_thrust", &Propeller::get_thrust)
        .def("get_torque", &Propeller::get_torque)
        .def("get_efficiency", &Propeller::get_efficiency)
        .def("get_delivered_power", &Propeller::get_delivered_power)
        .def("get_cavitation_number", &Propeller::get_cavitation_number)
        .def("get_cavitating_area_ratio", &Propeller::get_cavitating_area_ratio, py::arg("thrust"))
        .def("get_cavitation", &Propeller::get_cavitation)
    ;
}
