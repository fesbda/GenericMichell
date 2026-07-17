#ifndef PROPELLER_H
#define PROPELLER_H

#include "environment.h"
//#include "tools.h"
//#include "interpolation.h"
#include <array>

/// Implemented propeller series
enum class PropellerType
{
    BSeries,
    Gawn,
    KCA,
    SurfacePiercing //(3<=Z<=6, 0.5<=EAR<=0.8, 0.5<=P/D<=2, 0.1<=J<=2)
};

/// Propeller calculating Wageningen, Gawn, SPP & KCA series
class Propeller
{

private:

    /// Get minimum number of blades for a propeller series
    static int min_blades(PropellerType type)
    {
        switch (type) {
        case PropellerType::BSeries:
            return 2;
        case PropellerType::Gawn:
            return 3;
        case PropellerType::KCA:
            return 3;
        case PropellerType::SurfacePiercing:
            return 3;
        }
        return 0;
    }

    /// Get maximum number of blades for a propeller series
    static int max_blades(PropellerType type)
    {
        switch (type) {
        case PropellerType::BSeries:
            return 7;
        case PropellerType::Gawn:
            return 3;
        case PropellerType::KCA:
            return 3;
        case PropellerType::SurfacePiercing:
            return 6;
        }
        return 0;
    }

    static double min_area(PropellerType type)
    {
        switch (type) {
        case PropellerType::BSeries:
            return 0.3;
        case PropellerType::Gawn:
            return 0.2;
        case PropellerType::KCA:
            return 0.5;
        case PropellerType::SurfacePiercing:
            return 0.5;
        }
        return 0;
    }

    static double max_area(PropellerType type)
    {
        switch (type) {
        case PropellerType::BSeries:
            return 1.05;
        case PropellerType::Gawn:
            return 1.1;
        case PropellerType::KCA:
            return 1.1;
        case PropellerType::SurfacePiercing:
            return 0.8;
        }
        return 0;
    }

    static double min_pd(PropellerType type)
    {
        switch (type) {
        case PropellerType::BSeries:
            return 0.45;
        case PropellerType::Gawn:
            return 0.8;
        case PropellerType::KCA:
            return 0.8;
        case PropellerType::SurfacePiercing:
            return 0.5;
        }
        return 0;
    }

    static double max_pd(PropellerType type)
    {
        switch (type) {
        case PropellerType::BSeries:
            return 1.4;
        case PropellerType::Gawn:
            return 1.4;
        case PropellerType::KCA:
            return 1.8;
        case PropellerType::SurfacePiercing:
            return 2.0;
        }
        return 0;
    }

    // Thrust polynomial coefficients for B-series
    constexpr static std::size_t NT = 39;
    constexpr static std::array<double, NT> BSeries_CT = {0.00880496, -0.204554, 0.166351, 0.158114, -0.147581, -0.481497, 0.415437, 0.0144043, -0.0530054, 0.0143481, 0.0606826, -0.0125894, 0.0109689,
        -0.133698, 0.00638407, -0.00132718, 0.168496, -0.0507214, 0.0854559, -0.0504475, 0.010465, -0.00648272, -0.00841728, 0.0168424, -0.00102296, -0.0317791, 0.018604, -0.00410798, -0.000606848,
        -0.0049819, 0.0025983, -0.000560528, -0.00163652, -0.000328787, 0.000116502, 0.000690904, 0.00421749, 5.65229E-05, -0.00146564};
    constexpr static std::array<int, NT> BSeries_ST = {0, 1, 0, 0, 2, 1, 0, 0, 2, 0, 1, 0, 1, 0, 0, 2, 3, 0, 2, 3, 1, 2, 0, 1, 3, 0, 1, 0, 0, 1, 2, 3, 1, 1, 2, 0, 0, 3, 0};
    constexpr static std::array<int, NT> BSeries_TT = {0, 0, 1, 2, 0, 1, 2, 0, 0, 1, 1, 0, 0, 3, 6, 6, 0, 0, 0, 0, 6, 6, 3, 3, 3, 3, 0, 2, 0, 0, 0, 0, 2, 6, 6, 0, 3, 6, 3};
    constexpr static std::array<int, NT> BSeries_UT = {0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 2, 2, 2, 2, 2, 0, 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2};
    constexpr static std::array<int, NT> BSeries_VT = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

    // Torque polynomial coefficients for B-series
    constexpr static std::size_t NQ = 47;
    constexpr static std::array<double, NQ> BSeries_CQ = {0.00379368, 0.00886523, -0.032241, 0.00344778, -0.0408811, -0.108009, -0.0885381, 0.188561, -0.00370871, 0.00513696, 0.0209449, 0.00474319,
        -0.00723408, 0.00438388, -0.0269403, 0.0558082, 0.0161886, 0.00318086, 0.015896, 0.0471729, 0.0196283, -0.0502782, -0.030055, 0.0417122, -0.0397722, -0.00350024, -0.0106854, 0.0010903,
        -0.000313912, 0.0035985, -0.00142121, -0.00383637, 0.0126803, -0.00318278, 0.00334268, -0.00183491, 0.000112451, -2.97228E-05, 0.000269551, 0.00083265, 0.00155334, 0.000302683, -0.0001843,
        -0.000425399, 8.69243E-05, -0.0004659, 5.54194E-05};
    constexpr static std::array<int, NQ> BSeries_SQ = {0, 2, 1, 0, 0, 1, 2, 0, 1, 0, 1, 2, 2, 1, 0, 3, 0, 1, 0, 1, 3, 0, 3, 2, 0, 0, 3, 3, 0, 3, 0, 1, 0, 2, 0, 1, 3, 3, 1, 2, 0, 0, 0, 0, 3, 0, 1};
    constexpr static std::array<int, NQ> BSeries_TQ = {0, 0, 1, 2, 1, 1, 1, 2, 0, 1, 1, 1, 0, 1, 2, 0, 3, 3, 0, 0, 0, 1, 1, 2, 3, 6, 0, 3, 6, 0, 6, 0, 2, 3, 6, 1, 2, 6, 0, 0, 2, 6, 0, 3, 3, 6, 6};
    constexpr static std::array<int, NQ> BSeries_UQ = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 1, 1, 2, 2, 2, 2, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2};
    constexpr static std::array<int, NQ> BSeries_VQ = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

    // Thrust polynomial coefficients for Gawn series
    constexpr static std::array<double, NT> Gawn_CT = {-0.05586363, -0.21730109, 0.2605314, 0.158114, -0.147581, -0.481497, 0.37812278, 0.0144043, -0.0530054, 0.0143481, 0.0606826, -0.0125894,
        0.0109689, -0.133698, 0.00241157, -0.00053002, 0.168496, 0.02634542, 0.04360136, -0.03118493, 0.01249215, -0.00648272, -0.00841728, 0.0168424, -0.00102296, -0.0317791, 0.018604, -0.00410798,
        -0.000606848, -0.0049819, 0.0025963, -0.000560528, -0.00163652, -0.000328787, 0.000116502, 0.000690904, 0.00421749, 5.65229E-05, -0.00146564};
    constexpr static std::array<int, NT> Gawn_ST = {0, 1, 0, 0, 2, 1, 0, 0, 2, 0, 1, 0, 1, 0, 0, 2, 3, 0, 2, 3, 1, 2, 0, 1, 3, 0, 1, 0, 0, 1, 2, 3, 1, 1, 2, 0, 0, 3, 0};
    constexpr static std::array<int, NT> Gawn_TT = {0, 0, 1, 2, 0, 1, 2, 0, 0, 1, 1, 0, 0, 3, 6, 6, 0, 0, 0, 0, 6, 6, 3, 3, 3, 3, 0, 2, 0, 0, 0, 0, 2, 6, 6, 0, 3, 6, 3};
    constexpr static std::array<int, NT> Gawn_UT = {0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 2, 2, 2, 2, 2, 0, 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2};
    constexpr static std::array<int, NT> Gawn_VT = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

    // Torque polynomial coefficients for B-series
    constexpr static std::array<double, NQ> Gawn_CQ = {0.00515898, 0.01606668, -0.0441153, 0.00682223, -0.0408811, -0.07732967, -0.0885381, 0.16937502, -0.00370871, 0.00513696, 0.0209449, 0.00474319,
        -0.00723408, 0.00438388, -0.0269403, 0.0558082, 0.0161886, 0.00318086, 0.01290435, 0.02445084, 0.00700643, -0.02719046, -0.01664586, 0.0300449, -0.03369749, -0.00350024, -0.0106854,
        0.00110903, -0.000313912, 0.0035895, -0.00142121, -0.00383637, 0.0126803, -0.00318278, 0.00334268, -0.00183491, 0.000112451, -2.97228E-05, 0.000269551, 0.00083265, 0.00155334, 0.000302683,
        -0.0001843, -0.000425399, 8.69243E-05, -0.0004659, 5.54194E-05};
    constexpr static std::array<int, NQ> Gawn_SQ = {0, 2, 1, 0, 0, 1, 2, 0, 1, 0, 1, 2, 2, 1, 0, 3, 0, 1, 0, 1, 3, 0, 3, 2, 0, 0, 3, 3, 0, 3, 0, 1, 0, 2, 0, 1, 3, 3, 1, 2, 0, 0, 0, 0, 3, 0, 1};
    constexpr static std::array<int, NQ> Gawn_TQ = {0, 0, 1, 2, 1, 1, 1, 2, 0, 1, 1, 1, 0, 1, 2, 0, 3, 3, 0, 0, 0, 1, 1, 2, 3, 6, 0, 3, 6, 0, 6, 0, 2, 3, 6, 1, 2, 6, 0, 0, 2, 6, 0, 3, 3, 6, 6};
    constexpr static std::array<int, NQ> Gawn_UQ = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 1, 1, 2, 2, 2, 2, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2};
    constexpr static std::array<int, NQ> Gawn_VQ = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

    // Thrust polynomial coefficients for KCA series
    constexpr static std::size_t KCA_NT = 16;
    constexpr static std::array<double, KCA_NT> KCA_CT
        = {0.1193852, -0.6574682, 0.3493294, 0.4119366, -0.1991927, 5.8630510, -1.1077350, -0.1341679, 0.2628839, -0.5217023, 0.2970728, 6.1525800, -2.4708400, -4.0801660, 4.1542010, -1.1364520};
    constexpr static std::array<int, KCA_NT> KCA_ET = {0, 0, 0, 0, 0, -2, -2, 0, 0, 0, 0, -2, -2, -3, -3, -3};
    constexpr static std::array<int, KCA_NT> KCA_XT = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 1, 1, 2};
    constexpr static std::array<int, KCA_NT> KCA_YT = {0, 0, 1, 1, 2, 2, 2, 0, 0, 1, 2, 1, 2, 6, 6, 6};
    constexpr static std::array<int, KCA_NT> KCA_ZT = {0, 1, 0, 1, 1, 2, 3, 0, 1, 1, 0, 3, 3, 0, 1, 0};
    // Thrust cavitation for KCA series
    constexpr static std::size_t KCA_NTC = 20;
    constexpr static std::array<double, KCA_NTC> KCA_DTC = {6.688144, 3.579195, -5.700350, -1.359994, -8.111903, 4.770548, -2.313208, -1.387858, 4.992201, -7.161204, 1.721436, 2.322218, -1.156897,
        5.014178, -6.555364, 2.852867, -8.081759, 8.671852, -3.727835, 8.043970};
    constexpr static std::array<int, KCA_NTC> KCA_ETC = {-2, 0, 0, 0, 0, 1, -1, 1, 1, 1, 1, 1, 2, -2, -2, -1, -1, 1, 1, 1};
    constexpr static std::array<int, KCA_NTC> KCA_STC = {0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 0, 0, 1, 1, 3, 3, 3};
    constexpr static std::array<int, KCA_NTC> KCA_TTC = {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 0, 0, 0, 1, 2, 3, 3};
    constexpr static std::array<int, KCA_NTC> KCA_UTC = {0, 2, 3, 0, 1, 3, 0, 2, 3, 4, 0, 1, 2, 0, 1, 1, 0, 5, 0, 1};
    constexpr static std::array<int, KCA_NTC> KCA_VTC = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 1, 0, 0, 0};

    // Torque polynomial coefficients for KCA series
    constexpr static std::size_t KCA_NQ = 17;
    constexpr static std::array<double, KCA_NQ> KCA_CQ = {1.5411660, 0.1091688, -0.3102420, 0.1547428, -4.3706150, 0.2490295, -0.1594602, 8.5367470, -9.5121630, -9.3203070, 3.2878050, 5.4960340,
        -4.8650630, -0.1062500, 8.5299550, 1.1010230, -3.1517560};
    constexpr static std::array<int, KCA_NQ> KCA_EQ = {-3, 0, 0, 0, -2, 0, 0, -2, -2, -3, -2, -2, -2, 0, -2, -2, -3};
    constexpr static std::array<int, KCA_NQ> KCA_XQ = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2};
    constexpr static std::array<int, KCA_NQ> KCA_YQ = {0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 0, 1, 1, 2, 0, 2};
    constexpr static std::array<int, KCA_NQ> KCA_ZQ = {0, 1, 2, 3, 0, 2, 3, 0, 1, 2, 3, 1, 0, 1, 0, 3, 2};
    // Torque cavitation for KCA series
    constexpr static std::size_t KCA_NQC = 18;
    constexpr static std::array<double, KCA_NQC> KCA_DQC
        = {4.024475, 1.202447, -9.836070, -8.318840, 5.098177, -5.192839, 2.641109, -1.688934, 4.928417, 1.024274, -1.194521, 5.498736, -2.488235, -5.832879, 1.503955, -3.316121, 3.890792, 1.682032};
    constexpr static std::array<int, KCA_NQC> KCA_EQC = {-3, -1, -2, -1, 0, -1, 0, 1, -2, -2, -1, -2, -1, -1, -1, 0, 0, 1};
    constexpr static std::array<int, KCA_NQC> KCA_SQC = {0, 0, 1, 1, 1, 2, 2, 2, 0, 0, 0, 1, 1, 0, 0, 3, 3, 3};
    constexpr static std::array<int, KCA_NQC> KCA_TQC = {0, 0, 1, 1, 1, 1, 2, 2, 0, 0, 1, 0, 1, 0, 3, 3, 3, 3};
    constexpr static std::array<int, KCA_NQC> KCA_UQC = {0, 2, 0, 1, 3, 1, 0, 3, 1, 0, 1, 1, 0, 5, 0, 0, 1, 3};
    constexpr static std::array<int, KCA_NQC> KCA_VQC = {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 1, 1, 0, 0, 0, 0, 0};

    // prop
    std::shared_ptr<Environment> env;
    int n_blades = 3;
    double diameter = 0.4;
    double pitch = 0.4;
    double area_ratio = 0.6;
    double advance_speed = 1.0;
    double rps = 10.0;
    double shaft_draft = 1.0;
    double inclination = 0.0;
    PropellerType type = PropellerType::BSeries;

public:
    Propeller(std::shared_ptr<Environment> environment)
        : env(environment)
    { }
    ~Propeller() { }

    /// Get the version of the solver
    static std::string get_version() { return std::string(__DATE__); }

    /// Get environment properties
    std::shared_ptr<Environment> get_env() const { return env; }

    // GEOMETRY:

    /// Set the geometry/series of the prop
    void set_type(PropellerType value) { type = value; }

    /// Get the geometry/series of the prop
    PropellerType get_type() const { return type; }

    /// Set number of prop blades, Z_min <= Z <= Z_max
    void set_blades(int value) { n_blades = value; }

    /// Get number of prop blades, Z
    int get_blades() const { return n_blades; }

    /// Set the diameter of the prop, D
    void set_diameter(double value) { diameter = value; }

    /// Get the diameter of the prop, D
    double get_diameter() const { return diameter; }

    /// Set the average pitch of the prop, P
    void set_pitch(double value) { pitch = value; }

    /// Get the average pitch of the prop, P
    double get_pitch() const { return pitch; }

    /// Set the ratio of expanded surface area and area of disc, A_E / A_0
    void set_area_ratio(double value) { area_ratio = value; }

    /// Get the ratio of expanded surface area and area of disc, A_E / A_0
    double get_area_ratio() const { return area_ratio; }

    /// Get the developed area ratio, A_D / A_0, from expanded surface area
    double get_developed_area_ratio() const { return (std::sqrt(3200 * get_area_ratio() * get_blades() / 17.0 + 121 * sq(get_blades())) - 11 * get_blades()) / 8.0; }

    /// Set the immersed depth of the shaft/hub
    void set_shaft_draft(double value) { shaft_draft = value; }

    /// Get the immersed depth of the shaft/hub
    double get_shaft_draft() const { return shaft_draft; }

    /// Set shaft inclination in degrees
    void set_inclination(double value) { inclination = radians(value); }

    /// Get shaft inclination in degrees
    double get_inclination() const { return degrees(inclination); }

    // OPERATING POINT:

    /// Set the advancing speed of the prop in open water
    void set_advance_speed(double value) { advance_speed = value; }

    /// Get the advancing speed of the prop in open water
    double get_advance_speed() const { return advance_speed; }

    /// Set rotation speed in RPM
    void set_rpm(double value) { rps = value / 60.0; }

    /// Set the rotation speed RPM by inputting advancing coefficient J
    void set_rpm_by_j(double J) { rps = advance_speed / (J * diameter); }

    /// Get rotation speed in RPM
    double get_rpm() const { return rps * 60.0; }

    /// Get rotation speed in RPS
    double get_rps() const { return rps; }

    /// Get the advance coefficient, J
    double get_advance_coef() const { return get_advance_speed() / (rps * diameter); }

    /// Get the thrust coefficient delivered by the prop, K_T
    double get_thrust_coef() const
    {
        double const J = get_advance_coef();
        double const PD = (pitch / diameter);

        if (n_blades < min_blades(type) || n_blades > max_blades(type) || J < 0.3 || J > 1.8 || PD < min_pd(type) || PD > max_pd(type) || area_ratio < min_area(type) || area_ratio > max_area(type)) {
            return 0.0;
        }

        if (type == PropellerType::SurfacePiercing) {

            // A Review on the Hydrodynamic Characteristics of the SPP Concerning to the Available Experimental Data and Evaluating Regression Polynomial Functions
            // Montazeri and Ghassemi
            double const K_T = -0.133 - 0.286 * J + 0.507 * PD + 0.197 * cub(PD) + 0.164 * cub(J) + 0.564 * PD * J * n_blades - 0.488 * sq(PD) * J - 0.176 * sq(PD) * n_blades - 0.28 * sq(J) * n_blades
                - 0.601 * sq(area_ratio) * J - 0.021 * sq(n_blades) * PD + 0.175 * sq(area_ratio) * n_blades;
            return std::max(0.0, K_T);

            // Ferrando et al.
            // double K_T5 = -0.61986*J + 0.14553*PD + 0.72956*J*PD - 0.3049*sq(J) - 0.12523*sq(PD) + 0.28459;
            // double K_T4 = -0.691625*J + 0.794973*PD + 0.870696*J*PD - 0.395012*sq(J) - 0.515183*sq(PD);
            // return glm::mix(K_T4, K_T5, n_blades - 4);

        } else if (type == PropellerType::BSeries) {

            double K_T = 0.0;
            for (std::size_t i = 0; i < NT; i++) {
                K_T += BSeries_CT[i] * std::pow(J, BSeries_ST[i]) * std::pow(PD, BSeries_TT[i]) * std::pow(area_ratio, BSeries_UT[i]) * std::pow(n_blades, BSeries_VT[i]);
            }

            // Reynolds number fix
            double const LR = std::log10(glm::clamp(get_section_rn(0.75), 1e6, 1e9)) - 0.301;
            double const delta_K_T = 0.0003534850 - 0.0033375800 * area_ratio * sq(J) - 0.0047812500 * area_ratio * PD * J + 0.0002577920 * sq(LR) * area_ratio * sq(J)
                + 0.0000643192 * LR * pown<6>(PD) * sq(J) - 0.0000110636 * sq(LR) * pown<6>(PD) * sq(J) - 0.0000276305 * sq(LR) * n_blades * area_ratio * sq(J)
                + 0.0000954000 * LR * n_blades * area_ratio * PD * J + 0.0000032049 * LR * sq(n_blades) * area_ratio * cub(PD) * J;
            K_T += delta_K_T;

            return std::max(0.0, K_T);

        } else if (type == PropellerType::Gawn) {

            double K_T = 0.0;
            for (std::size_t i = 0; i < NT; i++) { K_T += Gawn_CT[i] * std::pow(J, Gawn_ST[i]) * std::pow(PD, Gawn_TT[i]) * std::pow(area_ratio, Gawn_UT[i]) * std::pow(n_blades, Gawn_VT[i]); }
            return std::max(0.0, K_T);

        } else if (type == PropellerType::KCA) {

            double const AdA0 = get_developed_area_ratio();
            double K_T = 0.0;
            for (std::size_t i = 0; i < KCA_NT; i++) { K_T += KCA_CT[i] * std::pow(10.0, KCA_ET[i]) * std::pow(AdA0, KCA_XT[i]) * std::pow(PD, KCA_YT[i]) * std::pow(J, KCA_ZT[i]); }

            // fix thrust for cavitation
            double dK_T = 0.0;
            double const T0 = K_T * env->get_density() * sq(rps) * sq(sq(diameter));
            if (get_cavitation(T0) > 0.0) {
                double const omega = get_cavitation_number();
                for (std::size_t i = 0; i < KCA_NTC; i++) {
                    dK_T += KCA_DTC[i] * std::pow(10.0, KCA_ETC[i]) * std::pow(AdA0, KCA_STC[i]) * std::pow(omega, KCA_TTC[i]) * std::pow(K_T, KCA_UTC[i]) * std::pow(PD, KCA_VTC[i]);
                }
                dK_T = std::max(dK_T, 0.0);
            }

            //std::cout << "KT " << K_T << " + dKT " << dK_T << std::endl;
            // FIXME(KCA cavitation sign — SOLVER_FIX_PLAN.md H14): as written, cavitation can only ADD to
            // thrust (dK_T>=0, K_T+dK_T), which is backwards — cavitation causes thrust breakdown, so the
            // net effect should be a REDUCTION. The KCA_DTC regression depends on the open-water K_T
            // (pow(K_T, KCA_UTC[i])), which is typical of a formula that returns the CAVITATING KTC directly
            // rather than a delta; in that case the correct form is `return max(0, dK_T)` (replace), not add.
            // The exact convention needs Radojcic (1988) "Regression analysis of the Gawn-Burrill KCA series".
            // Standalone module: does NOT feed the resistance/planing solver, so this does not affect any
            // validation number. Left as-is pending the primary source rather than flip an unverified sign.
            return std::max(0.0, K_T + dK_T);
        }

        return 0.0;
    }

    /// Get the thrust force delivered by the prop, T
    double get_thrust() const
    {
        double area = sq(diameter);

        // SPPs are partially immersed
        if (type == PropellerType::SurfacePiercing) {
            double const It = (shaft_draft + diameter * 0.5) / diameter;
            area *= 0.25 * (M_PI_2 - std::asin(glm::clamp((0.5 - It) / 0.5, -1.0, 1.0)));
        }

        // TODO make thrust smaller due to (partial) cavitation

        return std::cos(inclination) * get_thrust_coef() * env->get_density() * sq(rps) * sq(diameter) * area;
    }

    /// Get the torque coefficient delivered by the prop, K_Q
    double get_torque_coef() const
    {
        double const J = get_advance_coef();
        double const PD = (pitch / diameter);

        if (n_blades < min_blades(type) || n_blades > max_blades(type) || J < 0.3 || J > 1.8 || PD < min_pd(type) || PD > max_pd(type) || area_ratio < min_area(type) || area_ratio > max_area(type)) {
            return 0.0;
        }

        if (type == PropellerType::SurfacePiercing) {

            // A Review on the Hydrodynamic Characteristics of the SPP Concerning to the Available Experimental Data and Evaluating Regression Polynomial Functions
            // Montazeri and Ghassemi
            double const K_Q_10 = 0.12 - 1.12 * J + 0.12 * cub(PD) + 0.34 * cub(J) + 0.58 * PD * J * n_blades + 0.07 * sq(n_blades) * J - 0.08 * sq(n_blades) * PD + 0.68 * sq(area_ratio) * PD
                - 0.48 * sq(J) * n_blades - 0.45 * sq(PD) * area_ratio - 0.44 * sq(PD) * J + 1.02 * PD;
            return std::max(0.0, K_Q_10 * 0.1);

            // Ferrando et al.
            //...

        } else if (type == PropellerType::BSeries) {

            double K_Q = 0.0;
            for (std::size_t i = 0; i < NQ; i++) {
                K_Q += BSeries_CQ[i] * std::pow(J, BSeries_SQ[i]) * std::pow(PD, BSeries_TQ[i]) * std::pow(area_ratio, BSeries_UQ[i]) * std::pow(n_blades, BSeries_VQ[i]);
            }

            // Reynolds number fix
            double const LR = std::log10(glm::clamp(get_section_rn(0.75), 1e6, 1e9)) - 0.301;
            double const delta_K_Q = -0.00059141200 + 0.00696898000 * PD - 0.00006666540 * n_blades * pown<6>(PD) + 0.01608180000 * sq(area_ratio) - 0.00093809100 * LR * PD
                - 0.00059593000 * LR * sq(PD) + 0.00007820990 * sq(LR) * sq(PD) + 0.00000521990 * LR * n_blades * area_ratio * sq(J) - 0.00000088528 * sq(LR) * n_blades * area_ratio * PD * J
                + 0.00002301710 * LR * n_blades * pown<6>(PD) - 0.00000184341 * sq(LR) * n_blades * pown<6>(PD) - 0.00400252000 * LR * sq(area_ratio) + 0.00022091500 * sq(LR) * sq(area_ratio);
            K_Q += delta_K_Q;

            return std::max(0.0, K_Q);

        } else if (type == PropellerType::Gawn) {

            double K_Q = 0.0;
            for (std::size_t i = 0; i < NQ; i++) { K_Q += Gawn_CQ[i] * std::pow(J, Gawn_SQ[i]) * std::pow(PD, Gawn_TQ[i]) * std::pow(area_ratio, Gawn_UQ[i]) * std::pow(n_blades, Gawn_VQ[i]); }
            return std::max(0.0, K_Q);

        } else if (type == PropellerType::KCA) {

            double const AdA0 = get_developed_area_ratio();
            double K_Q = 0.0;

            // non-cavitating torque
            for (std::size_t i = 0; i < KCA_NQ; i++) { K_Q += KCA_CQ[i] * std::pow(10.0, KCA_EQ[i]) * std::pow(AdA0, KCA_XQ[i]) * std::pow(PD, KCA_YQ[i]) * std::pow(J, KCA_ZQ[i]); }

            // non-cavitating thrust, needed to fix K_Q for cavitation
            double K_T = 0.0;
            for (std::size_t i = 0; i < KCA_NT; i++) { K_T += KCA_CT[i] * std::pow(10.0, KCA_ET[i]) * std::pow(AdA0, KCA_XT[i]) * std::pow(PD, KCA_YT[i]) * std::pow(J, KCA_ZT[i]); }
            K_T = std::max(0.0, K_T);
            double const T0 = K_T * env->get_density() * sq(rps) * sq(sq(diameter));

            // fix torque for cavitation
            double dK_Q = 0.0;
            if (get_cavitation(T0) > 0.0) {
                double const omega = get_cavitation_number();
                for (std::size_t i = 0; i < KCA_NQC; i++) {
                    dK_Q += KCA_DQC[i] * std::pow(10.0, KCA_EQC[i]) * std::pow(AdA0, KCA_SQC[i]) * std::pow(omega, KCA_TQC[i]) * std::pow(K_T, KCA_UQC[i]) * std::pow(PD, KCA_VQC[i]);
                }
                dK_Q = std::max(dK_Q, 0.0);
            }

            //std::cout << "KQ " << K_Q << " + dKQ " << dK_Q << std::endl;
            // FIXME(KCA cavitation sign — SOLVER_FIX_PLAN.md H14): same suspect convention as get_thrust's
            // cavitation branch (dK_Q>=0 then added). Resolve add-vs-replace and sign against Radojcic (1988)
            // before trusting cavitating torque. Standalone module: no effect on the resistance solver.
            return std::max(0.0, K_Q + dK_Q);
        }

        return 0.0;
    }

    /// Get the torque delivered by the prop, Q
    double get_torque() const
    {
        double area = sq(diameter);

        if (type == PropellerType::SurfacePiercing) {
            double const It = (shaft_draft + diameter * 0.5) / diameter;
            area *= 0.25 * (M_PI_2 - std::asin(glm::clamp((0.5 - It) / 0.5, -1.0, 1.0)));
        }

        return get_torque_coef() * env->get_density() * sq(rps) * cub(diameter) * area;
    }

    /// Get the power delivered by the prop, P_D
    double get_delivered_power() const { return 2 * M_PI * rps * get_torque(); }

    /// Get the efficiency of the prop operating point
    double get_efficiency() const
    {
        double const KQ = get_torque_coef();
        if (KQ < 1e-8) {
            return 0.0;
        }
        return (get_advance_coef() * get_thrust_coef()) / (2 * M_PI * KQ);
    }

    /// Get profile speed for a blade section, as combination of advancing speed and rotation (tangential) speed
    double get_section_speed(double r_ratio) const { return std::sqrt(sq(get_advance_speed()) + sq(r_ratio * M_PI * rps * diameter)); }

    /// Get chord length for a blade section
    double get_section_chord(double r_ratio) const
    {
        double const c70 = 2.737 * std::pow(n_blades, -1.186) * diameter * area_ratio;
        double const c75 = 2.073 * area_ratio * diameter / n_blades;
        return glm::mix(c70, c75, (r_ratio - 0.70) / 0.05);
    }

    /// Get Reynolds number for a blade section
    double get_section_rn(double r_ratio) const { return get_section_speed(r_ratio) * get_section_chord(r_ratio) / env->get_viscosity(); }

    // CAVITATION:

    /// Get cavitation number, sigma
    double get_cavitation_number() const
    {
        double const p_difference = 101325 + get_env()->get_density() * get_env()->get_gravity() * shaft_draft - 2500; // p_ATM + rho*g*h - p_VAPOR
        double const denominator = (0.5 * get_env()->get_density() * sq(get_section_speed(0.7)));
        if (denominator <= 0.0)
            return 0.0;
        return p_difference / denominator;
    }

    /// Keller's equation for minimum AE/A0 to avoid cavitation
    double get_cavitating_area_ratio(double thrust) const
    {
        double const K = 0.1;
        double const rho = get_env()->get_density();
        double const p_V = 2500;
        double const p_0 = 101325 + rho * get_env()->get_gravity() * shaft_draft;
        double const ratio = (1.3 + 0.3 * get_blades()) * thrust / ((p_0 - p_V) * sq(get_diameter())) + K;
        return ratio;
    }

    double get_opt_thickness() const
    {
        // "Further Computer-Analyzed Data of the Wageningen B-Screw Series"
        double const tc_ratio = glm::clamp(0.3 * get_cavitation_number() - 0.012, 0.1, 0.6);
        double const t_opt = tc_ratio * get_section_chord(0.75);

        // Saunders, Hydrodynamics in ship design, Vol. 2, page 620
        double const PS = get_torque() * (rps * 6.283185307179586) / n_blades * 1.341; // PS = shaft horsepower per blade
        double const SC = 11000;                                                       // Sc = maximum allowable stres in PSI
        double const t_min = (0.3048) * (diameter * 3.28084)
            * (0.0028 + 0.21 * std::pow(((2375 - 1125 * pitch / diameter) * PS) / (4.123 * (60 * rps) * cub(diameter * 3.28084) * (SC + sq(diameter * 3.28084 * rps * 60) / 12.788)), 1.0 / 3.0));

        return std::max(t_min, t_opt);
    }

    double get_cavitation(double thrust) const
    {

        if (type == PropellerType::SurfacePiercing) {
            return 1.0; // of course the SPP gonna cavitate
        }

        // TODO add probability, not 1 or 0
        if (get_cavitating_area_ratio(thrust) > area_ratio) {
            return 1.0;
        }
        return 0.0;
    }
};

/// Central class for choosing the propeller
/*
class PropellerChooser
{
private:
    PchipInterpolator power;
    double thrust = 1.0;

public:
    PropellerChooser() { }
    ~PropellerChooser() { }

    void set_power_graph(std::vector<double> const& rpm, std::vector<double> const& P) { power = PchipInterpolator(rpm, P); }

    void set_torque_graph(std::vector<double> const& rpm, std::vector<double> const& Q)
    {

        std::vector<double> P;
        for (std::size_t i = 0; i < Q.size(); i++) { P.push_back(Q[i] * (rpm[i] * 0.10471975511965977)); }
        set_power_graph(rpm, P);
    }

    void set_target_thrust(double value) { thrust = value; }

    void find_optimum(Propeller& prop)
    {

        // vary N, P, D, Ae/Ao
    }
};
*/
#endif // PROPELLER_H
