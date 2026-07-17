#pragma once

/// Properties of the environment
class Environment {

private:
    double gravity{9.80665};
    double density{1025.0};
    double viscosity{1.14e-6};
    double depth{0.0};
    double air_density{1.225};   // kg/m^3, ~15 C at sea level (windage)

public:

    void set_gravity(double value) { gravity = std::abs(value); }
    double get_gravity() const { return gravity; }

    void set_density(double value) { density = std::abs(value); }
    double get_density() const { return  density; }

    /// Air density used for the above-water windage (air-resistance) term.
    void set_air_density(double value) { air_density = std::abs(value); }
    double get_air_density() const { return air_density; }

    void set_viscosity(double value) { viscosity = std::abs(value); }
    double get_viscosity() const { return  viscosity; }

    void set_depth(double value) { depth = std::abs(value); }
    double get_depth() const { return  depth; }
};
