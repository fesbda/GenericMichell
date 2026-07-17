#ifndef STREAMLINE_H
#define STREAMLINE_H

#include "../ext/xfoil/XFoil.h"
#include "hull.h"
#include "interpolation.h"

class Streamline {

private:

    constexpr static int max_input = IBX;
    constexpr static int n_side_panels = 66;
    // typical distribution of points over the coord
    constexpr static std::array<double, n_side_panels> norm_xs = {
        1e-5,         0.0005839, 0.0023342, 0.0052468, 0.0093149, 0.0145291,
        0.0208771, 0.0283441, 0.0369127, 0.0465628, 0.057272,  0.0690152,
        0.0817649, 0.0954915, 0.1101628, 0.1257446, 0.1422005, 0.1594921,
        0.1775789, 0.1964187, 0.2159676, 0.2361799, 0.2570083, 0.2784042,
        0.3003177, 0.3226976, 0.3454915, 0.3686463, 0.3921079, 0.4158215,
        0.4397317, 0.4637826, 0.4879181, 0.5120819, 0.5362174, 0.5602683,
        0.5841786, 0.6078921, 0.6313537, 0.6545085, 0.6773025, 0.6996823,
        0.7215958, 0.7429917, 0.7638202, 0.7840324, 0.8035813, 0.8224211,
        0.8405079, 0.8577995, 0.8742554, 0.8898372, 0.9045085, 0.9182351,
        0.9309849, 0.942728,  0.9534372, 0.9630873, 0.9716559, 0.9791229,
        0.9854709, 0.990685,  0.9947532, 0.9976658, 0.9994161, 1,
    };
    constexpr static std::array<double, n_side_panels> norm_ys = {
        0,         0.0042603, 0.0084289, 0.0125011, 0.0164706, 0.02033,
        0.0240706, 0.0276827, 0.0311559, 0.0344792, 0.0376414, 0.040631,
        0.0434371, 0.0460489, 0.0484567, 0.0506513, 0.0526251, 0.0543715,
        0.0558856, 0.057164,  0.0582048, 0.0590081, 0.0595755, 0.0599102,
        0.0600172, 0.0599028, 0.0595747, 0.0590419, 0.0583145, 0.0574033,
        0.05632,   0.0550769, 0.0536866, 0.052162,  0.0505161, 0.0487619,
        0.0469124, 0.0449802, 0.0429778, 0.0409174, 0.0388109, 0.03667,
        0.0345058, 0.0323294, 0.0301515, 0.0279828, 0.0258337, 0.0237142,
        0.0216347, 0.0196051, 0.0176353, 0.0157351, 0.0139143, 0.0121823,
        0.0105485, 0.0090217, 0.0076108, 0.0063238, 0.0051685, 0.0041519,
        0.0032804, 0.0025595, 0.0019938, 0.001587,  0.0013419, 0.00126};

    std::unique_ptr<XFoil> foil;
    bool one_sided;
    double x_fore;
    double x_length;
    //double vel;
    int fwd;
    double top;

    // results
    std::vector<double> panels_x[2];
    std::vector<double> panels_y[2];
//    std::vector<double> delta_star[2];
    std::vector<double> pressures[2];
    std::vector<double> taus[2];

public:

    void setup_one_side(std::vector<double> xs, std::vector<double> fs, int x_fwd) {

        foil = std::make_unique<XFoil>();

        x_fore = (x_fwd == 1) ? xs.back() : xs.front();
        double const x_aft = (x_fwd == 1) ? xs.front() : xs.back();
        x_length = std::abs(x_fore - x_aft);
        fwd = x_fwd;
        top = *std::max_element(fs.begin(), fs.end());

        LinearInterpolator form(xs, fs);

        double const ratio = std::abs(form(x_fore - 0.3*x_length*fwd) / x_length) / 0.06; // thickness ratio compared to NACA0012

//        std::cout << "Coordinates:" << std::endl;
        std::vector<double> fx;
        std::vector<double> fy;

        // upper surface
        for (int i = n_side_panels - 1; i >= 0; i--) {
            fx.push_back(norm_xs[i]);
            fy.push_back(ratio * norm_ys[i]);
//            std::cout << fx.back() << ";" << fy.back() << std::endl;
        }

        // bottom surface
        for (int i = 1; i < n_side_panels; i++) {
            fx.push_back(norm_xs[i]);
            double const form_y = (form(x_aft + (1.0 - fx.back()) * x_length * x_fwd) - top) / x_length;
            fy.push_back(glm::mix(form_y, -ratio*norm_ys[i], std::exp(-10*fx.back())));
//            std::cout << fx.back() << ";" << fy.back() << std::endl;
        }

        fy[0] += 1e-4;

        bool const success = foil->initXFoilGeometry(fx.size(), fx.data(), fy.data(), nullptr, nullptr);
        if (!success) {
            std::cout << "Failed to initialize streamline geometry" << std::endl;
        }

        one_sided = true;

    }

    /// Speed in m/s
    void calculate(std::shared_ptr<Environment> env, double speed) {

        if (!foil) {
            return;
        }

        double constexpr angle = 0.0;
        double const Rn = speed * x_length / env->get_viscosity();

        for (int i = 0; i < 2; i++) {
            panels_x[i].clear();
            panels_y[i].clear();
//            delta_star[i].clear();
            pressures[i].clear();
            taus[i].clear();
        }

        // init XFoil analysis with Ncrit = 0.0, and force trip on x = 0
        std::stringstream ss;
        foil->initXFoilAnalysis(Rn, angle, 0.0 /*mach*/,
                                0.0 /*ncrit*/, 0.0 /*trip top*/, 0.0 /*trip top*/,
                                1, 1, true /*viscous*/, ss);

        bool success = foil->specal();
        if (!success) {
            std::cout << "Could not converge setup for AoA " << angle << " and Rn " << Rn << std::endl;
        }

        foil->lwake = false;
        foil->lvconv = false;

        int its = 0;
        int constexpr max_its = 3;
        // iterate() returns false (loop continues) forever if viscal() keeps failing; cap the outer loop.
        int outer_guard = 0;
        int constexpr max_outer = 100;
        while (!iterate(its, max_its) && ++outer_guard < max_outer);
        success = foil->lvconv;

//        std::cout << ss.str() << std::endl;

        for (int i_side = 1; i_side <= 2; i_side++) {
//            std::cout << "#SIDE: " << i_side << std::endl;
            for (int i_bl = 2; i_bl <= foil->iblte[i_side]; i_bl++)
            {
                int const i_panel = foil->ipan[i_bl][i_side];

                // scale back to real size
                // TODO rotate to alpha
                panels_x[i_side-1].push_back(foil->x[i_panel]); //x_fore - fwd * foil->x[i_panel] * x_length);
                panels_y[i_side-1].push_back(foil->y[i_panel]); //top + foil->y[i_panel] * x_length);

                if (success) {
//                    delta_star[i_side-1].push_back(foil->dstr[i_bl][i_side] * x_length);
                    pressures[i_side-1].push_back(0.5 * env->get_density() * sq(speed) * foil->cpv[i_panel]);
                    taus[i_side-1].push_back(0.5 * env->get_density() * sq(speed) * foil->tau[i_bl][i_side]);
                } else {
                    pressures[i_side-1].push_back(0.5 * env->get_density() * sq(speed) * foil->cpi[i_panel]);
                    taus[i_side-1].push_back(0);
//                    delta_star[i_side-1].push_back(0);
                }
//                std::cout << panels_x[i_side-1].back() << ", "
//                          << panels_y[i_side-1].back() << ", "
//                          << pressures[i_side-1].back() << ", "
//                          << std::endl;

            }
        }
    }

    std::pair<double,double> get_lift(int side, bool account_friction) {

        int const px = (side == 1) ? -1 : 1;
        side--;
        auto const p = trapz_and_centre(panels_x[side], pressures[side]);
        // NOTE (dormant, account_friction=false today): the friction term's centroid is measured on
        // the panels_y axis but is then combined with the pressure x-centroid and mapped through the
        // x-transform below -- an axis mix. A correct fix needs trapz_and_centre to decouple the force
        // magnitude axis from the location axis; deferred until the friction/BL path is re-enabled.
        auto const f = account_friction ? trapz_and_centre(panels_y[side], taus[side]) : std::make_pair(0.0, 0.0);
        double const total = px*p.first + f.first;
        if (total == 0.0) {
            return {0, 0};
        }

        return {total*x_length, (px*p.first*p.second + f.first*f.second) / total * (-fwd*x_length) + x_fore};
    }

    std::pair<double,double> get_drag(int side, bool account_friction) {

        int const px = side == 1 ? 1 : -1;
        side--;
        auto const p = trapz_and_centre(panels_y[side], pressures[side]);
        auto const f = account_friction ? trapz_and_centre(panels_x[side], taus[side]) : std::make_pair(0.0, 0.0);
        double const total = px*p.first + f.first;
        if (total == 0.0) {
            return {0, 0};
        }

        return {total*x_length, (px*p.first*p.second + f.first*f.second) / total * x_length + top};
    }


private:

    bool iterate(int &m_Iterations, int s_IterLim) {

        if (!foil->viscal()) {
            foil->lvconv = false;
            std::cout << "CpCalc: local speed too large\nCompressibility corrections invalid" << std::endl;
            return false;
        }

        while (m_Iterations < s_IterLim && !foil->lvconv) {
            if (foil->ViscousIter()) {
                m_Iterations++;
            } else
                m_Iterations = s_IterLim;
        }

        if (!foil->ViscalEnd()) {
            foil->lvconv = false;  // point is unconverged
            foil->setBLInitialized(false);
            foil->lipan = false;
            return true;  // to exit loop
        }

        if (m_Iterations >= s_IterLim && !foil->lvconv) {
            foil->setBLInitialized(false);
            foil->lipan = false;
            foil->fcpmin();  // Is it of any use ?
            return true;
        }
        if (!foil->lvconv) {
            foil->fcpmin();  // Is it of any use ?
            return false;
        } else {
            // converged at last
            foil->fcpmin();  // Is it of any use ?
            return true;
        }
        return false;
    }

};

#endif // STREAMLINE_H
