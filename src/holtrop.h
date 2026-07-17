#ifndef HOLTROP_H
#define HOLTROP_H

#include "hull.h"

class Holtrop {
    std::shared_ptr<Hull> hull;
    std::shared_ptr<Environment> env;
public:
    ~Holtrop() {}
    Holtrop() {}
    Holtrop(std::shared_ptr<Hull> ship_hull) : hull(ship_hull) {
        env = hull->get_environment();
    }

    double get_viscous_pressure_resistance(double speed, double friction_known) {

        auto B = hull->get_beam_wl();
        auto L = std::max(0.8*hull->get_length_ref(), hull->get_length_wl());
        auto T = hull->get_draft();
        auto vol = hull->get_displacement();
        auto lcb = hull->get_fwd() * 100*(hull->get_cb().x - 0.5*(hull->get_ap() + hull->get_fp())) / hull->get_length_wl();
        auto Cp = std::max(0.4, hull->get_prismatic_coef());
        auto Lr = L * (1 - Cp + 0.06 * Cp * lcb / (4*Cp - 1));
        auto Cstern = 0.0;
        auto c14 = 1 + 0.011*Cstern;
        auto c22 = (Cp < 0.9) ? 0.0 : 0.4*(Cp - 0.9);
        auto form_factor = 0.93 + 0.487118 * c14 * std::pow(B/L, 1.06806) * std::pow(T/L, 0.46106) * std::pow(L/Lr, 0.121563) * std::pow(cub(L)/vol, 0.36486) * std::pow(1-Cp+c22, -0.604247);

        form_factor *= hull->get_form_factor_correction(hull->get_speed(0.1), speed);
        if (std::isnan(form_factor)) {
            form_factor = 1.0;
        }

        form_factor = std::max(form_factor, 1.0);

        auto R_exp = (form_factor-1) * friction_known;
        return R_exp;
    }

    double get_wave_resistance(double speed) {

        if (speed <= 0.01) {
            return 0.0;
        }

        auto vol = hull->get_displacement();
        auto B = hull->get_beam_wl();
        auto L = hull->get_length_wl();
        auto T = hull->get_draft();
        auto lcb = hull->get_fwd() * 100*(hull->get_cb().x - 0.5*(hull->get_ap() + hull->get_fp()) /*hull->get_x_beam_wl()*/) / hull->get_length_wl();
        auto Cp = std::max(0.4, hull->get_prismatic_coef());   // guard 4Cp-1 (div-by-zero at Cp=0.25)
        auto Lr = hull->get_length_wl() * (1 - Cp + 0.06*Cp*lcb/(4*Cp-1));
        auto Cwp = hull->get_waterplane_coef();
        auto CM = hull->get_max_section_coef();

        double C7;
        if (B/L < 0.11) {
            C7 = 0.229577 * std::pow(B/L, 1.0/3.0);
        } else if (B/L < 0.25) {
            C7 = B/L;
        } else {
            C7 = 0.5 - 0.0625*(L/B);
        }

        auto iE = 1 + 89 * std::exp(-std::pow(L/B, 0.80856) * std::pow(1-Cwp, 0.30484)
                    * std::pow(1-Cp-0.0225*lcb, 0.6367) * std::pow(Lr/B, 0.34574)
                    * std::pow(100*vol/cub(L), 0.16302));
        iE += std::abs(std::atan(hull->get_beam_wl() / hull->get_x_shoulder())) * 180/M_PI;
        iE /= 2.0;
        if (std::isnan(iE)) { iE = 89.0; }   // std::min does NOT filter NaN (pow of a negative base upstream)
        iE = std::min(iE, 89.0);   // keep 90-iE > 0 so pow(90-iE, -1.37565) stays finite

        auto C1 = 2223105*std::pow(C7, 3.78613)*std::pow(T/B, 1.07961)*std::pow(90.0-iE, (-1.37565));
        // std::pow(L/B-2, frac) is NaN when L/B<2 (beamy hulls that shouldn't be on the Holtrop path);
        // clamp the base to keep C17 finite.
        auto C17 = 6919.3*std::pow(CM, (-1.3346))*std::pow(vol/cub(L), 2.00977) * std::pow(std::max(L/B-2.0, 0.0), 1.40692);
        auto C3 = 0.0; //0.56*(Abulb**(1.5))/(B*T*(0.31*(Abulb**(0.5))+T_fore-Hbulb))
        auto C2 = std::exp(-1.89 * std::sqrt(C3));
        auto At = hull->get_transom_area();
        auto C5 = (At > 0.0) ? std::max(0.0, 1.0 - 0.8*At/(B*T*CM)) : 1.0;

        double lamda;
        if (L/B < 12.0) {
            lamda = (1.446*Cp)-(0.03*L/B);
        } else {
            lamda = (1.446*Cp)-0.36;
        }

        double C16;
        if (Cp < 0.8) {
            C16 = 8.07981*Cp-13.8673*sq(Cp)+6.984388*cub(Cp);
        } else {
            C16 = 1.73014-0.7067*Cp;
        }

        auto m1 = 0.0140407*L/T - 1.75254 * std::pow(vol, 1.0/3.0)/L - 4.79323*B/L - C16;
        auto m3 = -7.2035*std::pow(B/L, 0.326869)*std::pow(T/B, 0.605375);
        auto d = -0.9;

        auto C15 = 0.0;
        if (cub(L)/vol < 512.0) {
            C15 = -1.69385;
        } else if (cub(L)/vol < 1726.91) {
            C15 = -1.69385 + (L/std::pow(vol, 1.0/3.0) - 8.0) / 2.36;
        }

        auto Fn = hull->get_fn(speed);
//        auto m2 = C15 * sq(Cp) * std::exp(-0.1/sq(Fn)); // 1982
        auto m4 = C15 * 0.4 * std::exp(-0.034*std::pow(Fn, -3.29)); // 1984 // martina C16 ????

        double Rw = 0.0;
        if (Fn < 0.4) {
            Rw = (C1*C2*C5*vol * env->get_density() * env->get_gravity() * std::exp(m1*std::pow(Fn,d) + m4*std::cos(lamda/sq(Fn))));
        } else if (Fn < 0.55) {
            // Freeze m4 at each blend endpoint: the published 1984 method interpolates the COMPLETE
            // Fn=0.40 and Fn=0.55 wave-resistance formulas (m4 included), not m4 at the running Fn.
            double const m4_040 = C15 * 0.4 * std::exp(-0.034*std::pow(0.40, -3.29));
            double const m4_055 = C15 * 0.4 * std::exp(-0.034*std::pow(0.55, -3.29));
            auto Rw_040 = (C1*C2*C5*vol * env->get_density() * env->get_gravity() * std::exp(m1*std::pow(0.4,d) + m4_040*std::cos(lamda/sq(0.4))));
            auto Rw_055 = (C17*C2*C5*vol * env->get_density() * env->get_gravity() * std::exp(m3*std::pow(0.55,d) + m4_055*std::cos(lamda/sq(0.55))));
            Rw = Rw_040 + (10*Fn-4)*(Rw_055-Rw_040)/1.5;
        } else { // Fn > 0.55:
            Rw = (C17*C2*C5*vol * env->get_density() * env->get_gravity() * std::exp(m3*std::pow(Fn,d) + m4 * std::cos(lamda/sq(Fn))));
        }

        return Rw;
    }
};

#endif // HOLTROP_H
