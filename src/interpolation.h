#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>

struct LinearInterpolator {

    std::vector<double> xi;
    std::vector<double> yi;
    bool uniform;
    double dx_inv;

    LinearInterpolator() {}
    LinearInterpolator(std::vector<double> xs, std::vector<double> fs) {
        xi = xs;
        yi = fs;
        uniform = true;
        dx_inv = 0.0;
        if (xs.size() < 2) { return; }   // operator()/get_dy already short-circuit <2 points; avoid xs[1] OOB
        auto const tmp = xs[1] - xs[0];
        for (std::size_t i = 2; i < xs.size(); i++) {
            if (std::abs(xs[i] - xs[i-1] - tmp) > 1e-3*tmp) {
                uniform = false;
                break;
            }
        }
        dx_inv = 1.0 / tmp;
    }

    ~LinearInterpolator() {}

    double operator()(double x) const {

        if (yi.empty()) {
            return 0.0;
        } else if (yi.size() == 1) {
            return yi[0];
        }

        int i = uniform ? (int((x-xi[0])*dx_inv)+1) : int(std::upper_bound(xi.begin(), xi.end(), x) - xi.begin());
        if (i <= 0) { i = 1; }
        if (static_cast<std::size_t>(i) >= xi.size()) { i = int(xi.size()) - 1; }
        auto ratio = uniform ? (x-xi[i-1])*dx_inv : (x-xi[i-1])/(xi[i]-xi[i-1]);
        auto y = (1.0 - ratio) * yi[i-1] + ratio * yi[i];
        return y;
    }

    double get_dy(double x) const {

        if (yi.size() < 2) {
            return 0.0;
        }

        int i = uniform ? (int((x-xi[0])*dx_inv)+1) : int(std::upper_bound(xi.begin(), xi.end(), x) - xi.begin());
        if (i <= 0) { i = 1; }
        if (static_cast<std::size_t>(i) >= xi.size()) { i = int(xi.size()) - 1; }
        auto dy = yi[i] - yi[i-1];

        if (uniform) {
            return dy * dx_inv;
        }

        auto dx = xi[i] - xi[i-1];
        if (std::abs(dx) < 1e-18) {
            return 0;
        }
        return dy/dx;
    }

};

struct PchipInterpolator {

    std::vector<double> xi;
    std::vector<double> yi;

    PchipInterpolator() {}

    PchipInterpolator(std::vector<double> xs, std::vector<double> fs) {
        xi = xs;
        yi = fs;
    }

    ~PchipInterpolator() {}

    double operator()(double x) const {

        if (yi.empty()) {
            return 0.0;
        } else if (yi.size() == 1) {
            return yi[0];
        }

        auto N = xi.size();

        if ( x <= xi[0] || N<=2 ) {
            double dx = (x-xi[0])/(xi[1]-xi[0]);
            return (1.0-dx)*yi[0] + dx*yi[1];
        } else if ( x >= xi[N-1] ) {
            double dx = (x-xi[N-2])/(xi[N-1]-xi[N-2]);
            return (1.0-dx)*yi[N-2] + dx*yi[N-1];
        }
        auto i = size_t(std::upper_bound(xi.begin(), xi.end(), x) - xi.begin());
        double f1 = yi[i-1];
        double f2 = yi[i];
        double dx = (x-xi[i-1])/(xi[i]-xi[i-1]);
        // Compute the gradient in normalized coordinates [0,1]
        double g1=0, g2=0;
        if ( i<=1 ) {
            g1 = f2-f1;
        } else if ( ( f1<f2 && f1>yi[i-2] ) || ( f1>f2 && f1<yi[i-2] ) ) {
            // Compute the gradient by using a 3-point finite difference to f'(x)
            // Note: the real gradient is g1/(xi[i]-xi[i-1])
            double f0 = yi[i-2];
            double dx1 = xi[i-1]-xi[i-2];
            double dx2 = xi[i]-xi[i-1];
            double a1 = (dx2-dx1)/dx1;
            double a2 = dx1/(dx1+dx2);
            g1 = a1*(f1-f0)+a2*(f2-f0);
            double g_max = 2*dx2*std::min(fabs(f1-f0)/dx1,fabs(f2-f1)/dx2);
            g1 = ((g1>=0)?1:-1)*std::min(fabs(g1),g_max);
        }
        if ( i>=N-1 ) {
            g2 = f2-f1;
        } else if ( ( f2<f1 && f2>yi[i+1] ) || ( f2>f1 && f2<yi[i+1] ) ) {
            // Compute the gradient by using a 3-point finite difference to f'(x)
            // Note: the real gradient is g2/(xi[i]-xi[i-1])
            double f0 = yi[i+1];
            double dx1 = xi[i]-xi[i-1];
            double dx2 = xi[i+1]-xi[i];
            double a1 = -dx2/(dx1+dx2);
            double a2 = (dx2-dx1)/dx2;
            g2 = a1*(f1-f0)+a2*(f2-f0);
            double g_max = 2*dx1*std::min(fabs(f2-f1)/dx1,fabs(f0-f2)/dx2);
            g2 = ((g2>=0)?1:-1)*std::min(fabs(g2),g_max);
        }
        // Perform the interpolation
        double dx2 = dx*dx;
        double f = f1 + dx2*(2*dx-3)*(f1-f2) + dx*g1 - dx2*(g1+(1-dx)*(g1+g2));
        return f;
    }

};


class BiLinearInterpolator {
public:
    std::vector<double> xi;
    std::vector<double> yi;
    std::vector<double> fi;
    bool x_fwd = true;
    bool y_fwd = true;

    BiLinearInterpolator() {}

    BiLinearInterpolator(std::vector<double> xs, std::vector<double> ys, std::vector<double> fs) {
        xi = xs;
        yi = ys;
        fi = fs;
        x_fwd = xs.size() < 2 || (xs[1] - xs[0]) > 0.0;   // guard xs[1]/ys[1] OOB for degenerate axes
        y_fwd = ys.size() < 2 || (ys[1] - ys[0]) > 0.0;
    }

    ~BiLinearInterpolator() {}

    double get_max() const {
        return *std::max_element(fi.begin(), fi.end());
    }

    double get_min() const {
        return *std::min_element(fi.begin(), fi.end());
    }

    double operator()(double x, double y) const {

        if (fi.empty()) {
            return 0.0;
        } else if (fi.size() == 1) {
            return fi[0];
        }

        std::size_t N = xi.size();
        std::size_t M = yi.size();
        std::size_t i, j;

        if (x_fwd) {
            for (i = 0; (i < N) && (x > xi[i]); i++);
        } else {
            for (i = N-1; (i > 0) && (x > xi[i]); i--);
            i++;
        }

        if (y_fwd) {
            for (j = 0; (j < M) && (y > yi[j]); j++);
        } else {
            for (j = M-1; (j > 0) && (y > yi[j]); j--);
            j++;
        }

        if ( i <= 0 ) { i = 1; }
        if ( j <= 0 ) { j = 1; }
        if ( i >= N ) { i = N-1; }
        if ( j >= M ) { j = M-1; }
        // Clamp to [0,1] on BOTH sides: the bare std::min clamped the high edge but let a query
        // below the low edge extrapolate (dx<0). Clamp low too so out-of-range queries saturate.
        auto dx = std::max(std::min(( x - xi[i-1] ) / ( xi[i] - xi[i-1] ), 1.0), 0.0);
        auto dy = std::max(std::min(( y - yi[j-1] ) / ( yi[j] - yi[j-1] ), 1.0), 0.0);
        auto f1 = fi[i-1+(j-1)*N];
        auto f2 = fi[i+(j-1)*N];
        auto f3 = fi[i-1+j*N];
        auto f4 = fi[i+j*N];
        auto dx2 = 1.0-dx;
        auto dy2 = 1.0-dy;
        auto f = (dx*f2 + dx2*f1)*dy2 + (dx*f4 + dx2*f3)*dy;
        return f;
    }

};


#endif // INTERPOLATION_H
