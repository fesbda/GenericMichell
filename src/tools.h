#ifndef TOOLS_H
#define TOOLS_H

#define _USE_MATH_DEFINES

#include <iostream>
#include <iterator>
#include <vector>
#include <map>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <memory>

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/geometric.hpp>

template<typename T> T sq(T x) { return x*x; }
template<typename T> T cub(T x) { return x*x*x; }

template<int EX> double pown(double x) { return std::pow(x, EX); }
// Explicit full specializations are ordinary functions (not templates): they need `inline`
// to avoid an ODR clash if tools.h is included in more than one translation unit.
template<> inline double pown<0>(double /*x*/) { return 1; }
template<> inline double pown<1>(double x) { return x; }
template<> inline double pown<2>(double x) { return sq(x); }
template<> inline double pown<3>(double x) { return cub(x); }
template<> inline double pown<4>(double x) { return sq(sq(x)); }
template<> inline double pown<5>(double x) { return sq(x)*cub(x); }
template<> inline double pown<6>(double x) { return cub(x)*cub(x); }

inline int sign(double v)  { return (0.0 < v) - (v < 0.0); }
inline int sign(double v, double tolerance)  { return (tolerance < v) - (v < -tolerance); }

template<typename T> T clamp(T x, T a, T b) { return std::max(std::min(x,b), a); }

template<typename T> T lerp(T a, T b, T alpha) { return a + (b - a) * alpha; }

inline double radians(double degrees) { return degrees * (M_PI / 180.0); }

inline double degrees(double radians) { return radians * (180.0 / M_PI); }

inline double round(double num, int decimal_places) {
    return std::round(num * std::pow(10, decimal_places)) / std::pow(10, decimal_places);
}

inline double exp_approx(double x) {
    x = 1.0 + x / 1024;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x;
    return x;
}

template<typename T>
std::pair<T,T> trapz_and_centre(std::vector<T> xrange, std::vector<T> yrange) {
    if (xrange.size() < 2) {
        return {0, 0};
    }
    T moments = 0;
    T sum = 0;
    for (std::size_t i = 1; i < xrange.size(); i++) {
        T const mass = 0.5 * (yrange[i] + yrange[i-1]) * (xrange[i] - xrange[i-1]);
        sum += mass;
        moments += mass * (xrange[i] + xrange[i-1]) * 0.5;
    }
    if (std::abs(sum) < std::numeric_limits<T>::min()) {
        return {sum, 0};
    }
    return {sum, moments / sum};
}

template<typename T>
T trapz(std::vector<T> xrange, std::vector<T> yrange) {
    T sum = 0;
    for (std::size_t i = 1; i < xrange.size(); i++) {
        sum += (yrange[i] + yrange[i-1]) * (xrange[i] - xrange[i-1]);
    }
    return sum / T(2);
}

template<typename T>
T trapz(T dx, std::vector<T> yrange) {
    T sum = 0;
    for (std::size_t i = 1; i < yrange.size(); i++) {
        sum += (yrange[i] + yrange[i-1]) * (dx);
    }
    return sum / T(2);
}

template <typename T>
std::vector<T> linspace(T a, T b, std::size_t N) {
    T const h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);

    for (std::size_t i = 0; i < N-1; i++) {
        xs[i] = a + h * i;
    }
    xs[N-1] = b; // last one separately to ignore machine errors for h*(n-1)

    return xs;
}

template<class Func>
double bisection(Func func, double target, double a, double b, double min_step, double abs_tolerance, bool* found = nullptr) {
    abs_tolerance = std::abs(abs_tolerance);
    auto c = a;
    while (std::abs(b - a) >= min_step) {
        c = (a + b) * 0.5;
        auto f_c = func(c) - target;
        if (std::abs(f_c) < abs_tolerance){
            if (found) *found = true;
            return c;
        } else if (f_c * (func(a) - target) < 0) {
            b = c;
        } else {
            a = c;
        }
    }
    if (found) *found = false;
    return c;
}

// Interpolate-Truncate-Project bracketed root finder (Oliveira & Takahashi 2020). The bracket [a,b] must
// satisfy g(a) < 0 < g(b) with g increasing (pass ga=g(a), gb=g(b)). g_at(x) evaluates g at x; budget_ok()
// returns false once the caller's evaluation budget is exhausted (an early-stop guard, like a max-iteration
// cap). Returns the midpoint of the final bracket. Bracketed => never worse than bisection (unconditionally
// convergent) and superlinear via false-position interpolation with minmax truncation + projection. eps is
// the half-width tolerance on the bracket; k1, k2, n0 are the standard ITP tuning parameters.
template<class GFunc, class BudgetFunc>
double itp_root(double a, double b, double ga, double gb,
                double eps, double k1, double k2, int n0,
                GFunc g_at, BudgetFunc budget_ok) {
    int const n_max = (int)std::ceil(std::log2(std::max(1.0, (b - a) / (2.0 * eps)))) + n0;
    int j = 0;
    while ((b - a) > 2.0 * eps && budget_ok() && (gb - ga) > 0.0) {
        double const xh = 0.5 * (a + b);
        double const xf = (a * gb - b * ga) / (gb - ga);     // regula-falsi (false position) on g
        double const sigma = (xh - xf) >= 0.0 ? 1.0 : -1.0;
        double const delta = k1 * std::pow(b - a, k2);
        double const xt = (delta <= std::abs(xh - xf)) ? xf + sigma * delta : xh;
        double const r = eps * std::pow(2.0, double(n_max - j)) - 0.5 * (b - a);
        double const x = (std::abs(xt - xh) <= r) ? xt : xh - sigma * r;
        double const gx = g_at(x);
        if (gx > 0.0) { b = x; gb = gx; }
        else if (gx < 0.0) { a = x; ga = gx; }
        else { a = x; b = x; break; }
        ++j;
    }
    return 0.5 * (a + b);
}

template <typename V>
bool intersect_edge_plane_x(V a, V b, double x, V& out) {

    if ((a.x < x && b.x < x) || (a.x > x && b.x > x)) {
        return false;
    }

    if (std::abs(b.x - a.x) < 1e-12) {       // edge parallel to the plane: avoid 0/0
        out = (a + b) * 0.5;
    } else {
        auto ratio = (x - a.x) / (b.x - a.x);
        out = a + (b - a) * ratio;
    }
    out.x = x;
    return true;
}

template <typename V>
bool intersect_edge_plane_y(V a, V b, double y, V& out) {

    if ((a.y < y && b.y < y) || (a.y > y && b.y > y)) {
        return false;
    }

    if (std::abs(b.y - a.y) < 1e-12) {       // edge parallel to the plane: avoid 0/0
        out = (a + b) * 0.5;
    } else {
        auto ratio = (y - a.y) / (b.y - a.y);
        out = a + (b - a) * ratio;
    }
    out.y = y;
    return true;
}

template <typename V>
bool intersect_edge_plane_z(V a, V b, double z, V& out) {

    if ((a.z < z && b.z < z) || (a.z > z && b.z > z)) {
        return false;
    }

    if (std::abs(b.z - a.z) < 1e-6) {
        out = (a + b) * 0.5;
    }
    else {
        auto ratio = (z - a.z) / (b.z - a.z);
        out = a + (b - a) * ratio;
    }
    out.z = z;
    return true;
}

template <typename V>
bool intersect_edges_2d(V p0, V p1, V p2, V p3, V& intersection) {

    auto const a = p1 - p0;
    auto const b = p0 - p2;

    auto const s = a.x*b.y - a.y*b.x;
    if (s < 0) {
        return false;
    }

    auto const c = p3 - p2;
    auto t = c.x*b.y - c.y*b.x;
    if (t < 0) {
        return false;
    }

    auto const denom = a.x*c.y - c.x*a.y;
    if (s > denom || t > denom) {
        return false;
    }

    t /= denom;
    intersection = p0 + a*t;
    return true;
}

template <typename V, typename T>
bool intersect_ray_triangle(const V& from, const V& dir, const V vertices[3], T& out_distance)
{
    auto const e1 = vertices[1] - vertices[0];
    auto const e2 = vertices[2] - vertices[0];
    auto const p = glm::cross(dir, e2);
    auto const det = glm::dot(e1, p);
    auto const sd = sign(det, T(1e-28));

    if(!sd) // ray is parallel
        return false;

    auto const t = from - vertices[0];
    auto const adet = det * sd;              // |det| (sd = sign(det)); u,v are sign-corrected
    auto const u = glm::dot(t, p) * sd;
    if (u < T(0) || u > adet)                // bound vs |det|, not raw det (was wrong for det<0)
        return false;

    auto const q = glm::cross(t, e1);
    auto const v = glm::dot(dir, q) * sd;
    if (v < T(0) || (u + v) > adet)
        return false;

    out_distance = glm::dot(e2, q) / det;
    return (out_distance > T(0));
}

template<typename V>
auto edge_sign(V p1, V p2, V p3) -> decltype(V::x) {
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

template<typename V>
bool is_point_in_2d_triangle (V pt, V v1, V v2, V v3) {
    bool b1 = edge_sign(pt, v1, v2) < 0.0;
    bool b2 = edge_sign(pt, v2, v3) < 0.0;
    bool b3 = edge_sign(pt, v3, v1) < 0.0;
    return ((b1 == b2) && (b2 == b3));
}

inline double rely(double x, double interval_start, double interval_end, double too_far=1.0) {
    if (interval_end - interval_start < 1e-16)
        return 0.0; // interval sucks
    if (x < interval_start)
        return std::exp(-5*(interval_start - x)/(too_far*(interval_end - interval_start)));
    if (x > interval_end)
        return std::exp(-5*(x - interval_end)/(too_far*(interval_end - interval_start)));
    return 1.0; // it's inside the range
}

// integrate[F(X) * COS(T*X) * dX] from a to b, with odd number of values
inline double filon_cos(double* ftab, int n, double a, double b, double t) {

    double alpha;
    double beta;
    double gamma;
    int i;

    if ( a == b ) {
        return 0;
    }

    if ( n <= 1 || ( n % 2 ) != 1 ) {
        return 0;
    }

    double *x = new double[n];
    for ( i = 0; i < n; i++ ) {
        x[i] = (double(n-i-1)*a + double(i)*b) / double(n-1);
    }

    auto h = (b-a) / double(n-1);
    auto theta = t * h;
    auto sint = std::sin(theta);
    auto cost = std::cos(theta);

    if (6.0 * std::abs(theta) <= 1.0) {
        alpha = 2.0 * std::pow ( theta, 3 ) /   45.0
              - 2.0 * std::pow ( theta, 5 ) /  315.0
              + 2.0 * std::pow ( theta, 7 ) / 4725.0;
        beta  = 2.0                         /     3.0
              + 2.0 * std::pow ( theta, 2 ) /    15.0
              - 4.0 * std::pow ( theta, 4 ) /   105.0
              + 2.0 * std::pow ( theta, 6 ) /   567.0
              - 4.0 * std::pow ( theta, 8 ) / 22275.0;
        gamma = 4.0                         /      3.0
              - 2.0 * std::pow ( theta, 2 ) /     15.0
              +       std::pow ( theta, 4 ) /    210.0
              -       std::pow ( theta, 6 ) /  11340.0;
    } else {
        alpha = (std::pow(theta, 2) + theta * sint * cost - 2.0 * sint * sint) / std::pow(theta, 3);
        beta = ( 2.0 * theta + 2.0 * theta * cost * cost - 4.0 * sint * cost ) / std::pow(theta, 3);
        gamma = 4.0 * (sint - theta*cost) / std::pow(theta, 3);
    }

    auto c2n = + 0.5 * ftab[0] * std::cos(t * x[0]);
    for ( i = 2; i < n - 1; i = i + 2 ) {
        c2n = c2n + ftab[i] * std::cos(t * x[i]);
    }
    c2n = c2n + 0.5 * ftab[n-1] * std::cos ( t * x[n-1] );

    auto c2nm1 = 0.0;
    for ( i = 1; i <= n - 2; i = i + 2 ) {
        c2nm1 = c2nm1 + ftab[i] * std::cos ( t * x[i] );
    }

    auto value = h * (
              alpha * ( ftab[n-1] * std::sin ( t * x[n-1] )
            - ftab[0]   * std::sin ( t * x[0] ) )
            + beta * c2n
            + gamma * c2nm1 );

    delete [] x;
    return value;
}

inline double filon_cos(std::vector<double> ftab, double a, double b, double t) {
    return filon_cos(ftab.data(), int(ftab.size()), a, b, t);
}

// integrate[F(X) * SIN(T*X) * dX] from a to b, with odd number of values
inline double filon_sin(double* ftab, int n, double a, double b, double t) {

    double alpha;
    double beta;
    double gamma;
    int i;

    if ( a == b ) {
        return 0;
    }

    if ( n <= 1 || ( n % 2 ) != 1 ) {
        return 0;
    }

    double *x = new double[n];
    for ( i = 0; i < n; i++ ) {
        x[i] = (double(n-i-1)*a + double(i)*b) / double(n-1);
    }

    auto h = (b-a) / double(n-1);
    auto theta = t * h;
    auto sint = std::sin(theta);
    auto cost = std::cos(theta);

    if (6.0 * std::abs ( theta ) <= 1.0) {
        alpha = 2.0 * std::pow ( theta, 3 ) /   45.0
                - 2.0 * std::pow ( theta, 5 ) /  315.0
                + 2.0 * std::pow ( theta, 7 ) / 4725.0;

        beta =  2.0                    /     3.0
                + 2.0 * std::pow ( theta, 2 ) /    15.0
                - 4.0 * std::pow ( theta, 4 ) /   105.0
                + 2.0 * std::pow ( theta, 6 ) /   567.0
                - 4.0 * std::pow ( theta, 8 ) / 22275.0;

        gamma = 4.0                    /      3.0
                - 2.0 * std::pow ( theta, 2 ) /     15.0
                +       std::pow ( theta, 4 ) /    210.0
                -       std::pow ( theta, 6 ) /  11340.0;
    } else {
        alpha = ( std::pow ( theta, 2 ) + theta * sint * cost
                  - 2.0 * sint * sint ) / std::pow ( theta, 3 );

        beta = ( 2.0 * theta + 2.0 * theta * cost * cost
                 - 4.0 * sint * cost ) / std::pow ( theta, 3 );

        gamma = 4.0 * ( sint - theta * cost ) / std::pow ( theta, 3 );
    }

    auto s2n = + 0.5 * ftab[0] * std::sin ( t * x[0] );
    for ( i = 2; i < n - 1; i = i + 2 ) {
        s2n = s2n + ftab[i] * std::sin ( t * x[i] );
    }
    s2n = s2n + 0.5 * ftab[n-1] * std::sin ( t * x[n-1] );

    auto s2nm1 = 0.0;
    for ( i = 1; i <= n - 2; i = i + 2 ) {
        s2nm1 = s2nm1 + ftab[i] * std::sin ( t * x[i] );
    }

    auto value = h * (
                alpha * ( ftab[0]   * std::cos ( t * x[0] )
            - ftab[n-1] * std::cos ( t * x[n-1] ) )
            + beta * s2n
            + gamma * s2nm1 );

    delete [] x;

    return value;
}

inline double filon_sin(std::vector<double> ftab, double a, double b, double t) {
    return filon_sin(ftab.data(), int(ftab.size()), a, b, t);
}

// Exact integral over one panel [xa,xb] of the LINEAR interpolant between (xa,fa) and (xb,fb)
// times cos(t*x) (filon_cos_seg) or sin(t*x) (filon_sin_seg). Building block for the causal
// near-field accumulation (SOLVER_FIX_PLAN.md Batch 3): panels are summed to prefix quadratures
// A=∫F cos(t s)ds, B=∫F sin(t s)ds. Small-|theta| Taylor guard mirrors filon_cos (avoids the 1/θ,
// 1/θ² cancellation). In the u=(x-xa)/h ∈[0,1] frame with φ=t*xa, θ=t*h:
//   ∫cos = h[fa·∫₀¹cos(φ+θu)du + (fb-fa)·∫₀¹u cos(φ+θu)du]   (and sin analogously).
inline double filon_cos_seg(double fa, double fb, double xa, double xb, double t) {
    double const h = xb - xa;
    if (h == 0.0) { return 0.0; }
    double const phi = t * xa, theta = t * h;
    double const c = std::cos(phi), s = std::sin(phi);
    double C0, C1;   // ∫₀¹cos(φ+θu)du , ∫₀¹u cos(φ+θu)du
    if (6.0 * std::abs(theta) <= 1.0) {
        double const t2 = theta*theta;
        double const E  = 1.0 - t2/6.0 + t2*t2/120.0 - t2*t2*t2/5040.0;          // Σ(-1)ⁿθ²ⁿ/((2n)!(2n+1))
        double const O  = theta*(0.5 - t2/24.0 + t2*t2/720.0 - t2*t2*t2/40320.0);// Σ(-1)ⁿθ²ⁿ⁺¹/((2n+1)!(2n+2))
        double const E1 = 0.5 - t2/8.0 + t2*t2/144.0 - t2*t2*t2/5760.0;          // Σ(-1)ⁿθ²ⁿ/((2n)!(2n+2))
        double const O1 = theta*(1.0/3.0 - t2/30.0 + t2*t2/840.0 - t2*t2*t2/45360.0);// Σ(-1)ⁿθ²ⁿ⁺¹/((2n+1)!(2n+3))
        C0 = c * E  - s * O;
        C1 = c * E1 - s * O1;
    } else {
        double const cb = std::cos(phi + theta), sb = std::sin(phi + theta);
        C0 = (sb - s) / theta;
        C1 = sb / theta + (cb - c) / (theta*theta);
    }
    return h * (fa * C0 + (fb - fa) * C1);
}

inline double filon_sin_seg(double fa, double fb, double xa, double xb, double t) {
    double const h = xb - xa;
    if (h == 0.0) { return 0.0; }
    double const phi = t * xa, theta = t * h;
    double const c = std::cos(phi), s = std::sin(phi);
    double S0, S1;   // ∫₀¹sin(φ+θu)du , ∫₀¹u sin(φ+θu)du
    if (6.0 * std::abs(theta) <= 1.0) {
        double const t2 = theta*theta;
        double const E  = 1.0 - t2/6.0 + t2*t2/120.0 - t2*t2*t2/5040.0;
        double const O  = theta*(0.5 - t2/24.0 + t2*t2/720.0 - t2*t2*t2/40320.0);
        double const E1 = 0.5 - t2/8.0 + t2*t2/144.0 - t2*t2*t2/5760.0;
        double const O1 = theta*(1.0/3.0 - t2/30.0 + t2*t2/840.0 - t2*t2*t2/45360.0);
        S0 = s * E  + c * O;
        S1 = s * E1 + c * O1;
    } else {
        double const cb = std::cos(phi + theta), sb = std::sin(phi + theta);
        S0 = (c - cb) / theta;
        S1 = -cb / theta + (sb - s) / (theta*theta);
    }
    return h * (fa * S0 + (fb - fa) * S1);
}

#endif // TOOLS_H
