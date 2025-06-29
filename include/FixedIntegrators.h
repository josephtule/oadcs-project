#pragma once

#include "typedefs.h"
#include <utility>

// -----------------------------------------------------------------------------
// Fixed Step Integrators
// -----------------------------------------------------------------------------

// RK1/Euler
template <typename F, typename State> struct RK1Policy {
    static std::pair<f64, State>
    step(const F &f, double t, const State &x, double dt) {
        State k1 = f(t, x);
        State x_new = x + dt * k1;
        f64 t_new = t + dt;
        return {t_new, x_new};
    }
};

// RK2/Midpoint
template <typename F, typename State> struct RK2Policy {
    static std::pair<f64, State>
    step(const F &f, double t, const State &x, double dt) {
        State k1 = f(t, x);
        State k2 = f(t + dt / 2., x + dt * k1 / 2.);
        State x_new = x + dt * k2;
        f64 t_new = t + dt;
        return {t_new, x_new};
    }
};

// RK3
template <typename F, typename State> struct RK3Policy {
    static std::pair<f64, State>
    step(const F &f, double t, const State &x, double dt) {
        State k1 = f(t, x);
        State k2 = f(t + dt / 2., x + dt * k1 / 2.);
        State k3 = f(t + dt, x - dt * k1 + 2. * dt * k2);
        State x_new = x + dt * (k1 + 4. * k2 + k3) / 6.;
        f64 t_new = t + dt;
        return {t_new, x_new};
    }
};

// RK4
template <typename F, typename State> struct RK4Policy {
    static std::pair<f64, State>
    step(const F &f, double t, const State &x, double dt) {
        State k1 = f(t, x);
        State k2 = f(t + dt / 2., x + dt * k1 / 2);
        State k3 = f(t + dt / 2., x + dt * k2 / 2);
        State k4 = f(t + dt, x + dt * k3);
        State x_new = x + dt / 6. * (k1 + 2 * k2 + 2 * k3 + k4);
        f64 t_new = t + dt;
        return {t_new, x_new};
    }
};

// RK5
template <typename F, typename State> struct RK5Policy {
    static std::pair<f64, State>
    step(const F &f, double t, const State &x, double dt) {
        State k1 = f(t, x);
        State k2 = f(t + dt / 4., x + dt * k1 / 4.);
        State k3 = f(t + dt / 4., x + dt * (k1 / 8. + k2 / 8.));
        State k4 = f(t + dt / 2., x + dt * (-k2 / 2. + k3));
        State k5
            = f(t + 3. * dt / 4., x + dt * (3. * k1 / 16. + 9. * k4 / 16.));
        State k6
            = f(t + dt,
                x
                    + dt
                          * (-3. * k1 / 7. + 2. * k2 / 7. + 12. * k3 / 7.
                             - 12. * k4 / 7. + 8. * k5 / 7.));
        State x_new = x
                      + dt
                            * (7. * k1 / 90. + 32. * k3 / 90. + 12. * k4 / 90.
                               + 32. * k5 / 90. + 7. * k6 / 90.);
        f64 t_new = t + dt;
        return {t_new, x_new};
    }
};

// RK6
template <typename F, typename State> struct RK6Policy {
    static std::pair<f64, State>
    step(const F &f, double t, const State &x, double dt) {
        State k1 = f(t, x);
        State k2 = f(t + dt / 5., x + dt * k1 / 5.);
        State k3 = f(t + dt / 5., x + dt * (3. * k1 / 40. + 9. * k2 / 40.));
        State k4
            = f(t + dt / 2.,
                x + dt * (3. * k1 / 10. - 9. * k2 / 10. + 6. * k3 / 10.));
        State k5
            = f(t + 3. * dt / 4.,
                x
                    + dt
                          * (-11. * k1 / 54. + 5. * k2 / 2. - 70. * k3 / 27.
                             + 35. * k4 / 27.));
        State k6
            = f(t + dt,
                x
                    + dt
                          * (1631. * k1 / 55296. + 175. * k2 / 512.
                             + 575. * k3 / 13824. + 44275. * k4 / 110592.
                             + 253. * k5 / 4096.));
        State k7
            = f(t + dt,
                x
                    + dt
                          * (37. * k1 / 378. + 0. * k2 + 250. * k3 / 621.
                             + 125. * k4 / 594. + 0. * k5 + 512. * k6 / 1771.));
        State x_new
            = x
              + dt
                    * (37. * k1 / 378. + 0. * k2 + 250. * k3 / 621.
                       + 125. * k4 / 594. + 0. * k5 + 512. * k6 / 1771.);
        f64 t_new = t + dt;
        return {t_new, x_new};
    }
};

// Heun
template <typename F, typename State> struct HeunPolicy {
    static std::pair<f64, State>
    step(const F &f, double t, const State &x, double dt) {
        State k1 = f(t, x);
        State k2 = f(t + dt, x + dt * k1);
        State x_new = x + dt * (k1 + k2) / 2.;
        f64 t_new = t + dt;
        return {t_new, x_new};
    }
};

// Ralston
template <typename F, typename State> struct RalstonPolicy {
    static std::pair<f64, State>
    step(const F &f, double t, const State &x, double dt) {
        State k1 = f(t, x);
        State k2 = f(t + 3. * dt / 4., x + 3. * dt * k1 / 4.);
        State x_new = x + dt * (k1 / 3. + 2. * k2 / 3.);
        f64 t_new = t + dt;
        return {t_new, x_new};
    }
};