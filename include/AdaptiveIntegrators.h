#pragma once

#include "typedefs.h"

// Dormand-Prince 4(5) (ODE45)
template <typename State, typename F> struct DormandPrince45Policy {
    static std::tuple<f64, State, f64>
    step(const F &f, f64 t, const State &x, f64 dt, f64 tol) {
        // Coefficients
        const f64 c2 = 1.0 / 5.0;
        const f64 c3 = 3.0 / 10.0;
        const f64 c4 = 4.0 / 5.0;
        const f64 c5 = 8.0 / 9.0;
        // ks
        State k1 = f(t, x);
        State k2 = f(t + c2 * dt, x + dt * (1.0 / 5.0 * k1));
        State k3
            = f(t + c3 * dt, x + dt * ((3.0 / 40.0) * k1 + (9.0 / 40.0) * k2));
        State k4
            = f(t + c4 * dt,
                x
                    + dt
                          * ((44.0 / 45.0) * k1 - (56.0 / 15.0) * k2
                             + (32.0 / 9.0) * k3));
        State k5
            = f(t + c5 * dt,
                x
                    + dt
                          * ((19372.0 / 6561.0) * k1 - (25360.0 / 2187.0) * k2
                             + (64448.0 / 6561.0) * k3 - (212.0 / 729.0) * k4));
        State k6
            = f(t + dt,
                x
                    + dt
                          * ((9017.0 / 3168.0) * k1 - (355.0 / 33.0) * k2
                             + (46732.0 / 5247.0) * k3 + (49.0 / 176.0) * k4
                             - (5103.0 / 18656.0) * k5));
        // 5th order
        State y5 = x
                   + dt
                         * ((35.0 / 384.0) * k1 + (500.0 / 1113.0) * k3
                            + (125.0 / 192.0) * k4 - (2187.0 / 6784.0) * k5
                            + (11.0 / 84.0) * k6);
        // 4th order
        State y4
            = x
              + dt
                    * ((5179.0 / 57600.0) * k1 + (7571.0 / 16695.0) * k3
                       + (393.0 / 640.0) * k4 - (92097.0 / 339200.0) * k5
                       + (187.0 / 2100.0) * k6 + (1.0 / 40.0) * f(t + dt, y5));
        f64 err = (y5 - y4).norm();
        f64 safety = 0.9;
        f64 dt_new
            = dt * std::clamp(safety * std::pow(tol / err, 0.2), 0.2, 5.0);
        return {dt_new, y5, err};
    }
};