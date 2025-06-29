#pragma once

#include "typedefs.h"
#include <algorithm>
// #include <stdexcept>

#include "AdaptiveIntegrators.h"
#include "FixedIntegrators.h"
#include <vector>

const int N_DEFAULT = 1000;
const f64 TOL_DEFAULT = 1e-9;
const f64 DT_DEFAULT = 1e-3;

// Fixed Step Size
template <
    typename F,
    typename State,
    template <typename, typename> class Policy>
struct FixedStepIntegrator {
    // Members
    F f_;
    f64 dt0_;

    // Constructors
    FixedStepIntegrator(F f, f64 dt0 = DT_DEFAULT)
        : f_(std::move(f)), dt0_(dt0) {};
    FixedStepIntegrator(F f, f64 t0, f64 tf, int N_steps = N_DEFAULT)
        : f_(std::move(f)) {
        f64 t_span = tf - t0;
        dt0_ = t_span / static_cast<f64>(N_steps);
    };

    // Integrate from t0 to tf
    void integrate(
        double t0,
        double tf,
        const State &x0,
        std::vector<double> &times,
        std::vector<State> &states
    ) const {
        times.clear();
        states.clear();

        // Initialize states and storage
        double t = t0;
        State x = x0;
        times.push_back(t);
        states.push_back(x);

        // Integration Loop
        while (t < tf) {
            double dt = std::min(dt0_, tf - t);
            auto [t_new, x_new] = Policy<F, State>::step(f_, t, x, dt);
            t = t_new;
            x = x_new;
            times.push_back(t_new);
            states.push_back(x_new);
        }
    }
};

// Adaptive Step Size
template <
    typename State,
    typename F,
    template <typename, typename> class Policy>
struct AdaptiveStepIntegrator {
    // Members
    F f_;
    f64 dt_;
    f64 tol_;

    // Constructors
    AdaptiveStepIntegrator(F f, f64 dt0 = DT_DEFAULT, f64 tol = TOL_DEFAULT)
        : f_(std::move(f)), dt_(dt0), tol_(tol) {};

    // Integrate from t0 to tf
    void integrate(
        double t0,
        double tf,
        const State &x0,
        std::vector<double> &times,
        std::vector<State> &states
    ) {
        times.clear();
        states.clear();

        // Initialize states and storage
        double t = t0;
        State x = x0;
        times.push_back(t);
        states.push_back(x);

        // Integration Loop
        while (t < tf) {
            double dt = std::min(dt_, tf - t);

            State x_trial;
            f64 dt_new, err;

            // Search for new dt
            while (true) {
                std::tie(dt_new, x_trial, err)
                    = Policy<F, State>::step(f_, t, x, dt, tol_);
                if (err <= tol_) {
                    // New dt found
                    break;
                }
                dt = dt_new;
                if (dt < std::numeric_limits<double>::epsilon()) {
                    // dt lower than precision
                    throw std::runtime_error(
                        "Adaptive integrator step size underflow"
                    );
                }
            }
            x = x_trial;
            t += dt;
            dt_ = dt_new;
            times.push_back(t);
            states.push_back(x);
        }
    }
};
