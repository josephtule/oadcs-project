
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "CelestialBody.h"
#include "typedefs.h"

template <typename State, typename ForcePolicy> struct EOM {
    ForcePolicy forces;

    // Constructor
    explicit EOM(const ForcePolicy &fp) : forces(fp) {};

    State operator()(f64 t, State &x) const {
        State dxdt;
        dxdt << x.template segment<3>(0), forces.acceleration(x);

        return dxdt;
    }
};

template <typename State, typename... Policies> struct ForcePolicy3D {
    std::tuple<Policies...> policies;

    explicit ForcePolicy3D(const Policies &...ps) : policies(ps...) {};

    vec3 acceleration(const State &x) const {
        vec3 a_total = vec3::Zero();

        std::apply(
            [&](auto const &...p) { ((a_total += p.acceleration(x)), ...); },
            policies
        );

        return a_total;
    }
};

template <typename State> struct NewtonianGravityPolicy {
    f64 mu;

    explicit NewtonianGravityPolicy(f64 mu) : mu(mu) {};

    vec3 acceleration(const State &x) {
        vec3 r = x.template segment<3>(0);
        f64 r_mag = r.norm();

        vec3 a = -mu / (r_mag * r_mag * r_mag) * r;
        return a;
    }
};

template <typename State> struct ZonalGravityPolicy {
    f64 mu;
    std::vector<f64> J;
    f64 R_cb;
    int max_degree;

    explicit ZonalGravityPolicy(
        f64 mu,
        std::vector<f64> J,
        f64 R_cb,
        int max_degree
    )
        : mu(mu), J(J), R_cb(R_cb), max_degree(max_degree) {}

    vec3 acceleration(const State &x) {
        // Newontian Gravity
        NewtonianGravityPolicy<vec6> newton(mu);
        vec3 a = newton.acceleration(x);

        // Add Zonal Gravity Perturbations (only J2 for now)
        vec3 r = x.template segment<3>(0);
        f64 J2 = J[1];
        f64 r_mag = r.norm();
        f64 r_mag2 = r.squaredNorm();
        f64 Ror = R_cb / r_mag;
        f64 Ror2 = Ror * Ror;
        f64 muor2 = mu / r_mag2;
        f64 r0or = r(0) / r_mag;
        f64 r1or = r(1) / r_mag;
        f64 r2or = r(2) / r_mag;
        f64 r2or2 = r2or * r2or;

        f64 coef = -3. / 2. * J2 * muor2 * Ror2;
        vec3 vec_comp = vec3(
            (1. - 5. * r2or2) * r0or, //
            (1. - .5 * r2or2) * r1or, //
            (3. - 5. * r2or2) * r2or  //
        );

        a += coef * vec_comp;
        return a;
    }
};

template <typename State> struct SphericalHarmonicGravityPolicy {
    f64 mu;
    f64 R_cb;
    EGMCoefficients &egm;
    int max_degree, max_order;

    explicit SphericalHarmonicGravityPolicy(
        
    )
        : mu(mu), R_cb(R_cb), max_degree(egm.max_degree),
          max_order(egm.max_order), egm(egm) {}
};
