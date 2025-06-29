#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

#include "FixedIntegrators.h"
#include "integrator.h"
#include "translational.h"
#include "typedefs.h"

int main() {
    // Gravitational parameter for Earth (km^3/s^2)
    const f64 mu = 398600.4418;

    // Define the two-body ODE: state x = [r; v]
    // auto f = [mu](f64 t, const vec6 &x) {
    //     vec3 r = x.segment<3>(0);
    //     vec3 v = x.segment<3>(3);
    //     f64 r_norm = r.norm();

    //     vec3 a = -mu / (r_norm * r_norm * r_norm) * r;
    //     vec6 dx;
    //     dx.segment<3>(0) = v;
    //     dx.segment<3>(3) = a;
    //     return dx;
    // };
    auto f = [mu](f64 t, const vec6 &x) {
        vec6 dxdt = gravity_newton(t, x, mu);
        return dxdt;
    };

    // Initial circular orbit at radius = 7000 km
    vec3 r0, v0;
    f64 r0_mag = 7000.;
    r0 << r0_mag, 0, 0;
    v0 << 0, std::sqrt(mu / r0_mag), 0;
    vec6 x0;
    x0 << r0, v0;

    // Inital state print
    std::cout << "Initial State" << std::endl;
    std::cout << "Initial position (km): " << r0.transpose() << std::endl;
    std::cout << "Initial velocity (km/s): " << v0.transpose() << std::endl;

    // Orbital period
    f64 period = 2.0 * M_PI * std::sqrt(r0_mag * r0_mag * r0_mag / mu);
    f64 t0 = 0.;
    // f64 tf = 10;
    f64 tf = period;

    std::cout << "Orbital period: " << period << " s" << std::endl;

    // Number of steps for fixed integrator
    int Nsteps = 1000;

    // Test RK4 fixed-step integrator
    FixedStepIntegrator<decltype(f), vec6, RK4Policy> rk4(f, 0.1);
    std::vector<f64> times;
    std::vector<vec6> states;
    rk4.integrate(t0, tf, x0, times, states);

    // Final state
    vec6 xf = states.back();
    vec3 rf = xf.segment<3>(0);
    vec3 vf = xf.segment<3>(3);

    std::cout << "\nFixed RK4 Results\n";
    std::cout << "Final position (km): " << rf.transpose() << std::endl;
    std::cout << "Final velocity (km/s): " << vf.transpose() << std::endl;
    std::cout << "Radius error (km): " << (rf.norm() - r0.norm()) << std::endl;

    // std::cout << times.size() << std::endl;
    // for (int i = 0; i < times.size(); i++) {
    //     vec3 r_curr = states[i].segment(0, 3);
    //     std::cout << "Time (s): " << times[i]
    //               << ", Position: " << r_curr.transpose() << std::endl;
    // }

    // Export time & positions
    std::ofstream ofs("orbit.csv");
    ofs << "t_km,x_km,y_km,z_km\n";
    for (int i = 0; i < times.size(); ++i) {
        const auto &x6 = states[i];
        f64 t = times[i];
        f64 x = x6(0), y = x6(1), z = x6(2);
        ofs << t << "," << x << "," << y << "," << z << "\n";
    }
    ofs.close();
    std::cout << "Wrote orbit.csv (" << times.size() << " lines)\n";

    return 0;
}
