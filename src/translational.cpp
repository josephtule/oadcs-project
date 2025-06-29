#include "translational.h"
#include "typedefs.h"

vec6 gravity_newton(f64 t, vec6 x, f64 mu) {
    vec3 r = x.segment(0, 3);
    vec3 v = x.segment(3, 3);
    f64 r_mag = r.norm();

    vec3 g = -mu / (r_mag * r_mag * r_mag) * r;
    
    vec6 dxdt;
    dxdt << v, g;

    return dxdt;
}