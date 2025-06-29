#pragma once

#include "typedefs.h"
enum class GravityModel {
    pointmass,
    zonal,
    spherical_harmonic,
    none,
};


vec6 gravity_newton(f64 t, vec6 x, f64 mu);