#pragma once

#include <string>

#include "typedefs.h"
#include "units.h"

struct Body {
    int id;
    std::string name;

    // state
    vec3 pos, vel; // translation
    vec4 ep;       // euler param/quaternions
    vec3 omega;    // angular velocity

    // units
    UnitsLinear u_linear = UnitsLinear::KILOMETER;
    UnitsAngle u_angle = UnitsAngle::RADIANS;
    UnitsMass u_mass = UnitsMass::KILOGRAMS;

    // gravity model
    // GravityModel gravity_model = GravityModel::none;

    // flags
    bool update_position = true;
    bool update_attitude = false;

    // methods
    Body();
};