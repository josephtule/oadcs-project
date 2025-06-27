#pragma once

#include <string>

#include "attitude.h"
#include "quaternionAtt.h"
#include "translational.h"
#include "typedefs.h"
#include "units.h"


struct Body {
    int id;
    std::string name;

    // state
    vec3 pos, vel; // translation
    QuaternionAtt quat;

    // units
    UnitsLinear u_linear = UnitsLinear::KILOMETER;
    UnitsAngle u_angle = UnitsAngle::RADIANS;
    UnitsMass u_mass = UnitsMass::KILOGRAMS;

    // gravity model
    GravityModel gravity_model = GravityModel::none;

    // flags
    bool update_position = true;
    bool update_attitude = false;

    // methods
    Body();
    Body(vec3 pos, vec3 vel, Attitude &attitude,std::string name) : pos(pos), vel(vel), name(name) {
        
    }
};