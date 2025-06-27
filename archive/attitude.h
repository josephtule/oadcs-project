#pragma once

#include <string>

#include "typedefs.h"
#include "units.h"

// enum struct AttitudeRepresentation {
//     quaternion,
//     euler_angle,
//     angle_axis,
//     rodrigues,
//     modified_rodrigues,
// };

struct Attitude {
    // AttitudeRepresentation att_rep = AttitudeRepresentation::quaternion;
    virtual ~Attitude() = default;

    // dcm conversions
    virtual mat3 to_dcm() const = 0;
    virtual void from_dcm(const mat3 &dcm) = 0;
    // NOTE: quaternions are the default attitude representation and
    // computations will be done with them
    // quaternions in the form [q1,q2,q3,q4] = [x,y,z,w]
    virtual vec4 to_quaternion() const = 0;
    virtual void from_quaternion(const vec4 &quat) = 0;

    virtual std::string to_string() const = 0; // for display
};

mat3 skew(vec3 v);

enum struct RotationalAxis { x = 1, y = 2, z = 3 };
inline mat3
single_axis_rotation(f64 angle, RotationalAxis axis, UnitsAngle units_in);
