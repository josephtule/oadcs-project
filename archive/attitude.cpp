#include "attitude.h"

mat3 skew(vec3 v) {
    mat3 s;
    s << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;

    return s;
}

inline mat3 single_axis_rotation(
    f64 angle,
    RotationalAxis axis,
    UnitsAngle units_in = UnitsAngle::RADIANS
) {

    mat3 R = mat3::Identity();

    if (units_in != UnitsAngle::RADIANS) {
        angle = convertAngle(angle, units_in, UnitsAngle::RADIANS);
    }

    f64 c = cos(angle);
    f64 s = sin(angle);

    switch (axis) {
    case RotationalAxis::x: {
        R << 1, 0, 0, 0, c, -s, 0, s, c;
        break;
    }
    case RotationalAxis::y: {
        R << c, 0, s, 0, 1, 0, -s, 0, c;
        break;
    }
    case RotationalAxis::z: {
        R << c, -s, 0, s, c, 0, 0, 0, 1;
        break;
    }
    }

    return R;
}
