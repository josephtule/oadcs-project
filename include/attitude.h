#pragma once

#include "typedefs.h"
#include "units.h"
#include <array>
#include <vector>

using std::array;
using std::vector;

enum struct RotationalAxis { x = 1, y = 2, z = 3 };

inline mat3
single_axis_rotation(f64 angle, RotationalAxis axis, UnitsAngle units_in);

inline mat3 skew(vec3 v);

// euler param (quaternion)
inline mat3 ep_to_dcm(const vec4 &ep);
inline vec4 dcm_to_ep(const mat3 &R);
inline vec4 ep_kde(vec4 &ep, vec3 &omega);
inline vec3 omega_ode(vec3 omega, vec3 torque, mat3 I, mat3 I_inv);

// euler angles
inline mat3
ea_to_dcm(const vec3 &angles, const vec3 &sequence, UnitsAngle units_in);
inline vec4
ea_to_ep(const vec3 &angles, const vec3 &sequence, UnitsAngle units_in);

// principle rotation
inline mat3 pr_to_dcm(const vec3 &axis, const f64 &angle, UnitsAngle units_in);
inline vec4 pr_to_ep(const vec3 &axis, const f64 &angle, UnitsAngle units_in);

// classical rodrigues param
inline mat3 crp_to_dcm(const vec3 &crp);
inline vec4 crp_to_ep(const vec3 &crp);

// modified rodrigues param
inline mat3 mrp_to_dcm(const vec3 &crp);
inline vec4 mrp_to_ep(const vec3 &crp);
