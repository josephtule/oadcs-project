#include "Attitude.h"
#include "typedefs.h"
#include "units.h"

inline mat3 skew(vec3 v) {
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

// euler params
inline vec4 dcm_to_ep(const mat3 &R) {
    vec4 q;

    f64 tr = R.trace();

    if (tr > 0.0) {
        f64 S = std::sqrt(tr + 1.0) * 2.0; // 4*q4
        q[3] = 0.25 * S;
        q[0] = (R(2, 1) - R(1, 2)) / S;
        q[1] = (R(0, 2) - R(2, 0)) / S;
        q[2] = (R(1, 0) - R(0, 1)) / S;
    } else if (R(0, 0) > R(1, 1) && R(0, 0) > R(2, 2)) {
        f64 S = std::sqrt(1.0 + R(0, 0) - R(1, 1) - R(2, 2)) * 2.0; // 4*q1
        q[0] = 0.25 * S;
        q[1] = (R(0, 1) + R(1, 0)) / S;
        q[2] = (R(0, 2) + R(2, 0)) / S;
        q[3] = (R(2, 1) - R(1, 2)) / S;
    } else if (R(1, 1) > R(2, 2)) {
        f64 S = std::sqrt(1.0 + R(1, 1) - R(0, 0) - R(2, 2)) * 2.0; // 4*q2
        q[0] = (R(0, 1) + R(1, 0)) / S;
        q[1] = 0.25 * S;
        q[2] = (R(1, 2) + R(2, 1)) / S;
        q[3] = (R(0, 2) - R(2, 0)) / S;
    } else {
        f64 S = std::sqrt(1.0 + R(2, 2) - R(0, 0) - R(1, 1)) * 2.0; // 4*q3
        q[0] = (R(0, 2) + R(2, 0)) / S;
        q[1] = (R(1, 2) + R(2, 1)) / S;
        q[2] = 0.25 * S;
        q[3] = (R(1, 0) - R(0, 1)) / S;
    }

    return q.normalized();
}
inline mat3 ep_to_dcm(const vec4 &ep) {
    const f64 q1 = ep.x();
    const f64 q2 = ep.y();
    const f64 q3 = ep.z();
    const f64 q4 = ep.w();

    const f64 q1s = q1 * q1;
    const f64 q2s = q2 * q2;
    const f64 q3s = q3 * q3;
    const f64 q4s = q4 * q4;

    const vec3 q_13 = ep.segment(0, 3);
    const f64 q_13_norm = q_13.norm();

    mat3 eye3 = mat3::Identity();

    // mat3 R_ = (q4 * q4 - q_13_norm * q_13_norm) * eye3 - 2 * q4 *
    // skew(q_13)
    //           + 2 * q_13 * q_13.transpose();
    mat3 R;
    R <<
        // row 0
        q1s - q2s - q3s - q4s,
        2 * (q1 * q2 + q3 + q4), 2 * (q1 * q3 - q2 * q4),
        // row 1
        2 * (q2 * q1 - q3 * q4), -q1s + q2s - q3s + q4s,
        2 * (q2 * q3 + q1 * q4),
        // row 2
        2 * (q3 * q1 + q2 * q4), 2 * (q3 * q2 - q1 * q4),
        -q1s - q2s + q3s + q4s;
    return R;
}
inline vec4 ep_kde(vec4 &ep, vec3 &omega) {

    f64 ep1 = ep(0);
    f64 ep2 = ep(1);
    f64 ep3 = ep(2);
    f64 ep4 = ep(3);
    f64 o1 = omega(0);
    f64 o2 = omega(1);
    f64 o3 = omega(2);

    vec4 depdt = 0.5
                 * vec4(
                     o3 * ep2 - o2 * ep3 + o1 * ep[3],
                     -o3 * ep1 + o1 * ep3 + o2 * ep[3],
                     o2 * ep1 - o1 * ep2 + o3 * ep[3],
                     -o1 * ep1 - o2 * ep2 - o3 * ep3
                 );

    return depdt;
}
inline vec3
omega_ode(vec3 omega, vec3 torque, mat3 I, mat3 I_inv, bool I_diag = true) {
    vec3 domegadt;
    f64 w1 = omega(0);
    f64 w2 = omega(1);
    f64 w3 = omega(2);
    if (I_diag) {
        domegadt = vec3(
            (I(1, 1) - I(2, 2)) / I(0, 0) * w2 * w3,
            (I(2, 2) - I(0, 0)) / I(1, 1) * w1 * w3,
            (I(0, 0) - I(1, 1)) / I(2, 2) * w1 * w2
        );
        domegadt += vec3(
            torque(0) * I_inv(0, 0),
            torque(1) * I_inv(1, 1),
            torque(2) * I_inv(2, 2)
        );
    } else {
        // angular momentum
        vec3 h = vec3(
            I(0, 0) * w1 + I(0, 1) * w2 + I(0, 2) * w3,
            I(1, 0) * w1 + I(1, 1) * w2 + I(1, 2) * w3,
            I(2, 0) * w1 + I(2, 1) * w2 + I(2, 2) * w3
        );
        vec3 dhdt_tf = vec3(
            h(2) * w2 - h(1) * w3, h(0) * w3 - h(2) * w1, h(1) * w1 - h(0) * w2
        );
        vec3 dhdt = torque - dhdt_tf;
        // domegadt = I_inv * (torque - dhdt_tf);
        domegadt = vec3(
            I_inv(0, 0) * dhdt(0) + I_inv(0, 1) * dhdt(1)
                + I_inv(0, 2) * dhdt(2),
            I_inv(1, 0) * dhdt(0) + I_inv(1, 1) * dhdt(1)
                + I_inv(1, 2) * dhdt(2),
            I_inv(2, 0) * dhdt(0) + I_inv(2, 1) * dhdt(1)
                + I_inv(2, 2) * dhdt(2)
        );
    }

    return domegadt;
}

// euler angles
inline mat3 ea_to_dcm(
    const vec3 &angles,
    const vec3 &sequence,
    UnitsAngle units_in = UnitsAngle::RADIANS
) {
    vec3 angles_ = angles;
    if (units_in != UnitsAngle::RADIANS) {
        for (int i = 0; i < 3; i++) {
            angles_(i) = convertAngle(angles(i), units_in, UnitsAngle::RADIANS);
        }
    }

    mat3 R = mat3::Identity();
    for (int i = 2; i >= 0; --i) {
        R = R * single_axis_rotation(angles_(i), RotationalAxis(sequence(i)));
    }
    return R;
}
inline vec4 ea_to_ep(
    const vec3 &angles,
    const vec3 &sequence,
    UnitsAngle units_in = UnitsAngle::RADIANS
) {
    mat3 R = ea_to_dcm(angles, sequence, units_in);
    return dcm_to_ep(R);
}

// principle rotation
inline mat3 pr_to_dcm(
    const vec3 &axis,
    const f64 &angle,
    UnitsAngle units_in = UnitsAngle::RADIANS
) {
    f64 theta = angle;
    if (units_in != UnitsAngle::RADIANS) {
        theta = convertAngle(angle, units_in, UnitsAngle::RADIANS);
    }

    f64 sigma = 1 - cos(angle);
    f64 c = cos(theta);
    f64 s = sin(theta);
    f64 e1 = axis(0);
    f64 e2 = axis(1);
    f64 e3 = axis(2);

    mat3 R;
    R <<
        // row 1
        e1 * e1 * sigma + c,
        e1 * e2 * sigma + e3 * s, e1 * e3 * sigma - e2 * s,
        // row 2
        e2 * e1 * sigma - e3 * s, e2 * e2 * sigma + c, e2 * e3 * sigma + e1 * s,
        // row 3
        e3 * e1 * sigma + e2 * s, e3 * e2 * sigma - e1 * s, e3 * e3 * sigma + c;
    return R;
}
inline vec4 pr_to_ep(
    const vec3 &axis,
    const f64 &angle,
    UnitsAngle units_in = UnitsAngle::RADIANS
) {
    mat3 R = pr_to_dcm(axis, angle, units_in);
    return dcm_to_ep(R);
}

// classical rodrigues param
inline mat3 crp_to_dcm(const vec3 &crp) {
    f64 q1 = crp(0);
    f64 q2 = crp(1);
    f64 q3 = crp(2);
    f64 q1s = q1 * q1;
    f64 q2s = q2 * q2;
    f64 q3s = q3 * q3;
    f64 qdotq = crp.dot(crp);
    mat3 R;
    R <<
        // row 1
        1 + q1s - q2s - q3s,
        2 * (q1 * q2 + q3), 2 * (q1 * q3 - q2),
        // row 2
        2 * (q2 * q1 - q3), 1 - q1s + q2s - q3s, 2 * (q2 * q3 + q1),
        // row 3
        2 * (q3 * q1 + q2), 2 * (q3 * q2 - q1), 1 - q1s - q2s + q3s;
    R = R * 1 / (1 + qdotq);

    return R;
}
inline vec4 crp_to_ep(const vec3 &crp) {
    mat3 R = crp_to_dcm(crp);
    return dcm_to_ep(R);
}

// modified rodrigues param
inline mat3 mrp_to_dcm(const vec3 &mrp) {
    f64 s1 = mrp(0);
    f64 s2 = mrp(1);
    f64 s3 = mrp(2);
    f64 s1s = s1 * s1;
    f64 s2s = s2 * s2;
    f64 s3s = s3 * s3;
    f64 sdots = mrp.dot(mrp);
    f64 denom = (1 + sdots) * (1 + sdots);
    f64 oms2 = 1 - sdots;
    mat3 R;
    R <<
        // row 1
        4 * (s1s - s2s - s3s) + oms2 * oms2,
        8 * s1 * s2 + 4 * s3 * oms2, 8 * s1 * s3 - 4 * s2 * oms2,
        // row 2
        8 * s2 * s1 - 4 * s3 * oms2, 4 * (-s1s + s2s - s3s) + oms2 * oms2,
        8 * s2 * s3 + 4 * s1 * oms2,
        // row 3
        8 * s3 * s1 + 4 * s2 * oms2, 8 * s3 * s2 - 4 * s1 * oms2,
        4 * (-s1s - s2s + s3s) + oms2 * oms2;
    R = R / denom;
    return R;
}
inline vec4 mrp_to_ep(const vec3 &mrp) {
    mat3 R = mrp_to_dcm(mrp);
    return dcm_to_ep(R);
}
