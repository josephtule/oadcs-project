#pragma once

#include "attitude.h"
#include <iostream>

struct QuaternionAtt : public Attitude {

    vec4 quat; // [x,y,z,w]
    vec3 omega; // angular velocity

    // identity quaternion/euler parameter
    QuaternionAtt() : quat(vec4(0, 0, 0, 1)) {}
    QuaternionAtt(const vec4 &quat_init) { from_quaternion(quat_init); }

    vec4 to_quaternion() const override { return quat; }
    void from_quaternion(const vec4 &quat_new) override {
        // quat = quat_new.normalized();
        quat = quat_new;
    }

    static vec4 dcm_to_quat(const mat3 &R) {
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
    static mat3 quat_to_dcm(const vec4 &q) {
        const f64 q1 = q.x();
        const f64 q2 = q.y();
        const f64 q3 = q.z();
        const f64 q4 = q.w();

        const f64 q1s = q1 * q1;
        const f64 q2s = q2 * q2;
        const f64 q3s = q3 * q3;
        const f64 q4s = q4 * q4;

        const vec3 q_13 = q.segment(0, 3);
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

    mat3 to_dcm() const override {
        return quat_to_dcm(quat);
    }
    void from_dcm(const mat3 &R) override { quat = dcm_to_quat(R); }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "Quat [x,y,z,w]: [" << quat[0] << ", " << quat[1] << ", "
            << quat[2] << ", " << quat[3] << "]";
        return oss.str();
    }
};