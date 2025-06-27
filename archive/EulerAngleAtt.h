#pragma once

#include "attitude.h"
#include "QuaternionAtt.h"

struct EulerAngleAtt : public Attitude {
    vec3 angles;                 // angles (radians)
    std::array<int, 3> sequence; // rotation sequence

    // constructors
    EulerAngleAtt() : angles(0, 0, 0), sequence({3, 2, 1}) {}
    EulerAngleAtt(
        const vec3 &angles,
        std::array<int, 3> sequence,
        UnitsAngle units_in
    )
        : sequence(sequence) {
        if (units_in != UnitsAngle::RADIANS) {
            for (int i = 0; i < 3; i++) {
                this->angles[i]
                    = convertAngle(angles[i], units_in, UnitsAngle::RADIANS);
            }
        }
    }

    // methods
    mat3 to_dcm() const override {
        mat3 R = mat3::Identity();
        for (int i = 2; i >= 0; --i) {
            R = R
                * single_axis_rotation(angles[i], RotationalAxis(sequence[i]));
        }
        return R;
    }
    void from_dcm(const mat3 &R) override {
        const int a = sequence[0];
        const int b = sequence[1];
        const int c = sequence[2];

        f64 &t1 = angles[0];
        f64 &t2 = angles[1];
        f64 &t3 = angles[2];

        const f64 eps = 1e-8;

        if (sequence == std::array<int, 3>{1, 2, 3}) {
            t2 = std::asin(-R(0, 2));
            if (std::abs(std::cos(t2)) > eps) {
                t1 = std::atan2(R(1, 2), R(2, 2));
                t3 = std::atan2(R(0, 1), R(0, 0));
            } else {
                t1 = 0;
                t3 = std::atan2(-R(1, 0), R(1, 1));
            }
        } else if (sequence == std::array<int, 3>{1, 3, 2}) {
            t2 = std::asin(R(0, 1));
            if (std::abs(std::cos(t2)) > eps) {
                t1 = std::atan2(-R(2, 1), R(1, 1));
                t3 = std::atan2(-R(0, 2), R(0, 0));
            } else {
                t1 = 0;
                t3 = std::atan2(R(2, 0), R(2, 2));
            }
        } else if (sequence == std::array<int, 3>{2, 1, 3}) {
            t2 = std::asin(-R(1, 0));
            if (std::abs(std::cos(t2)) > eps) {
                t1 = std::atan2(R(0, 0), R(2, 0));
                t3 = std::atan2(R(1, 2), R(1, 1));
            } else {
                t1 = 0;
                t3 = std::atan2(-R(0, 2), R(0, 1));
            }
        } else if (sequence == std::array<int, 3>{2, 3, 1}) {
            t2 = std::asin(R(1, 2));
            if (std::abs(std::cos(t2)) > eps) {
                t1 = std::atan2(-R(0, 2), R(2, 2));
                t3 = std::atan2(-R(1, 0), R(1, 1));
            } else {
                t1 = 0;
                t3 = std::atan2(R(0, 1), R(0, 0));
            }
        } else if (sequence == std::array<int, 3>{3, 1, 2}) {
            t2 = std::asin(-R(2, 1));
            if (std::abs(std::cos(t2)) > eps) {
                t1 = std::atan2(R(0, 1), R(1, 1));
                t3 = std::atan2(R(2, 0), R(2, 2));
            } else {
                t1 = 0;
                t3 = std::atan2(-R(0, 2), R(0, 0));
            }
        } else if (sequence == std::array<int, 3>{3, 2, 1}) {
            t2 = std::asin(-R(2, 0));
            if (std::abs(std::cos(t2)) > eps) {
                t1 = std::atan2(R(2, 1), R(2, 2));
                t3 = std::atan2(R(1, 0), R(0, 0));
            } else {
                t1 = 0;
                t3 = std::atan2(-R(0, 1), R(1, 1));
            }

            // Proper Euler sequences
        } else if (sequence == std::array<int, 3>{1, 2, 1}) {
            t2 = std::acos(R(0, 0));
            if (std::abs(std::sin(t2)) > eps) {
                t1 = std::atan2(R(1, 0), -R(2, 0));
                t3 = std::atan2(R(0, 1), R(0, 2));
            } else {
                t1 = 0;
                t3 = std::atan2(-R(2, 1), R(1, 1));
            }
        } else if (sequence == std::array<int, 3>{1, 3, 1}) {
            t2 = std::acos(R(0, 0));
            if (std::abs(std::sin(t2)) > eps) {
                t1 = std::atan2(R(2, 0), R(1, 0));
                t3 = std::atan2(R(0, 2), -R(0, 1));
            } else {
                t1 = 0;
                t3 = std::atan2(R(1, 2), R(2, 2));
            }
        } else if (sequence == std::array<int, 3>{2, 1, 2}) {
            t2 = std::acos(R(1, 1));
            if (std::abs(std::sin(t2)) > eps) {
                t1 = std::atan2(R(0, 1), -R(2, 1));
                t3 = std::atan2(R(1, 0), R(1, 2));
            } else {
                t1 = 0;
                t3 = std::atan2(-R(2, 0), R(0, 0));
            }
        } else if (sequence == std::array<int, 3>{2, 3, 2}) {
            t2 = std::acos(R(1, 1));
            if (std::abs(std::sin(t2)) > eps) {
                t1 = std::atan2(R(0, 1), R(2, 1));
                t3 = std::atan2(R(1, 2), -R(1, 0));
            } else {
                t1 = 0;
                t3 = std::atan2(R(2, 0), R(0, 0));
            }
        } else if (sequence == std::array<int, 3>{3, 1, 3}) {
            t2 = std::acos(R(2, 2));
            if (std::abs(std::sin(t2)) > eps) {
                t1 = std::atan2(R(0, 2), -R(1, 2));
                t3 = std::atan2(R(2, 0), R(2, 1));
            } else {
                t1 = 0;
                t3 = std::atan2(-R(1, 0), R(0, 0));
            }
        } else if (sequence == std::array<int, 3>{3, 2, 3}) {
            t2 = std::acos(R(2, 2));
            if (std::abs(std::sin(t2)) > eps) {
                t1 = std::atan2(R(1, 2), R(0, 2));
                t3 = std::atan2(R(2, 1), -R(2, 0));
            } else {
                t1 = 0;
                t3 = std::atan2(R(0, 1), R(1, 1));
            }
        } else {
            throw std::runtime_error("Unsupported Euler angle sequence");
        }
    }

    vec4 to_quaternion() const override {
        mat3 R = to_dcm();
        vec4 q = QuaternionAtt::dcm_to_quat(R);

        return vec4::Zero();
    }
    void from_quaternion(const vec4 &R) override {
        
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "Euler Angles [t_1, t_2, t_3]: [" << angles[0]
            << ", " << angles[1] << ", " << angles[2] << "]" << "\n"
            << "Sequence [s_1, s_2, s_3]: [" << sequence[0] << ", "
            << sequence[1] << ", " << sequence[2] << "]";
        return oss.str();
    }
};
