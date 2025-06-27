#pragma once

#include "QuaternionAtt.h"
#include "attitude.h"


struct DirectionCosineAtt : public Attitude {

    mat3 dcm;

    // identity
    DirectionCosineAtt() : dcm(mat3::Identity()) {}
    DirectionCosineAtt(const mat3 &R_init) { from_dcm(R_init); }

    mat3 to_dcm() const override { return dcm; }
    void from_dcm(const mat3 R_init) {
        // ensure orthonormal
        Eigen::JacobiSVD<mat3> svd(
            R_init, Eigen::ComputeFullU | Eigen::ComputeFullV
        );
        dcm = svd.matrixU() * svd.matrixV().transpose();
    }

    vec4 to_quaternion() const override {
        return QuaternionAtt::dcm_to_quat(dcm);
    }
    void from_quaternion(const vec4 &q) override {
        dcm = QuaternionAtt::quat_to_dcm(q);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "DCM:\n" << dcm;
        return oss.str();
    }
};