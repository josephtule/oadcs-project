#pragma once

#include "Eigen/Geometry"
#include <Eigen/Dense>
#include <cstdint>
#include <vector>

namespace eig = Eigen;

// Type aliases
using f32 = float;
using f64 = double;
using i8 = int8_t;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;
using uint = unsigned int;
using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

template <class T> using vec = std::vector<T>;

// Eigen types
// Vectors
using vec2 = eig::Vector2d;
using vec3 = eig::Vector3d;
using vec4 = eig::Vector4d;
using vec5 = eig::Vector<f64, 5>;
using vec6 = eig::Vector<f64, 6>;
using vec7 = eig::Vector<f64, 7>;
using vec8 = eig::Vector<f64, 8>;
using vec9 = eig::Vector<f64, 9>;
using vec10 = eig::Vector<f64, 10>;
using vec11 = eig::Vector<f64, 11>;
using vec12 = eig::Vector<f64, 12>;
using vecx = eig::VectorXd;
// Matrices
using mat2 = eig::Matrix2d;
using mat3 = eig::Matrix3d;
using mat4 = eig::Matrix4d;
using mat5 = eig::Matrix<f64, 5, 5>;
using mat6 = eig::Matrix<f64, 6, 6>;
using mat7 = eig::Matrix<f64, 7, 7>;
using mat8 = eig::Matrix<f64, 8, 8>;
using mat9 = eig::Matrix<f64, 9, 9>;
using mat10 = eig::Matrix<f64, 10, 10>;
using mat11 = eig::Matrix<f64, 11, 11>;
using mat12 = eig::Matrix<f64, 12, 12>;
using matx = eig::Matrix<f64, eig::Dynamic, eig::Dynamic>;
// Quaternions
using quate = eig::Quaterniond;

