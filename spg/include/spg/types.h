#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace spg
{
using Real = double;

using Int2 = std::array<int, 2>;
using Int3 = std::array<int, 3>;
using Int4 = std::array<int, 4>;

template <int DIM>
using Vector = Eigen::Matrix<Real, DIM, 1>;
using Vector2 = Eigen::Matrix<Real, 2, 1>;
using Vector3 = Eigen::Matrix<Real, 3, 1>;
using VectorX = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

template <typename T, int DIM>
using VectorT = Eigen::Matrix<T, DIM, 1>;
template <typename T>
using Vector3T = Eigen::Matrix<T, 3, 1>;
template <typename T>
using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <int DIM0, int DIM1 = DIM0>
using Matrix = Eigen::Matrix<Real, DIM0, DIM1>;
using Matrix2 = Eigen::Matrix<Real, 2, 2>;
using Matrix3 = Eigen::Matrix<Real, 3, 3>;
using Matrix4 = Eigen::Matrix<Real, 4, 4>;
using Matrix6 = Eigen::Matrix<Real, 6, 6>;
using MatrixX = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T, int DIM0, int DIM1 = DIM0>
using MatrixT = Eigen::Matrix<T, DIM0, DIM1>;
template <typename T>
using Matrix3T = Eigen::Matrix<T, 3, 3>;

using Triplet = Eigen::Triplet<Real>;
using SparseMatrix = Eigen::SparseMatrix<Real>;
}  // namespace spg