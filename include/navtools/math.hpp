/**
|========================================== math.hpp ==============================================|
|                                                                                                  |
|   @file     include/navtools/math.hpp                                                            |
|   @brief    Common mathematical operations.                                                      |
|   @ref      Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems          |
|               - (2013) Paul D. Groves                                                            |
|   @date     July 2024                                                                            |
|                                                                                                  |
|==================================================================================================|
*/

// TODO: move conversions

#ifndef NAVTOOLS_MATH_HPP
#define NAVTOOLS_MATH_HPP

#include <Eigen/Dense>
#include <cmath>

#include "navtools/constants.hpp"

namespace navtools {

//* ===== Skew-Symmmetric Matrix =============================================================== *//

//! === SKEW ===
/// @brief      Converts vector into its skew symmetric form
/// @param v    3x1 input vector
/// @param M    3x3 skew symmetric matrix
/// @returns    3x3 skew symmetric form of v
template <typename T = double>
void Skew(Eigen::Matrix<T, 3, 3> &M, const Eigen::Vector<T, 3> &v) {
  M << 0.0, -v(2), v(1), v(2), 0.0, -v(0), -v(1), v(0), 0.0;
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> Skew(const Eigen::Vector<T, 3> &v) {
  Eigen::Matrix<T, 3, 3> M;
  Skew(M, v);
  return M;
}

//! === DESKEW ===
/// @brief      Converts skew symmetric form into its vector
/// @param M    3x3 input skew symmetric matrix
/// @param v    3x1 vector
/// @returns    3x1 vector form of M
template <typename T = double>
void DeSkew(Eigen::Vector<T, 3> &v, const Eigen::Matrix<T, 3, 3> &M) {
  v << M(2, 1), M(0, 2), M(1, 0);
}
template <typename T = double>
Eigen::Vector<T, 3> DeSkew(const Eigen::Matrix<T, 3, 3> &M) {
  Eigen::Vector<T, 3> v;
  DeSkew(v, M);
  return v;
}

//* ===== Sorting ============================================================================== *//

//! === MINPOSITIVEIDX ===
/// @brief find the index of the minimum value in vector greater than 0
template <typename T = double>
size_t MinPositiveIdx(const std::vector<T> &arr, const T &set_point = 0.0) {
  size_t result = 0;
  T tmp;
  T max_val = std::numeric_limits<T>::max();
  for (size_t i = 0; i < arr.size(); i++) {
    tmp = arr[i] - set_point;
    if ((0 < tmp) && (tmp < max_val)) {
      max_val = tmp;
      result = i;
    }
  }
  return result;
}

//* ===== Modulus Operations ===================================================================
//*//
//! === CIRCFMOD ===
/// @brief      Modulus of floating point number
/// @param x    user input and output
/// @param y    value to take modulus about
/// @returns    modulus of number
template <typename T>
constexpr void CircMod(T &x, const T y) {
  x -= std::floor(x / y) * y;
}
// template <typename T>
// constexpr T CircMod(T x, const T y) {
//     return x - std::floor(x / y) * y;
// }

//* ===== Pi Wrapping ==========================================================================
//*//

//! === CIRCMOD2PI ===
/// @brief      Wraps angles from [0, 2*pi]
/// @param x    Inputs containing angles in radians
/// @returns    Wrapped/normalized angles [radians]
template <typename T = double>
void CircMod2Pi(T &x) {
  CircMod(x, TWO_PI<T>);
}

//! === WRAPPITOPI ===
/// @brief      Wraps angles from [-pi, pi]
/// @param x    Inputs containing angles in radians
/// @returns    Wrapped/normalized angles [radians]
template <typename T = double>
void WrapPiToPi(T &x) {
  CircMod2Pi(x);
  if (x > PI<T>) {
    x -= TWO_PI<T>;
  }
}

//! === WRAPEULERANGLES ===
/// @brief      Auto wrap euler angles depending on pitch angle
/// @param x    Euler angles [radians]
/// @returns    Correctly wrapped Euler angles
template <typename T = double>
void WrapEulerAngles(Eigen::Vector<T, 3> &x) {
  if (x(1) > HALF_PI<T>) {
    x(0) += PI<T>;
    x(1) = PI<T> - x(1);
    x(2) += PI<T>;
  } else if (x(1) < -HALF_PI<T>) {
    x(0) += PI<T>;
    x(1) = -PI<T> - x(1);
    x(2) += PI<T>;
  }
  WrapPiToPi(x(0));
  WrapPiToPi(x(2));
}

//! === RAD2DEG ===
/// @brief      Convert radians to degrees
/// @param  x   radians
/// @returns    degrees
template <typename T>
void rad2deg(T &x) {
  x *= RAD2DEG<T>;
}
template <typename T>
constexpr T rad2deg(const T &x) {
  return x * RAD2DEG<T>;
}

//! === DEG2RAD ===
/// @brief      Convert degrees to radians
/// @param  x   degrees
/// @returns    radiasn
template <typename T>
void deg2rad(T &x) {
  x *= DEG2RAD<T>;
}
template <typename T>
constexpr T deg2rad(const T &x) {
  return x * DEG2RAD<T>;
}

//* ===== Matrix Products/Normalization ========================================================
//*//

//! === QUATMAT ===
/// @brief      convert quaternion into its 4x4 matrix view
template <typename T>
Eigen::Matrix<T, 4, 4> quatmat(const Eigen::Vector<T, 4> &q) {
  return Eigen::Matrix<T, 4, 4>{
      {q(0), -q(1), -q(2), -q(3)},
      {q(1), q(0), -q(3), q(2)},
      {q(2), q(3), q(0), -q(1)},
      {q(3), -q(2), q(1), q(0)}};
}

//! === QUATDOT ===
/// @brief      quaternion product
template <typename T = double>
Eigen::Vector<T, 4> quatdot(const Eigen::Vector<T, 4> &p, const Eigen::Vector<T, 4> &q) {
  // p o q
  return Eigen::Vector<T, 4>{
      p(0) * q(0) - p(1) * q(1) - p(2) * q(2) - p(3) * q(3),
      p(0) * q(1) + p(1) * q(0) + p(2) * q(3) - p(3) * q(2),
      p(0) * q(2) - p(1) * q(3) + p(2) * q(0) + p(3) * q(1),
      p(0) * q(3) + p(1) * q(2) - p(2) * q(1) + p(3) * q(0)};
}
template <typename T = double>
Eigen::Vector<T, 4> quatdot(const Eigen::Vector<T, 3> &a, const Eigen::Vector<T, 4> &q) {
  // a o q
  return Eigen::Vector<T, 4>{
      -a(0) * q(1) - a(1) * q(2) - a(2) * q(3),
      a(0) * q(0) + a(1) * q(3) - a(2) * q(2),
      -a(0) * q(3) + a(1) * q(0) + a(2) * q(1),
      a(0) * q(2) - a(1) * q(1) + a(2) * q(0)};
}
template <typename T = double>
Eigen::Vector<T, 4> quatdot(const Eigen::Vector<T, 4> &p, const Eigen::Vector<T, 3> &a) {
  // p o a
  return Eigen::Vector<T, 4>{
      -p(1) * a(0) - p(2) * a(1) - p(3) * a(2),
      p(0) * a(0) + p(2) * a(2) - p(3) * a(1),
      p(0) * a(1) - p(1) * a(2) + p(3) * a(0),
      p(0) * a(2) + p(1) * a(1) - p(2) * a(0)};
}

//! === QUATCONJ/QUATINV ===
/// @brief      quaternion conjugate/inverse
template <typename T = double>
void quatconj(Eigen::Vector<T, 4> &q) {
  q(1) *= -1.0;
  q(2) *= -1.0;
  q(3) *= -1.0;
}
template <typename T = double>
void quatinv(Eigen::Vector<T, 4> &q) {
  q(1) *= -1.0;
  q(2) *= -1.0;
  q(3) *= -1.0;
}

//! === QUATNORM ===
/// @brief      quaternion normalization
template <typename T = double>
void quatnorm(Eigen::Vector<T, 4> &q) {
  q /= q.norm();
}

//! === DCMNORM ===
/// @brief      DCM normalization
template <typename T = double>
void dcmnorm(Eigen::Matrix<T, 3, 3> &R) {
  // Groves 5.79, 5.80
  Eigen::Vector<T, 3> c1 = R.col(0);
  Eigen::Vector<T, 3> c2 = R.col(1);
  Eigen::Vector<T, 3> c3 = R.col(2);

  T c1c2 = R(0, 0) * R(0, 1) + R(1, 0) * R(1, 1) + R(2, 0) * R(2, 1);
  T c1c3 = R(0, 0) * R(0, 2) + R(1, 0) * R(1, 2) + R(2, 0) * R(2, 2);
  T c2c3 = R(0, 1) * R(0, 2) + R(1, 1) * R(1, 2) + R(2, 1) * R(2, 2);

  c1 -= 0.5 * (c1c2 * c2 + c1c3 * c3);
  c2 -= 0.5 * (c1c2 * c1 + c2c3 * c3);
  c3 -= 0.5 * (c1c3 * c1 + c2c3 * c2);
  c1 /= c1.norm();
  c2 /= c2.norm();
  c3 /= c3.norm();

  R.col(0) = c1;
  R.col(1) = c2;
  R.col(2) = c3;
}

//* ===== Matrix Exponential ===================================================================
//*//

//! === RODRIGUES ===
/// @brief      Rodrigues formula for the approximation of a matrix exponential
/// @param vec          size 3 vector
/// @param vec_norm     2-norm of vec
/// @returns    matrix exponential
template <typename T = double>
Eigen::Matrix<T, 3, 3> Rodrigues(const Eigen::Vector<T, 3> &vec) {
  T vec_norm = vec.norm();
  Eigen::Matrix<T, 3, 3> skew_sym = Skew(vec / vec_norm);
  return Eigen::Matrix<T, 3, 3>::Identity() + (std::sin(vec_norm) * skew_sym) +
         ((1.0 - std::cos(vec_norm)) * skew_sym * skew_sym);
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> Rodrigues(const Eigen::Vector<T, 3> &vec, const T &vec_norm) {
  Eigen::Matrix<T, 3, 3> skew_sym = Skew<T>(vec / vec_norm);
  return Eigen::Matrix<T, 3, 3>::Identity() + (std::sin(vec_norm) * skew_sym) +
         ((1.0 - std::cos(vec_norm)) * skew_sym * skew_sym);
}

//! === RODRIGUES4 ===
/// @brief      Rodrigues formula for the 4-th order approximation of a matrix exponential
/// @param vec          size 3 vector
/// @param vec_norm     2-norm of vec
/// @returns    matrix exponential
template <typename T = double>
Eigen::Matrix<T, 3, 3> Rodrigues4(const Eigen::Vector<T, 3> &vec) {
  T vec_norm = vec.norm();
  Eigen::Matrix<T, 3, 3> skew_sym = Skew(vec);
  T norm_squared = vec_norm * vec_norm;
  return Eigen::Matrix<T, 3, 3>::Identity() + ((1.0 - (norm_squared / 6.0)) * skew_sym) +
         ((0.5 - (norm_squared / 24.0)) * skew_sym * skew_sym);
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> Rodrigues4(const Eigen::Vector<T, 3> &vec, const T &vec_norm) {
  Eigen::Matrix<T, 3, 3> skew_sym = Skew(vec);
  T norm_squared = vec_norm * vec_norm;
  return Eigen::Matrix<T, 3, 3>::Identity() + ((1.0 - (norm_squared / 6.0)) * skew_sym) +
         ((0.5 - (norm_squared / 24.0)) * skew_sym * skew_sym);
}

//! === VEC2EXPM ===
/// @brief      Converts vector into its matrix exponential approximation form
/// @param vec  size 3 vector
/// @returns    3x3 matrix exponential
template <typename T = double>
Eigen::Matrix<T, 3, 3> vec2expm(const Eigen::Vector<T, 3> &vec) {
  T vec_norm = vec.norm();
  if (vec_norm < 0.02) {
    return Rodrigues4(vec, vec_norm);
  } else {
    return Rodrigues(vec, vec_norm);
  }
}

//! === EXPM2VEC ===
/// @brief      Converts matrix exponential approximation into its vector form
/// @param mat  3x3 matrix exponential
/// @returns    size 3 vector
template <typename T = double>
Eigen::Vector<T, 3> expm2vec(const Eigen::Matrix<T, 3, 3> &mat) {
  T phi = std::acos((mat.trace() - 1.0) / 2.0);
  if (phi == 0.0) {
    return Eigen::Vector<T, 3>::Zero();
  }
  return phi * DeSkew(mat - mat.transpose()) / (2.0 * std::sin(phi));
}

//* ===== Signal to Noise ======================================================================
//*//
// https://insidegnss.com/wp-content/uploads/2018/01/novdec10-Solutions.pdf

//! === WATT2DB ===
/// @brief      Convert unit of power to decibels
/// @param w    Input power
/// @param db   Output decibels
/// @returns    Decibels
template <typename T = double>
void watt2db(T &db, const T &w) {
  db = 10.0 * std::log10(w);
}
template <typename T = double>
T watt2db(const T &w) {
  T db;
  watt2db(db, w);
  return db;
}

//! === DB2WATT ===
/// @brief      Convert decibels to unit of power
/// @param db   Input decibels
/// @param w    Output power
/// @returns    unit of power
template <typename T = double>
void db2watt(T &w, const T &db) {
  w = std::pow(10.0, db / 10.0);
}
template <typename T = double>
T db2watt(const T &db) {
  T w;
  db2watt(w, db);
  return w;
}

//! === VOLT2DB ===
/// @brief      Convert unit of sqrt-power (volts) to decibels
/// @param v    Input sqrt-power
/// @param dB   Output decibels
/// @returns    Decibels
template <typename T = double>
void volt2db(T &db, const T &v) {
  db = 20.0 * std::log10(v);
}
template <typename T = double>
T volt2db(const T &v) {
  T db;
  volt2db(db, v);
  return db;
}

//! === DB2VOLT ===
/// @brief      Convert decibels to unit of sqrt-power
/// @param db   Input decibels
/// @param v    Output sqrt-power
/// @returns    unit of sqrt-power
template <typename T = double>
void db2volt(T &v, const T db) {
  v = std::pow(10.0, db / 20.0);
}
template <typename T = double>
T db2volt(const T &db) {
  T v;
  db2volt(v, db);
  return v;
}

//! === CN02SNR ===
/// @brief      Convert Carrier-to-Noise Ratio into raw Signal-to-Noise Ratio
/// @param cn0      Carrier-to-noise density ratio [dB/Hz]
/// @param fe_bw    Receiver front-end bandwidth [Hz]
/// @param temp     Additional noise figure/loses from the given temperature [K]
/// @param eta      Additional noise figure/loses [dB]
/// @param snr      Signal to noise ratio [dB]
template <typename T = double>
void cn02snr(T &snr, const T &cn0, const T &fe_bw, const T &temp = 290.0, const T &eta = 0.0) {
  snr = cn0 - watt2db(fe_bw) - watt2db(1.0 + temp / 290.0) - eta;
}
template <typename T = double>
T cn02snr(const T &cn0, const T &fe_bw, const T &temp = 290.0, const T &eta = 0.0) {
  T snr;
  cn02snr(snr, cn0, fe_bw, temp, eta);
  return snr;
}

//! === SNR2CN0 ===
/// @brief      Convert Carrier-to-Noise Ratio into raw Signal-to-Noise Ratio
/// @param snr      Signal to noise ratio [dB]
/// @param fe_bw    Receiver front-end bandwidth [Hz]
/// @param temp     Additional noise figure/loses from the given temperature [K]
/// @param eta      Additional noise figure/loses [dB]
/// @param cn0      Carrier-to-noise density ratio [dB/Hz]
template <typename T = double>
void snr2cn0(T &cn0, const T &snr, const T &fe_bw, const T &temp = 290.0, const T &eta = 0.0) {
  cn0 = snr + watt2db(fe_bw) + watt2db(1.0 + temp / 290.0) + eta;
}
template <typename T = double>
T snr2cn0(const T &snr, const T &fe_bw, const T &temp = 290.0, const T &eta = 0.0) {
  T cn0;
  snr2cn0(cn0, snr, fe_bw, temp, eta);
  return cn0;
}

}  // namespace navtools

#endif
