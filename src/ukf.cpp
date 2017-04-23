#include "ukf.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

namespace {
///* Process noise standard deviation longitudinal acceleration in m/s^2
const double kStdA{30};

///* Process noise standard deviation yaw acceleration in rad/s^2
const double kStdYawdd{30};

///* Laser measurement noise standard deviation position1 in m
const double kStdLaspx{0.15};

///* Laser measurement noise standard deviation position2 in m
const double kStdLaspy{0.15};

///* Radar measurement noise standard deviation radius in m
const double kStdRadr{0.3};

///* Radar measurement noise standard deviation angle in rad
const double kStdRadphi{0.03};

///* Radar measurement noise standard deviation radius change in m/s
const double kStdRadrd{0.3};

// UKF Hyperparameters, described in 10.5.1 Sigma Point Computation chapter
// of Kalman and Bayesian Filters in Python by Roger R Labbe Jr book,
// January 14, 2017 edition

const double kAlpha{1};

const double kBeta{2.0};

double getKappa(int n) {
  // Set to obtain positive W0Mean and W0Covariance
  return 10 - n;
}

double getLambda(int n) {
  return kAlpha * kAlpha * (n + getKappa(n)) - n;
}

double getW0Mean(int n) {
  const double lambda = getLambda(n);
  return lambda / (lambda + n);
}

double getW0Covariance(int n) {
  return getW0Mean(n) + 1 - kAlpha * kAlpha + kBeta;
}

double getW1ThroughN(int n) {
  const double lambda = getLambda(n);
  return 0.5 / (lambda + n);
}

/**
 * Angle normalization
 */
double normPi(double angleRad) {
  angleRad = std::fmod(angleRad, 2 * M_PI);
  if (angleRad > M_PI) {
    angleRad -= 2 * M_PI;
  }

  if (angleRad < -M_PI) {
    angleRad += 2 * M_PI;
  }

  return angleRad;
}
}  // namespace

UKF::UKF()
  : x_(static_cast<int>(UKF::kNx))
  , P_(static_cast<int>(UKF::kNx), static_cast<int>(UKF::kNx))
  , Xsig_pred_(static_cast<int>(UKF::kNxAug),
               static_cast<int>(UKF::kNxAug * 2 + 1)) {
  x_.setZero();
  P_.setIdentity();
  P_ *= 100000;
  Xsig_pred_.setZero();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage& meas_package) {
  if (!is_initialized_) {
    switch (meas_package.sensor_type_) {
      case MeasurementPackage::LASER: {
        x_.head(2) = meas_package.raw_measurements_;
        P_(kPosX, kPosX) = kStdLaspx * kStdLaspx;
        P_(kPosY, kPosY) = kStdLaspy * kStdLaspy;
        break;
      }

      case MeasurementPackage::RADAR: {
        const double ro = meas_package.raw_measurements_[0];
        const double phi = meas_package.raw_measurements_[1];
        const double ro_dot = meas_package.raw_measurements_[2];

        x_(kPosX) = ro * std::cos(phi);
        x_(kPosY) = ro * std::sin(phi);
        x_(kVelocity) = ro_dot;
        x_(kYaw) = phi;
        x_(kYawRate) = 0.0;

        P_(kPosX, kPosX) = kStdRadr * kStdRadr;
        P_(kPosY, kPosY) = kStdRadr * kStdRadr;
        P_(kVelocity, kVelocity) = kStdRadrd * kStdRadrd;
        break;
      }

      default:
        assert(false);  // Unknown sensor type
        break;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  const double delta_t = (meas_package.timestamp_ - time_us_) * 1e-6;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  switch (meas_package.sensor_type_) {
    case MeasurementPackage::LASER: {
      UpdateLidar(meas_package);
      break;
    }

    case MeasurementPackage::RADAR: {
      UpdateRadar(meas_package);
      break;
    }

    default:
      assert(false);  // Unknown sensor type
      break;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  assert(delta_t >= 0);

  // Bound covariance matrix coefficients to prevent NaN-s and Inf-s
  Eigen::JacobiSVD<MatrixXd> svd(P_, Eigen::ComputeFullU);
  MatrixXd diagP(P_.rows(), P_.cols());
  diagP.setZero();
  diagP.diagonal() = svd.singularValues().cwiseMax(1e-10).cwiseMin(1e10);
  P_ = (svd.matrixU() * diagP * svd.matrixU().transpose()).eval();

  MatrixXd P_aug(static_cast<int>(kNxAug), static_cast<int>(kNxAug));
  P_aug.setZero();
  P_aug.topLeftCorner(kNx, kNx) = P_;
  P_aug(kNoiseAccel, kNoiseAccel) = kStdA * kStdA;
  P_aug(kNoiseYawAccel, kNoiseYawAccel) = kStdYawdd * kStdYawdd;

  VectorXd x_aug(static_cast<int>(kNxAug));
  x_aug.setZero();
  x_aug.head<kNx>() = x_;

  const MatrixXd A = P_aug.llt().matrixL();
  Xsig_pred_.col(0) = x_aug;

  const double coeff = std::sqrt(getLambda(kNxAug) + kNxAug);
  for (int i = 0; i < kNxAug; ++i) {
    Xsig_pred_.col(i + 1) = x_aug + coeff * A.col(i);
    Xsig_pred_.col(i + 1 + kNxAug) = x_aug - coeff * A.col(i);
  }

  using std::sin;
  using std::cos;

  for (int i = 0, cols = Xsig_pred_.cols(); i < cols; ++i) {
    const VectorXd src = Xsig_pred_.col(i);
    MatrixXd::ColXpr dst = Xsig_pred_.col(i);

    if (std::abs(src[kYawRate]) > 1e-5) {
      const double v_yr = src[kVelocity] / src[kYawRate];

      dst[kPosX] +=
          v_yr * (sin(src[kYaw] + src[kYawRate] * delta_t) - sin(src[kYaw]));

      dst[kPosY] +=
          v_yr * (-cos(src[kYaw] + src[kYawRate] * delta_t) + cos(src[kYaw]));

    } else {
      dst[kPosX] += src[kVelocity] * cos(src[kYaw]) * delta_t;
      dst[kPosY] += src[kVelocity] * sin(src[kYaw]) * delta_t;
    }

    dst[kPosX] += 0.5 * delta_t * delta_t * cos(src[kYaw]) * src[kNoiseAccel];
    dst[kPosY] += 0.5 * delta_t * delta_t * sin(src[kYaw]) * src[kNoiseAccel];

    dst[kVelocity] += delta_t * src[kNoiseAccel];

    dst[kYaw] = normPi(
        src[kYaw] + delta_t * src[kYawRate]
            + 0.5 * delta_t * delta_t * src[kNoiseYawAccel]);

    dst[kYawRate] += delta_t * src[kNoiseYawAccel];
  }

  const double otherW = getW1ThroughN(kNxAug);
  const double w0Mean = getW0Mean(kNxAug);
  const double w0Cov = getW0Covariance(kNxAug);

  // 0th item is added in the end for slightly better numeric robustness.
  // It usually has the largest magnitude.

  x_.setZero();
  double sinYaw = 0.0;
  double cosYaw = 0.0;
  for (int i = 1, cols = Xsig_pred_.cols(); i < cols; ++i) {
    x_ += otherW * Xsig_pred_.col(i).head<kNx>();

    const double yaw = Xsig_pred_.col(i)(kYaw);
    sinYaw += otherW * sin(yaw);
    cosYaw += otherW * cos(yaw);
  }
  x_ += w0Mean * Xsig_pred_.col(0).head<kNx>();
  const double yaw = Xsig_pred_.col(0)(kYaw);
  sinYaw += w0Mean * sin(yaw);
  cosYaw += w0Mean * cos(yaw);

  if (std::abs(sinYaw) > 1e-10 || std::abs(cosYaw) > 1e-10) {
    x_[kYaw] = std::atan2(sinYaw, cosYaw);
  } else {
    x_[kYaw] = Xsig_pred_.col(0)(kYaw);
  }

  P_.setZero();

  for (int i = 1, cols = Xsig_pred_.cols(); i < cols; ++i) {
    VectorXd diff = Xsig_pred_.col(i).head<kNx>() - x_;
    diff[kYaw] = normPi(diff[kYaw]);
    P_ += otherW * diff * diff.transpose();
  }

  VectorXd diff = Xsig_pred_.col(0).head<kNx>() - x_;
  diff[kYaw] = normPi(diff[kYaw]);
  P_ += w0Cov * diff * diff.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& meas_package) {
  assert(MeasurementPackage::LASER == meas_package.sensor_type_);

  // Since LASER is a linear sensor, Kalman Filter produces same result as UKF

  MatrixXd R(2, 2);
  R.setIdentity();
  R(0, 0) = kStdLaspx * kStdLaspx;
  R(1, 1) = kStdLaspy * kStdLaspy;

  MatrixXd H(2, 5);
  H.setZero();
  H(0, 0) = H(1, 1) = 1.0;

  const MatrixXd S = (H * P_ * H.transpose() + R);
  const MatrixXd K = (P_ * H.transpose() * S.inverse());
  const VectorXd y = meas_package.raw_measurements_ - H * x_;

  x_ += K * y;
  x_(kYaw) = normPi(x_(kYaw));

  // More stable Joseph’s form covariance update from
  // https://en.wikipedia.org/wiki/Kalman_filter
  // instead of P_ -= (K * H * P_).eval();
  MatrixXd identity(P_.rows(), P_.cols());
  identity.setIdentity();

  P_ = ((identity - K * H) * P_ * (identity - K * H).transpose()
      + K * R * K.transpose()).eval();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage& meas_package) {
  assert(MeasurementPackage::RADAR == meas_package.sensor_type_);

  MatrixXd Zsig(3, static_cast<int>(2 * kNxAug + 1));

  using std::sin;
  using std::cos;

  for (int i = 0, cols = Xsig_pred_.cols(); i < cols; ++i) {
    const MatrixXd::ColXpr src = Xsig_pred_.col(i);
    MatrixXd::ColXpr dst = Zsig.col(i);

    dst[0] = std::sqrt(src[kPosX] * src[kPosX] + src[kPosY] * src[kPosY]);
    if(dst[0] < 1e-10) {
      return;
    }

    dst[1] = std::atan2(src[kPosY], src[kPosX]);

    dst[2] = (src[kPosX] * src[kVelocity] * std::cos(src[kYaw])
        + src[kPosY] * src[kVelocity] * std::sin(src[kYaw])) / dst[0];
  }

  const double otherW = getW1ThroughN(kNxAug);
  const double w0Mean = getW0Mean(kNxAug);
  const double w0Cov = getW0Covariance(kNxAug);

  double sinPhi = 0.0;
  double cosPhi = 0.0;

  VectorXd z_pred(3);
  z_pred.setZero();
  for (int i = 1, cols = Zsig.cols(); i < cols; ++i) {
    z_pred += otherW * Zsig.col(i);

    const double phi = Zsig.col(i)(1);
    sinPhi += otherW * sin(phi);
    cosPhi += otherW * cos(phi);
  }
  z_pred += w0Mean * Zsig.col(0);

  const double phi = Zsig.col(0)(1);
  sinPhi += w0Mean * sin(phi);
  cosPhi += w0Mean * cos(phi);
  if (std::abs(sinPhi) > 1e-10 || std::abs(cosPhi) > 1e-10) {
    z_pred[1] = std::atan2(sinPhi, cosPhi);
  } else {
    z_pred[1] = Zsig.col(0)(1);
  }

  MatrixXd S(3, 3);
  S.setIdentity();
  S(0, 0) = kStdRadr * kStdRadr;
  S(1, 1) = kStdRadphi * kStdRadphi;
  S(2, 2) = kStdRadrd * kStdRadrd;
  for (int i = 1, cols = Zsig.cols(); i < cols; ++i) {
    VectorXd diff = Zsig.col(i) - z_pred;
    diff[1] = normPi(diff[1]);
    S += otherW * diff * diff.transpose();
  }

  VectorXd diff = Zsig.col(0) - z_pred;
  diff[1] = normPi(diff[1]);
  S += w0Cov * diff * diff.transpose();

  MatrixXd Tc(static_cast<int>(kNx), 3);
  Tc.setZero();
  for(int i = 0, cols = Zsig.cols(); i < cols; ++i) {
    VectorXd diffX = Xsig_pred_.col(i).head<kNx>() - x_;
    diffX[kYaw] = normPi(diffX[kYaw]);
    VectorXd diffZ = Zsig.col(i) - z_pred;
    diff[1] = normPi(diff[1]);

    const double w = 0 == i ? w0Cov : otherW;
    Tc += w * diffX * diffZ.transpose();
  }

  const MatrixXd K = Tc * S.inverse();
  VectorXd diffZ = meas_package.raw_measurements_ - z_pred;
  diffZ[1] = normPi(diffZ[1]);
  x_ += K * diffZ;
  x_(kYaw) = normPi(x_(kYaw));
  P_ -= K * S * K.transpose();
}
