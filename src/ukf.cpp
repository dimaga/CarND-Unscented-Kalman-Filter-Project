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

const double kAlpha{1.0};

const double kBeta{2.0};

double getKappa(int n) {
  return 3 - n;
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
  return 0.5 / (n + lambda);
}

/**
 * Angle normalization
 */
double normPi(double angleRad) {
  while (angleRad > M_PI) {
    angleRad -= 2 * M_PI;
  }

  while (angleRad < -M_PI) {
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
  P_ *= 10000;
  Xsig_pred_.setZero();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
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

  const double delta_t = (time_us_ - meas_package.timestamp_) * 1e-6;
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

  for (int i = 0, cols = Xsig_pred_.cols(); i < cols; ++i) {
    const Eigen::VectorXd src = Xsig_pred_.col(i);
    Eigen::MatrixXd::ColXpr dst = Xsig_pred_.col(i);

    using std::sin;
    using std::cos;

    if (std::abs(src[kYawRate]) > 1e-5) {
      dst[kPosX] = src[kPosX]
          + src[kVelocity] / src[kYawRate]
              * (sin(src[kYaw] + src[kYawRate] * delta_t) - sin(src[kYaw]))
          + 0.5 * delta_t * delta_t * cos(src[kYaw]) * src[kNoiseAccel];

      dst[kPosY] = src[kPosY]
          + src[kVelocity] / src[kYawRate]
              * (-cos(src[kYaw] + src[kYawRate] * delta_t) + cos(src[kYaw]))
          + 0.5 * delta_t * delta_t * sin(src[kYaw]) * src[kNoiseAccel];
    } else {
      dst[kPosX] = src[kPosX]
          + (src[kVelocity] * cos(src[kYaw])
          + 0.5 * delta_t * cos(src[kYaw]) * src[kNoiseAccel]) * delta_t;

      dst[kPosY] = src[kPosY]
          + (src[kVelocity] * sin(src[kYaw])
          + 0.5 * delta_t * sin(src[kYaw]) * src[kNoiseAccel]) * delta_t;
    }

    dst[kVelocity] = src[kVelocity] + delta_t * src[kNoiseAccel];

    dst[kYaw] = normPi(
        src[kYaw] + delta_t * src[kYawRate]
            + 0.5 * delta_t * delta_t * src[kNoiseYawAccel]);

    dst[kYawRate] = src[kYawRate] + delta_t * src[kNoiseYawAccel];
  }

  const double otherW = getW1ThroughN(kNxAug);
  const double w0Mean = getW0Mean(kNxAug);

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
  }

  P_.setZero();

  for (int i = 1, cols = Xsig_pred_.cols(); i < cols; ++i) {
    VectorXd diff = Xsig_pred_.col(i).head<kNx>() - x_;
    diff[kYaw] = normPi(diff[kYaw]);
    P_ += otherW * diff * diff.transpose();
  }

  VectorXd diff = Xsig_pred_.col(0).head<kNx>() - x_;
  diff[kYaw] = normPi(diff[kYaw]);
  P_ += getW0Covariance(kNxAug) * diff * diff.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  assert(MeasurementPackage::LASER == meas_package.sensor_type_);
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  assert(MeasurementPackage::RADAR == meas_package.sensor_type_);
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
