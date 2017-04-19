#include "ukf.h"
#include <iostream>
#include <cassert>
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
}  // namespace

UKF::UKF()
  : x_(static_cast<int>(UKF::kNx))
  , P_(static_cast<int>(UKF::kNx), static_cast<int>(UKF::kNx))
  , Xsig_pred_(static_cast<int>(UKF::kNx),
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
        x_(kYawAngle) = phi;
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
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
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
