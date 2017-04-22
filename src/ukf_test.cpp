#include <gtest/gtest.h>
#include <cmath>
#include "test_utils.h"
#include "ukf.h"

TEST(UkfProcessMeasurement, FirstLaser) {
  UKF ukf;
  MeasurementPackage measurement;
  measurement.timestamp_ = 1;
  measurement.sensor_type_ = measurement.LASER;
  measurement.raw_measurements_.resize(2);
  measurement.raw_measurements_ << 1.0, -2.0;
  ukf.ProcessMeasurement(measurement);

  ASSERT_TRUE(ukf.is_initialized_);
  ASSERT_EQ(1, ukf.time_us_);
  ASSERT_PRED2(IsEigenEqual(), measurement.raw_measurements_, ukf.x_.head(2));
}

TEST(UkfProcessMeasurement, FirstRadar) {
  UKF ukf;
  MeasurementPackage measurement;
  measurement.timestamp_ = 2;
  measurement.sensor_type_ = measurement.RADAR;
  measurement.raw_measurements_.resize(3);
  measurement.raw_measurements_ << 1.0, M_PI / 4, 3.0;
  ukf.ProcessMeasurement(measurement);

  ASSERT_TRUE(ukf.is_initialized_);
  ASSERT_EQ(2, ukf.time_us_);
  ASSERT_NEAR(1.0 / std::sqrt(2), ukf.x_(UKF::kPosX), 1e-10);
  ASSERT_NEAR(1.0 / std::sqrt(2), ukf.x_(UKF::kPosY), 1e-10);
  ASSERT_NEAR(3.0, ukf.x_(UKF::kVelocity), 1e-10);
  ASSERT_NEAR(M_PI / 4, ukf.x_(UKF::kYaw), 1e-10);
  ASSERT_NEAR(0.0, ukf.x_(UKF::kYawRate), 1e-10);
}

TEST(UkfPrediction, NoMotion) {
  UKF ukf;
  ukf.x_ << 10.0, 20.0, 0.0, 0.1, 0.0;
  ukf.P_.setIdentity();
  ukf.P_ *= 100.0;

  ukf.Prediction(1.0);

  Eigen::VectorXd expected(5);
  expected << 10.0, 20.0, 0.0, 0.1, 0.0;

  ASSERT_PRED2(IsEigenEqual(), expected, ukf.x_);
}

TEST(UkfPrediction, StraightLine) {
  UKF ukf;
  ukf.x_ << 1.0, 2.0, 3.0, 0.0, 0.0;
  ukf.P_.setIdentity();
  ukf.P_ *= 0.0001;

  ukf.Prediction(1.0);

  Eigen::VectorXd expected(5);
  expected << 4.0, 2.0, 3.0, 0.0, 0.0;

  ASSERT_PRED2(IsEigenEqual(), expected, ukf.x_);
  ASSERT_PRED2(IsEigenEqual(), ukf.P_, ukf.P_.transpose());
  ASSERT_GT(ukf.P_.eigenvalues()[0].real(), 0.0001);
}

TEST(UkfPrediction, FullCircleInOneUpdate) {
  UKF ukf;
  ukf.x_ << 1.0, 2.0, 2 * M_PI, 0.0, 2 * M_PI;
  ukf.P_.setIdentity();
  ukf.P_ *= 0.0001;

  ukf.Prediction(1.0);

  Eigen::VectorXd expected(5);
  expected << 1.0, 2.0, 2 * M_PI, 0.0, 2 * M_PI;

  ASSERT_PRED2(IsEigenEqual(), expected, ukf.x_);
  ASSERT_PRED2(IsEigenEqual(), ukf.P_, ukf.P_.transpose());
  ASSERT_GT(ukf.P_.eigenvalues()[0].real(), 0.0001);
}

TEST(UkfPrediction, QuaterCircleInOneUpdate) {
  UKF ukf;
  ukf.x_ << 1.0, 0.0, 2 * M_PI, M_PI / 2, 2 * M_PI;
  ukf.P_.setIdentity();
  ukf.P_ *= 0.0001;

  ukf.Prediction(0.25);

  Eigen::VectorXd expected(5);
  expected << 0.0, 1.0, 2 * M_PI, M_PI, 2 * M_PI;

  ASSERT_PRED2(IsEigenEqual(), expected, ukf.x_);
  ASSERT_PRED2(IsEigenEqual(), ukf.P_, ukf.P_.transpose());
  ASSERT_GT(ukf.P_.eigenvalues()[0].real(), 0.0001);
}

TEST(UkfProcessMeasurement, LaserTrajectory) {
  UKF ukf;
  MeasurementPackage measurement;
  measurement.sensor_type_ = MeasurementPackage::LASER;
  measurement.raw_measurements_.resize(2);

  measurement.timestamp_ = 0;

  const double kRadius = 30;

  for (int degrees = 0; degrees <= 360; ++degrees) {
    if (degrees < 180) {
      const double radians = degrees * M_PI / 180.0;

      measurement.raw_measurements_
          << kRadius * std::cos(radians), kRadius * std::sin(radians);

      ukf.ProcessMeasurement(measurement);

      measurement.timestamp_ += 1000000;

    } else {
      ukf.Prediction(1.0);
    }
  }

  Eigen::VectorXd expected(2);
  expected << kRadius, 0.0;
  ASSERT_PRED2(IsEigenEqual(), expected, ukf.x_.head(2));
}
