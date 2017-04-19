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
  ASSERT_NEAR(M_PI / 4, ukf.x_(UKF::kYawAngle), 1e-10);
  ASSERT_NEAR(0.0, ukf.x_(UKF::kYawRate), 1e-10);
}
