#include <gtest/gtest.h>
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
  ASSERT_EQ(measurement.timestamp_, ukf.time_us_);
  ASSERT_PRED2(IsEigenEqual(), measurement.raw_measurements_, ukf.x_.head(2));
}
