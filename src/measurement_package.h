#ifndef SRC_MEASUREMENT_PACKAGE_H_
#define SRC_MEASUREMENT_PACKAGE_H_

#include "Eigen/Dense"
#include <cstdint>

class MeasurementPackage {
 public:
  std::int64_t timestamp_;

  enum SensorType {
    LASER,
    RADAR
  } sensor_type_;

  Eigen::VectorXd raw_measurements_;
};

#endif  // SRC_MEASUREMENT_PACKAGE_H_
