#ifndef SRC_UKF_H_
#define SRC_UKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <cstdint>
#include <cstddef>

class UKF {
  using MatrixXd = Eigen::MatrixXd;
  using VectorXd = Eigen::VectorXd;

 public:
  enum {
    kPosX,
    kPosY,
    kVelocity,
    kYaw,
    kYawRate,

    kNx,

    kNoiseAccel = kNx,
    kNoiseYawAccel,
    kNxAug
  };

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_{false};

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_{true};

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_{true};

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  std::int64_t time_us_{0};

  ///* the current NIS for radar
  double NIS_radar_{0};

  ///* the current NIS for laser
  double NIS_laser_{0};


  UKF();
  virtual ~UKF() = default;

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(const MeasurementPackage& meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const MeasurementPackage& meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const MeasurementPackage& meas_package);
};

#endif  // SRC_UKF_H_
