#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <stdlib.h>
#include "Eigen/Dense"
#include "ukf.h"
#include "tools.h"
#include "ground_truth_package.h"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::cout;
using std::cerr;

void check_arguments(int argc, char *argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (1 == argc) {
    cerr << usage_instructions << endl;
  } else if (2 == argc) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (3 == argc) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void check_files(const ifstream &in_file,
                 const string &in_name,
                 const ofstream &out_file,
                 const string &out_name) {
  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }

  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char *argv[]) {
  check_arguments(argc, argv);

  const string in_file_name_ = argv[1];
  ifstream in_file_(in_file_name_.c_str(), ifstream::in);

  const string out_file_name_ = argv[2];
  ofstream out_file_(out_file_name_.c_str(), ofstream::out);

  check_files(in_file_, in_file_name_, out_file_, out_file_name_);

  vector<MeasurementPackage> measurement_pack_list;
  vector<GroundTruthPackage> gt_pack_list;

  string line;

  // prep the measurement packages (each line represents a measurement at a
  // timestamp)
  while (getline(in_file_, line)) {
    string sensor_type;
    MeasurementPackage meas_package;
    GroundTruthPackage gt_package;
    std::istringstream iss(line);
    std::int64_t timestamp;

    // reads first element from the current line
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0) {
      // LASER MEASUREMENT

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      double x;
      double y;
      iss >> x;
      iss >> y;
      meas_package.raw_measurements_ << x, y;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    } else if (sensor_type.compare("R") == 0) {
      // RADAR MEASUREMENT

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      double ro;
      double phi;
      double ro_dot;
      iss >> ro;
      iss >> phi;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, phi, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }

    // read ground truth data to compare later
    double x_gt;
    double y_gt;
    double vx_gt;
    double vy_gt;
    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;
    gt_package.gt_values_ = VectorXd(4);
    gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
    gt_pack_list.push_back(gt_package);
  }

  // Create a UKF instance
  UKF ukf;

  // used to compute the RMSE later
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  // start filtering from the second frame (the speed is unknown in the first
  // frame)

  size_t number_of_measurements = measurement_pack_list.size();

  // column names for output file
  out_file_ << "px" << "\t";
  out_file_ << "py" << "\t";
  out_file_ << "v" << "\t";
  out_file_ << "yaw_angle" << "\t";
  out_file_ << "yaw_rate" << "\t";
  out_file_ << "px_measured" << "\t";
  out_file_ << "py_measured" << "\t";
  out_file_ << "px_true" << "\t";
  out_file_ << "py_true" << "\t";
  out_file_ << "vx_true" << "\t";
  out_file_ << "vy_true" << "\t";
  out_file_ << "NIS" << "\n";

  for (size_t k = 0; k < number_of_measurements; ++k) {
    const auto &measurement = measurement_pack_list[k];

    // Call the UKF-based fusion
    ukf.ProcessMeasurement(measurement_pack_list[k]);

    // output the estimation
    out_file_ << ukf.x_(UKF::kPosX) << "\t";
    out_file_ << ukf.x_(UKF::kPosY) << "\t";
    out_file_ << ukf.x_(UKF::kVelocity) << "\t";
    out_file_ << ukf.x_(UKF::kYaw) << "\t";
    out_file_ << ukf.x_(UKF::kYawRate) << "\t";

    using std::cos;
    using std::sin;

    // output the measurements
    if (MeasurementPackage::LASER == measurement.sensor_type_) {
      // output the estimation
      out_file_ << measurement.raw_measurements_(0) << "\t";
      out_file_ << measurement.raw_measurements_(1) << "\t";
    } else if (MeasurementPackage::RADAR == measurement.sensor_type_) {
      // output the estimation in the cartesian coordinates
      double ro = measurement.raw_measurements_(0);
      double phi = measurement.raw_measurements_(1);
      out_file_ << ro * cos(phi) << "\t";  // p1_meas
      out_file_ << ro * sin(phi) << "\t";  // ps_meas
    }

    // output the ground truth packages
    const auto &ground_truth_item = gt_pack_list[k].gt_values_;

    out_file_ << ground_truth_item(0) << "\t";
    out_file_ << ground_truth_item(1) << "\t";
    out_file_ << ground_truth_item(2) << "\t";
    out_file_ << ground_truth_item(3) << "\n";

    // output the NIS values

    if (MeasurementPackage::LASER == measurement.sensor_type_) {
      out_file_ << ukf.NIS_laser_ << "\n";
    } else if (MeasurementPackage::RADAR == measurement.sensor_type_) {
      out_file_ << ukf.NIS_radar_ << "\n";
    }


    // convert ukf x vector to cartesian to compare to ground truth
    VectorXd ukf_x_cartesian_ = VectorXd(4);

    const double x_estimate_ = ukf.x_(UKF::kPosX);
    const double y_estimate_ = ukf.x_(UKF::kPosY);
    const double vx_estimate_ = ukf.x_(UKF::kVelocity) * cos(ukf.x_(UKF::kYaw));
    const double vy_estimate_ = ukf.x_(UKF::kVelocity) * sin(ukf.x_(UKF::kYaw));

    ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;

    estimations.push_back(ukf_x_cartesian_);
    ground_truth.push_back(ground_truth_item);
  }

  // compute the accuracy (RMSE)
  cout << "Accuracy - RMSE:" << endl;
  cout << EvaluateRmse(estimations, ground_truth) << endl;

  // close files
  if (out_file_.is_open()) {
    out_file_.close();
  }

  if (in_file_.is_open()) {
    in_file_.close();
  }

  cout << "Done!" << endl;
  return 0;
}
