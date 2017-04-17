#include <iostream>
#include <vector>
#include <cmath>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;


VectorXd EvaluateRmse(const vector<VectorXd> &estimations,
                      const vector<VectorXd> &ground_truth) {
  assert(!estimations.empty());
  assert(estimations.size() == ground_truth.size());

  VectorXd resSq(estimations.front().size());
  resSq.setZero();

  for (std::size_t i = 0; i < estimations.size(); ++i) {
    const auto& est = estimations.at(i);
    const auto& gt = ground_truth.at(i);
    const auto residual = (est - gt).array();
    resSq.array() += residual * residual;
  }

  resSq /= static_cast<double>(estimations.size());

  return resSq.cwiseSqrt();
}
