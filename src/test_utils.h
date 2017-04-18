#ifndef SRC_TEST_UTILS_H_
#define SRC_TEST_UTILS_H_

#include "Eigen/Dense"

class IsEigenEqual {
 public:
  explicit IsEigenEqual(double eps = 1e-4)
      : eps_{eps} {}


  IsEigenEqual& operator=(const IsEigenEqual&) = delete;


  template<typename Derived1, typename Derived2>
  bool operator()(const Eigen::MatrixBase<Derived1>& m1,
                  const Eigen::MatrixBase<Derived2>& m2) {
    return m1.isApprox(m2, eps_);
  }

 private:
  const double eps_;
};

#endif  // SRC_TEST_UTILS_H_
