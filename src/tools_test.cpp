#include <gtest/gtest.h>
#include "test_utils.h"
#include "tools.h"

TEST(Tools, EvaluateRmse) {
  using Eigen::Vector3d;

  const Vector3d rmse = EvaluateRmse({
                                         Vector3d{0.5, 0.2, 0.0},
                                         Vector3d{1.0, 0.3, -5.0}
                                     }, {
                                         Vector3d{-0.5, -1.8, 3.0},
                                         Vector3d{2.0, 2.3, -2.0}
                                     });

  ASSERT_PRED2(IsEigenEqual(), (Vector3d{1.0, 2.0, 3.0}), rmse);
}
