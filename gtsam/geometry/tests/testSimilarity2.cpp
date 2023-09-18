/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file   testSimilarity2.cpp
 * @brief  Unit tests for Similarity2 class
 * @author Varun Agrawal
 */

#include <CppUnitLite/TestHarness.h>
#include <gtsam/base/Testable.h>
#include <gtsam/base/numericalDerivative.h>
#include <gtsam/base/testLie.h>
#include <gtsam/geometry/Similarity2.h>

#include <functional>

using namespace std::placeholders;
using namespace gtsam;
using namespace std;

GTSAM_CONCEPT_TESTABLE_INST(Similarity2)

static const Point2 P(0.2, 0.7);
static const Rot2 R = Rot2::fromAngle(0.3);
static const double s = 4;

//******************************************************************************
TEST(Similarity2, Concepts) {
  BOOST_CONCEPT_ASSERT((IsGroup<Similarity2>));
  BOOST_CONCEPT_ASSERT((IsManifold<Similarity2>));
  BOOST_CONCEPT_ASSERT((IsLieGroup<Similarity2>));
}

//******************************************************************************
TEST(Similarity2, Constructors) {
  Similarity2 sim2_Construct1;
  Similarity2 sim2_Construct2(s);
  Similarity2 sim2_Construct3(R, P, s);
  Similarity2 sim2_Construct4(R.matrix(), P, s);
}

//******************************************************************************
TEST(Similarity2, Getters) {
  Similarity2 sim2_default;
  EXPECT(assert_equal(Rot2(), sim2_default.rotation()));
  EXPECT(assert_equal(Point2(0, 0), sim2_default.translation()));
  EXPECT_DOUBLES_EQUAL(1.0, sim2_default.scale(), 1e-9);
}

//******************************************************************************
TEST(Similarity2, AdjointMap) {
  Similarity2 T(R, P, s);
  Vector4 xi(1, 2, 3e-4, 5e-6);
  Similarity2 expected = T * Similarity2::Expmap(xi) * T.inverse();
  Vector xiprime = T.AdjointMap() * xi;
  EXPECT(assert_equal(expected, Similarity2::Expmap(xiprime), 1e-6));

  Vector4 xi2(4, -3, 2e-4, -1e-4);
  Similarity2 expected2 = T * Similarity2::Expmap(xi2) * T.inverse();
  Vector xiprime2 = T.AdjointMap() * xi2;
  EXPECT(assert_equal(expected2, Similarity2::Expmap(xiprime2), 1e-6));
}

//******************************************************************************
int main() {
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
//******************************************************************************
