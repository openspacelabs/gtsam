/**
 * @file testInverseDepthBearingFactor.cpp
 * @brief Unit tests for the anchored inverse-depth bearing factor
 * @author OpenSpace
 */

#include <CppUnitLite/TestHarness.h>
#include <gtsam/base/TestableAssertions.h>
#include <gtsam/base/numericalDerivative.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/nonlinear/factorTesting.h>
#include <gtsam/slam/InverseDepthBearingFactor.h>

#include <cmath>

using namespace gtsam;
using namespace gtsam::symbol_shorthand;

static const SharedNoiseModel kModel = noiseModel::Isotropic::Sigma(2, 0.01);
static const Pose3 kAnchor(Rot3::RzRyRx(0.1, -0.2, 0.3), Point3(1.0, 2.0, 0.5));
static const Pose3 kObserver(Rot3::RzRyRx(-0.3, 0.1, 0.2), Point3(2.5, 1.0, 0.7));
static const Unit3 kAnchorBearing(0.3, -0.5, 0.8);

/* ************************************************************************* */
TEST(InverseDepthBearingFactor, AgreesWithPointGeometry) {
  // A landmark 4 m along the anchor bearing, seen from the observer along the bearing GTSAM's own
  // point geometry gives: the factor predicts that bearing and has zero error for it.
  const double rho = 0.25;
  const Point3 landmark = InverseDepthBearingFactor::Landmark(kAnchor, kAnchorBearing, rho);
  const Point3 expected =
      kAnchor.translation() + 4.0 * kAnchor.rotation().rotate(kAnchorBearing.point3());
  EXPECT(assert_equal(expected, landmark, 1e-9));
  const Unit3 measured = kObserver.bearing(landmark);
  InverseDepthBearingFactor factor(X(0), X(1), L(0), kAnchorBearing, measured, kModel);
  EXPECT(assert_equal(measured, factor.predict(kAnchor, kObserver, rho), 1e-9));
  EXPECT(assert_equal(Vector2::Zero(), factor.evaluateError(kAnchor, kObserver, rho), 1e-9));
  // The anchor itself sees the landmark along the anchor bearing at any depth.
  InverseDepthBearingFactor self(X(0), X(0), L(0), kAnchorBearing, kAnchorBearing, kModel);
  EXPECT(assert_equal(Vector2::Zero(), self.evaluateError(kAnchor, kAnchor, rho), 1e-9));
  EXPECT(assert_equal(Vector2::Zero(), self.evaluateError(kAnchor, kAnchor, 0.0), 1e-9));
}

/* ************************************************************************* */
TEST(InverseDepthBearingFactor, ErrorIsSineOfAngle) {
  // A measurement 2 degrees off the prediction: the error's norm is sin(2 deg), whatever the
  // tangent basis.
  const double rho = 0.25;
  InverseDepthBearingFactor probe(X(0), X(1), L(0), kAnchorBearing, Unit3(1, 0, 0), kModel);
  const Unit3 predicted = probe.predict(kAnchor, kObserver, rho);
  const double angle = 2.0 * M_PI / 180.0;
  const Unit3 measured = Rot3::AxisAngle(Point3(predicted.basis().col(0)), angle).rotate(predicted);
  InverseDepthBearingFactor factor(X(0), X(1), L(0), kAnchorBearing, measured, kModel);
  const Vector2 e = factor.evaluateError(kAnchor, kObserver, rho);
  EXPECT_DOUBLES_EQUAL(std::sin(angle), e.norm(), 1e-9);
}

/* ************************************************************************* */
TEST(InverseDepthBearingFactor, Jacobians) {
  const double rho = 0.25;
  const Point3 landmark = InverseDepthBearingFactor::Landmark(kAnchor, kAnchorBearing, rho);
  // A measurement a little off the truth, so the error is not zero where the derivatives are taken.
  const Unit3 measured = kObserver.bearing(landmark + Point3(0.05, -0.03, 0.02));
  InverseDepthBearingFactor factor(X(0), X(1), L(0), kAnchorBearing, measured, kModel);
  Values values;
  values.insert(X(0), kAnchor);
  values.insert(X(1), kObserver);
  values.insert(L(0), rho);
  EXPECT_CORRECT_FACTOR_JACOBIANS(factor, values, 1e-6, 1e-6);
}

/* ************************************************************************* */
TEST(InverseDepthBearingFactor, JacobiansAtInfinity) {
  // At zero inverse depth the bearing depends on the rotations only: the derivatives against both
  // positions vanish, those against the rotations do not, and all of them match finite differences.
  InverseDepthBearingFactor probe(X(0), X(1), L(0), kAnchorBearing, Unit3(1, 0, 0), kModel);
  const Unit3 predicted = probe.predict(kAnchor, kObserver, 0.0);
  const Unit3 measured = Rot3::AxisAngle(Point3(predicted.basis().col(1)), 0.01).rotate(predicted);
  InverseDepthBearingFactor factor(X(0), X(1), L(0), kAnchorBearing, measured, kModel);
  Values values;
  values.insert(X(0), kAnchor);
  values.insert(X(1), kObserver);
  values.insert(L(0), 0.0);
  EXPECT_CORRECT_FACTOR_JACOBIANS(factor, values, 1e-6, 1e-6);
  Matrix H1, H2, H3;
  factor.evaluateError(kAnchor, kObserver, 0.0, H1, H2, H3);
  EXPECT(assert_equal(Matrix::Zero(2, 3), Matrix(H1.rightCols(3)), 1e-12));
  EXPECT(assert_equal(Matrix::Zero(2, 3), Matrix(H2.rightCols(3)), 1e-12));
  EXPECT(H1.leftCols(3).norm() > 0.1);
  EXPECT(H2.leftCols(3).norm() > 0.1);
  // The anchor's own observation never moves, so its factor has no derivative at all.
  InverseDepthBearingFactor self(X(0), X(0), L(0), kAnchorBearing, kAnchorBearing, kModel);
  self.evaluateError(kAnchor, kAnchor, 0.0, H1, H2, H3);
  EXPECT(assert_equal(Matrix::Zero(2, 1), H3, 1e-12));
}

/* ************************************************************************* */
TEST(InverseDepthBearingFactor, RecoversDepthsAndPose) {
  // Three poses, four landmarks anchored at the first at known depths, exact bearings. Start the
  // inverse depths at half their value and the third pose perturbed: the optimizer recovers both.
  const std::vector<Pose3> poses = {
      kAnchor, kObserver, Pose3(Rot3::RzRyRx(0.2, 0.05, -0.1), Point3(0.5, 3.0, 1.2))};
  const std::vector<Unit3> bearings = {Unit3(0.3, -0.5, 0.8), Unit3(-0.7, 0.2, 0.1),
                                       Unit3(0.1, 0.9, -0.4), Unit3(0.5, 0.5, 0.5)};
  const std::vector<double> depths = {3.0, 6.0, 2.0, 10.0};
  NonlinearFactorGraph graph;
  Values initial;
  graph.addPrior(X(0), poses[0], noiseModel::Isotropic::Sigma(6, 1e-3));
  graph.addPrior(X(1), poses[1], noiseModel::Isotropic::Sigma(6, 1e-3));
  KeyVector inverseDepths;
  for (size_t k = 0; k < bearings.size(); ++k) {
    const Point3 landmark =
        InverseDepthBearingFactor::Landmark(poses[0], bearings[k], 1.0 / depths[k]);
    for (size_t j = 1; j < poses.size(); ++j)
      graph.emplace_shared<InverseDepthBearingFactor>(X(0), X(j), L(k), bearings[k],
                                                      poses[j].bearing(landmark), kModel);
    initial.insert(L(k), 0.5 / depths[k]);
    inverseDepths.push_back(L(k));
  }
  initial.insert(X(0), poses[0]);
  initial.insert(X(1), poses[1]);
  initial.insert(X(2), poses[2].retract((Vector(6) << 0.02, -0.01, 0.03, 0.2, -0.1, 0.15).finished()));
  // Inverse depths carry little information (a short baseline at one bearing), so they are
  // eliminated first: the Cholesky elimination declares a clique underconstrained when the last
  // two pivots of its frontal block differ by more than 2^12, as a scalar with that little
  // information does next to a pose held by a prior (see choleskyPartial).
  LevenbergMarquardtParams params;
  params.setOrdering(Ordering::ColamdConstrainedFirst(graph, inverseDepths));
  const Values result = LevenbergMarquardtOptimizer(graph, initial, params).optimize();
  for (size_t k = 0; k < bearings.size(); ++k)
    EXPECT_DOUBLES_EQUAL(1.0 / depths[k], result.at<double>(L(k)), 1e-6);
  EXPECT(assert_equal(poses[2], result.at<Pose3>(X(2)), 1e-5));
  EXPECT_DOUBLES_EQUAL(0.0, graph.error(result), 1e-9);
}

/* ************************************************************************* */
TEST(InverseDepthBearingFactor, SensorPoses) {
  // Bearings from lenses displaced from the body: the factor agrees with the point geometry of
  // the composed sensor poses, and its derivatives against the body poses match finite
  // differences.
  const Pose3 anchorSensor(Rot3::RzRyRx(0.0, 0.0, M_PI / 2), Point3(0.03, 0.0, 0.0));
  const Pose3 observerSensor(Rot3::RzRyRx(0.0, 0.0, -M_PI / 2), Point3(-0.03, 0.0, 0.01));
  const double rho = 0.25;
  const Point3 landmark =
      InverseDepthBearingFactor::Landmark(kAnchor.compose(anchorSensor), kAnchorBearing, rho);
  const Unit3 measured = kObserver.compose(observerSensor).bearing(landmark);
  InverseDepthBearingFactor factor(X(0), X(1), L(0), kAnchorBearing, measured, kModel,
                                   anchorSensor, observerSensor);
  EXPECT(assert_equal(landmark, factor.landmark(kAnchor, rho), 1e-9));
  EXPECT(assert_equal(measured, factor.predict(kAnchor, kObserver, rho), 1e-9));
  EXPECT(assert_equal(Vector2::Zero(), factor.evaluateError(kAnchor, kObserver, rho), 1e-9));
  // Without the sensor poses the same measurement is off
  InverseDepthBearingFactor bare(X(0), X(1), L(0), kAnchorBearing, measured, kModel);
  EXPECT(bare.evaluateError(kAnchor, kObserver, rho).norm() > 1e-3);
  EXPECT(!factor.equals(bare));
  EXPECT(factor.equals(*factor.clone()));
  // Derivatives, off the truth and at infinity
  for (const double r : {rho, 0.0}) {
    const Unit3 off = Rot3::AxisAngle(Point3(measured.basis().col(0)), 0.01).rotate(measured);
    InverseDepthBearingFactor f(X(0), X(1), L(0), kAnchorBearing, off, kModel, anchorSensor,
                                observerSensor);
    Values values;
    values.insert(X(0), kAnchor);
    values.insert(X(1), kObserver);
    values.insert(L(0), r);
    EXPECT_CORRECT_FACTOR_JACOBIANS(f, values, 1e-6, 1e-6);
  }
  // Only one side displaced
  InverseDepthBearingFactor half(X(0), X(1), L(0), kAnchorBearing, measured, kModel,
                                 anchorSensor);
  Values values;
  values.insert(X(0), kAnchor);
  values.insert(X(1), kObserver);
  values.insert(L(0), rho);
  EXPECT_CORRECT_FACTOR_JACOBIANS(half, values, 1e-6, 1e-6);
}

/* ************************************************************************* */
TEST(InverseDepthBearingFactor, EqualsAndClone) {
  InverseDepthBearingFactor factor(X(0), X(1), L(0), kAnchorBearing, Unit3(0.1, 0.2, 0.9), kModel);
  InverseDepthBearingFactor other(X(0), X(1), L(0), kAnchorBearing, Unit3(0.1, 0.3, 0.9), kModel);
  EXPECT(factor.equals(*factor.clone()));
  EXPECT(!factor.equals(other));
  EXPECT(assert_equal(kAnchorBearing, factor.anchorBearing()));
  EXPECT(assert_equal(Unit3(0.1, 0.2, 0.9), factor.measured()));
}

/* ************************************************************************* */
int main() {
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */
