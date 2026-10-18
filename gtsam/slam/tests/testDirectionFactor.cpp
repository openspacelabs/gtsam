/**
 * @file testDirectionFactor.cpp
 * @brief Unit tests for the building-direction factor
 * @author OpenSpace
 */

#include <CppUnitLite/TestHarness.h>
#include <gtsam/base/TestableAssertions.h>
#include <gtsam/base/numericalDerivative.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/PriorFactor.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/nonlinear/factorTesting.h>
#include <gtsam/slam/DirectionFactor.h>

#include <cmath>
#include <random>

using namespace gtsam;
using namespace gtsam::symbol_shorthand;

static const double kBuilding = 20.0 * M_PI / 180.0;  // the building's direction, world frame

/// Great-circle normals, in the body frame of a camera with rotation R, of horizontal edges
/// along the building's two axes seen from points around the camera, with a little noise.
static Matrix edgeNormals(const Rot3& R, size_t n, unsigned seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> u(-4.0, 4.0), z(0.5, 1.5), coin(0.0, 1.0);
  std::normal_distribution<double> noise(0.0, 0.3 * M_PI / 180.0);
  Matrix out(n, 3);
  for (size_t k = 0; k < n; ++k) {
    const double a = kBuilding + (coin(gen) < 0.5 ? M_PI / 2 : 0.0);
    const Vector3 d(std::cos(a), std::sin(a), 0.0);
    const Vector3 p(u(gen), u(gen), (coin(gen) < 0.5 ? -1 : 1) * z(gen));
    Vector3 nw = p.cross(d).normalized();
    nw += Vector3(noise(gen), noise(gen), noise(gen));
    out.row(k) = (R.inverse().matrix() * nw.normalized()).transpose();
  }
  return out;
}

static const SharedNoiseModel kModel = noiseModel::Isotropic::Sigma(1, 1.0 * M_PI / 180.0);

/* ************************************************************************* */
TEST(DirectionFactor, MeasuresTheBuildingAndTurnsWithTheFrame) {
  const Rot3 R = Rot3::RzRyRx(0.05, -0.03, 0.7);
  const Matrix normals = edgeNormals(R, 40, 1);
  DirectionFactor factor(X(0), A(0), normals, Vector::Constant(40, 0.1), kModel);
  // Under the true rotation the edges give the building's direction, to the noise
  EXPECT_DOUBLES_EQUAL(0.0, DirectionFactor::wrap90(factor.meanDirection(R) - kBuilding), 0.01);
  EXPECT_DOUBLES_EQUAL(0.0, factor.evaluateError(Pose3(R, Point3(1, 2, 3)), kBuilding)(0), 0.01);
  // A frame turned 5 degrees about the vertical measures the direction 5 degrees over
  const double turn = 5.0 * M_PI / 180.0;
  const Rot3 turned = Rot3::Rz(turn) * R;
  EXPECT_DOUBLES_EQUAL(turn, factor.evaluateError(Pose3(turned, Point3()), kBuilding)(0), 0.01);
  // ...and the error wraps modulo 90 degrees
  const Rot3 quarter = Rot3::Rz(M_PI / 2 + turn) * R;
  EXPECT_DOUBLES_EQUAL(turn, factor.evaluateError(Pose3(quarter, Point3()), kBuilding)(0), 0.01);
}

/* ************************************************************************* */
TEST(DirectionFactor, Jacobians) {
  const Rot3 R = Rot3::RzRyRx(0.1, -0.2, 0.7);
  const Matrix normals = edgeNormals(R, 30, 2);
  Vector weights(30);
  for (int k = 0; k < 30; ++k) weights(k) = 0.05 + 0.01 * k;
  DirectionFactor factor(X(0), A(0), normals, weights, kModel);
  for (const double yaw : {0.0, 0.08, -0.3}) {
    Values values;
    values.insert(X(0), Pose3(Rot3::Rz(yaw) * R, Point3(0.5, -1.0, 2.0)));
    values.insert(A(0), kBuilding + 0.02);
    EXPECT_CORRECT_FACTOR_JACOBIANS(factor, values, 1e-7, 1e-6);
  }
  // The mean direction's own derivative against the rotation, against finite differences
  Matrix13 H;
  factor.meanDirection(R, H);
  const Matrix expected = numericalDerivative11<double, Rot3>(
      [&factor](const Rot3& r) { return factor.meanDirection(r); }, R, 1e-7);
  EXPECT(assert_equal(expected, Matrix(H), 1e-6));
}

/* ************************************************************************* */
TEST(DirectionFactor, TurnsAFrameBackAndSetsTheDirection) {
  // Two keyframes with the heading 5 and -3 degrees off, the direction free: the solve turns
  // both to the direction they share, and the direction lands on the building's.
  const Rot3 R0 = Rot3::RzRyRx(0.02, 0.01, 0.4), R1 = Rot3::RzRyRx(-0.01, 0.03, 2.1);
  const Matrix n0 = edgeNormals(R0, 50, 3), n1 = edgeNormals(R1, 50, 4);
  NonlinearFactorGraph graph;
  graph.emplace_shared<DirectionFactor>(X(0), A(0), n0, Vector::Constant(50, 0.1), kModel);
  graph.emplace_shared<DirectionFactor>(X(1), A(0), n1, Vector::Constant(50, 0.1), kModel);
  // Roll, pitch and position held (gravity and the odometry would); heading free
  const Vector6 sig = (Vector6() << 1e-4, 1e-4, 10.0, 1e-4, 1e-4, 1e-4).finished();
  const Pose3 start0(Rot3::Rz(5.0 * M_PI / 180) * R0, Point3(0, 0, 0));
  const Pose3 start1(Rot3::Rz(-3.0 * M_PI / 180) * R1, Point3(1, 0, 0));
  graph.addPrior(X(0), start0, noiseModel::Diagonal::Sigmas(sig));
  graph.addPrior(X(1), start1, noiseModel::Diagonal::Sigmas(sig));
  // The direction has a weak prior on the building's value: the two frames alone would fix
  // only their difference
  graph.addPrior(A(0), kBuilding, noiseModel::Isotropic::Sigma(1, 0.5 * M_PI / 180));
  Values initial;
  initial.insert(X(0), start0);
  initial.insert(X(1), start1);
  initial.insert(A(0), kBuilding);
  const Values result = LevenbergMarquardtOptimizer(graph, initial).optimize();
  EXPECT(assert_equal(R0, result.at<Pose3>(X(0)).rotation(), 0.01));
  EXPECT(assert_equal(R1, result.at<Pose3>(X(1)).rotation(), 0.01));
  EXPECT_DOUBLES_EQUAL(kBuilding, result.at<double>(A(0)), 0.01);
}

/* ************************************************************************* */
TEST(DirectionFactor, EqualsAndClone) {
  const Rot3 R = Rot3::RzRyRx(0.1, -0.2, 0.7);
  DirectionFactor factor(X(0), A(0), edgeNormals(R, 10, 5), Vector::Constant(10, 0.1), kModel);
  DirectionFactor other(X(0), A(0), edgeNormals(R, 10, 6), Vector::Constant(10, 0.1), kModel);
  EXPECT(factor.equals(*factor.clone()));
  EXPECT(!factor.equals(other));
  EXPECT_LONGS_EQUAL(10, factor.normals().rows());
}

/* ************************************************************************* */
int main() {
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */
