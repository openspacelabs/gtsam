"""
GTSAM Copyright 2010-2019, Georgia Tech Research Corporation,
Atlanta, Georgia 30332-0415
All Rights Reserved

See LICENSE for the license information

InverseDepthBearingFactor unit tests.
Author: OpenSpace
"""

import unittest

import gtsam
import numpy as np
from gtsam import InverseDepthBearingFactor, Point3, Pose3, Rot3, Unit3
from gtsam.symbol_shorthand import L, X
from gtsam.utils.test_case import GtsamTestCase


def bearing(pose: Pose3, point: Point3) -> Unit3:
    """The bearing from a pose to a world point, in the pose's frame."""
    return Unit3(pose.transformTo(point))


class TestInverseDepthBearingFactor(GtsamTestCase):
    """A landmark anchored at a pose, along a fixed bearing, at an inverse depth."""

    def setUp(self):
        self.anchor = Pose3(Rot3.RzRyRx(0.1, -0.2, 0.3), Point3(1.0, 2.0, 0.5))
        self.observer = Pose3(Rot3.RzRyRx(-0.3, 0.1, 0.2), Point3(2.5, 1.0, 0.7))
        self.bearing = Unit3(np.array([0.3, -0.5, 0.8]))
        self.model = gtsam.noiseModel.Isotropic.Sigma(2, 0.01)

    def test_agrees_with_point_geometry(self):
        """The factor predicts the bearing GTSAM's point geometry gives, with zero error."""
        rho = 0.25
        landmark = InverseDepthBearingFactor.Landmark(self.anchor, self.bearing, rho)
        measured = bearing(self.observer, landmark)
        factor = InverseDepthBearingFactor(X(0), X(1), L(0), self.bearing, measured, self.model)
        self.gtsamAssertEquals(factor.landmark(self.anchor, rho), landmark, 1e-9)
        self.gtsamAssertEquals(factor.predict(self.anchor, self.observer, rho), measured, 1e-9)
        np.testing.assert_allclose(factor.evaluateError(self.anchor, self.observer, rho),
                                   np.zeros(2), atol=1e-9)
        # 2 degrees off: the error's norm is the sine of the angle
        off = Rot3.AxisAngle(Point3(measured.basis()[:, 0]), np.radians(2.0)).rotate(measured)
        factor = InverseDepthBearingFactor(X(0), X(1), L(0), self.bearing, off, self.model)
        err = factor.evaluateError(self.anchor, self.observer, rho)
        self.assertAlmostEqual(np.linalg.norm(err), np.sin(np.radians(2.0)), places=9)

    def test_sensor_poses(self):
        """Bearings from lenses displaced from the body poses."""
        anchor_sensor = Pose3(Rot3(), Point3(0.03, 0.0, 0.0))
        observer_sensor = Pose3(Rot3(), Point3(-0.03, 0.0, 0.01))
        rho = 0.25
        landmark = InverseDepthBearingFactor.Landmark(
            self.anchor.compose(anchor_sensor), self.bearing, rho)
        measured = bearing(self.observer.compose(observer_sensor), landmark)
        factor = InverseDepthBearingFactor(X(0), X(1), L(0), self.bearing, measured, self.model,
                                           anchor_sensor, observer_sensor)
        np.testing.assert_allclose(factor.evaluateError(self.anchor, self.observer, rho),
                                   np.zeros(2), atol=1e-9)
        bare = InverseDepthBearingFactor(X(0), X(1), L(0), self.bearing, measured, self.model)
        self.assertGreater(np.linalg.norm(bare.evaluateError(self.anchor, self.observer, rho)),
                           1e-3)

    def test_optimize(self):
        """Inverse depths and a free pose are recovered; inverse depths eliminated first."""
        poses = [self.anchor, self.observer,
                 Pose3(Rot3.RzRyRx(0.2, 0.05, -0.1), Point3(0.5, 3.0, 1.2))]
        bearings = [Unit3(np.array(b)) for b in
                    [[0.3, -0.5, 0.8], [-0.7, 0.2, 0.1], [0.1, 0.9, -0.4], [0.5, 0.5, 0.5]]]
        depths = [3.0, 6.0, 2.0, 10.0]
        graph, initial = gtsam.NonlinearFactorGraph(), gtsam.Values()
        prior = gtsam.noiseModel.Isotropic.Sigma(6, 1e-3)
        graph.addPriorPose3(X(0), poses[0], prior)
        graph.addPriorPose3(X(1), poses[1], prior)
        for k, (b, depth) in enumerate(zip(bearings, depths)):
            landmark = InverseDepthBearingFactor.Landmark(poses[0], b, 1.0 / depth)
            for j in (1, 2):
                graph.add(InverseDepthBearingFactor(X(0), X(j), L(k), b,
                                                    bearing(poses[j], landmark), self.model))
            initial.insert(L(k), 0.5 / depth)
        initial.insert(X(0), poses[0])
        initial.insert(X(1), poses[1])
        initial.insert(X(2), poses[2].retract(np.array([0.02, -0.01, 0.03, 0.2, -0.1, 0.15])))
        params = gtsam.LevenbergMarquardtParams()
        params.setOrdering(gtsam.Ordering.ColamdConstrainedFirstNonlinearFactorGraph(
            graph, [L(k) for k in range(len(depths))]))
        result = gtsam.LevenbergMarquardtOptimizer(graph, initial, params).optimize()
        for k, depth in enumerate(depths):
            self.assertAlmostEqual(result.atDouble(L(k)), 1.0 / depth, places=6)
        self.gtsamAssertEquals(result.atPose3(X(2)), poses[2], 1e-5)
        self.assertAlmostEqual(graph.error(result), 0.0, places=9)


if __name__ == "__main__":
    unittest.main()
