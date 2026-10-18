"""
GTSAM Copyright 2010-2019, Georgia Tech Research Corporation,
Atlanta, Georgia 30332-0415
All Rights Reserved

See LICENSE for the license information

DirectionFactor unit tests.
Author: OpenSpace
"""

import unittest

import gtsam
import numpy as np
from gtsam import DirectionFactor, Pose3, Rot3
from gtsam.symbol_shorthand import A, X
from gtsam.utils.test_case import GtsamTestCase


def edge_normals(R: Rot3, building: float, n: int, seed: int) -> np.ndarray:
    """Body-frame great-circle normals of horizontal edges along a building's two axes."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        a = building + (np.pi / 2 if rng.uniform() < 0.5 else 0.0)
        d = np.array([np.cos(a), np.sin(a), 0.0])
        p = np.r_[rng.uniform(-4, 4, 2), rng.choice([-1, 1]) * rng.uniform(0.5, 1.5)]
        nw = np.cross(p, d)
        nw /= np.linalg.norm(nw)
        out.append(R.inverse().matrix() @ nw)
    return np.array(out)


class TestDirectionFactor(GtsamTestCase):
    def test_measures_the_building_and_turns_with_the_frame(self):
        building = np.radians(20.0)
        R = Rot3.RzRyRx(0.05, -0.03, 0.7)
        normals = edge_normals(R, building, 40, 1)
        model = gtsam.noiseModel.Isotropic.Sigma(1, np.radians(1.0))
        factor = DirectionFactor(X(0), A(0), normals, np.full(40, 0.1), model)
        self.assertAlmostEqual(DirectionFactor.wrap90(factor.meanDirection(R) - building), 0.0, places=6)
        self.assertAlmostEqual(factor.evaluateError(Pose3(R, np.zeros(3)), building)[0], 0.0, places=6)
        turn = np.radians(5.0)
        turned = Pose3(Rot3.Rz(turn) * R, np.zeros(3))
        self.assertAlmostEqual(factor.evaluateError(turned, building)[0], turn, places=6)

    def test_optimizes_the_heading_and_the_direction(self):
        building = np.radians(20.0)
        R0, R1 = Rot3.RzRyRx(0.02, 0.01, 0.4), Rot3.RzRyRx(-0.01, 0.03, 2.1)
        model = gtsam.noiseModel.Isotropic.Sigma(1, np.radians(1.0))
        graph, initial = gtsam.NonlinearFactorGraph(), gtsam.Values()
        for k, R in enumerate((R0, R1)):
            graph.add(DirectionFactor(X(k), A(0), edge_normals(R, building, 50, k + 3), np.full(50, 0.1), model))
        sig = np.array([1e-4, 1e-4, 10.0, 1e-4, 1e-4, 1e-4])  # heading free, the rest held
        start = [Pose3(Rot3.Rz(np.radians(5.0)) * R0, np.zeros(3)), Pose3(Rot3.Rz(np.radians(-3.0)) * R1, np.array([1.0, 0, 0]))]
        for k, s in enumerate(start):
            graph.addPriorPose3(X(k), s, gtsam.noiseModel.Diagonal.Sigmas(sig))
            initial.insert(X(k), s)
        graph.add(gtsam.PriorFactorDouble(A(0), building, gtsam.noiseModel.Isotropic.Sigma(1, np.radians(0.5))))
        initial.insert(A(0), building)
        result = gtsam.LevenbergMarquardtOptimizer(graph, initial).optimize()
        self.gtsamAssertEquals(result.atPose3(X(0)).rotation(), R0, 0.01)
        self.gtsamAssertEquals(result.atPose3(X(1)).rotation(), R1, 0.01)
        self.assertAlmostEqual(result.atDouble(A(0)), building, delta=0.01)


if __name__ == "__main__":
    unittest.main()
