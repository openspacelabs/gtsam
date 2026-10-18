/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file DirectionFactor.h
 * @brief A keyframe's horizontal edges tie its rotation to a building direction that is a
 * variable: the Atlanta/Manhattan world's directions as state.
 * @author OpenSpace
 */

#pragma once

#include <gtsam/base/Matrix.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/nonlinear/NonlinearFactor.h>

#include <cmath>

namespace gtsam {

/**
 * Straight edges seen in a keyframe (walls, joists, door frames) are each known by the normal
 * of their great circle in the body frame. Under the keyframe's rotation R, an edge that is
 * horizontal runs along (R n) x z, whose angle in the horizontal plane is a direction of the
 * building, modulo 90 degrees in a rectilinear building. This factor ties the pose's rotation to
 * a direction variable (radians, modulo pi/2) through the weighted circular mean of the
 * keyframe's edges' directions:
 *
 *   a_k   = atan2(-(R n_k)_x, (R n_k)_y)                     the edge's horizontal direction
 *   mean  = arg( sum_k w_k exp(4 i a_k) ) / 4                 the keyframe's direction
 *   error = wrap90(mean - direction)                          in (-pi/4, pi/4]
 *
 * The edges given should already be the ones on that direction (within a few degrees of its
 * peak) with their tilts in range: nearly vertical edges say nothing about heading, and edges
 * near camera height have their direction set by noise. Derivatives are analytic.
 * @ingroup slam
 */
class GTSAM_EXPORT DirectionFactor : public NoiseModelFactorN<Pose3, double> {
 protected:
  Matrix normals_;  ///< N x 3 unit great-circle normals of the edges, body frame
  Vector weights_;  ///< N weights (the edges' lengths)

 public:
  typedef NoiseModelFactorN<Pose3, double> Base;
  typedef DirectionFactor This;
  typedef std::shared_ptr<This> shared_ptr;

  using Base::evaluateError;

  /// Default constructor, for serialization
  DirectionFactor() {}

  /**
   * Constructor
   * @param poseKey key of the keyframe's pose
   * @param directionKey key of the building direction (a double, radians modulo pi/2)
   * @param normals N x 3 unit great-circle normals of the keyframe's edges on that direction
   * @param weights N weights for the edges
   * @param model noise model on the one error coordinate (radians)
   */
  DirectionFactor(Key poseKey, Key directionKey, const Matrix& normals, const Vector& weights,
                  const SharedNoiseModel& model)
      : Base(model, poseKey, directionKey), normals_(normals), weights_(weights) {
    if (normals_.cols() != 3 || normals_.rows() != weights_.size())
      throw std::invalid_argument("DirectionFactor: normals must be N x 3 with N weights");
  }

  ~DirectionFactor() override {}

  /// @return a deep copy of this factor
  NonlinearFactor::shared_ptr clone() const override {
    return std::static_pointer_cast<NonlinearFactor>(NonlinearFactor::shared_ptr(new This(*this)));
  }

  /// print
  void print(const std::string& s = "",
             const KeyFormatter& keyFormatter = DefaultKeyFormatter) const override {
    std::cout << s << "DirectionFactor, " << normals_.rows() << " edges\n";
    Base::print("", keyFormatter);
  }

  /// equals
  bool equals(const NonlinearFactor& p, double tol = 1e-9) const override {
    const This* e = dynamic_cast<const This*>(&p);
    return e && Base::equals(p, tol) && equal_with_abs_tol(normals_, e->normals_, tol) &&
           equal_with_abs_tol(weights_, e->weights_, tol);
  }

  /// An angle wrapped to (-pi/4, pi/4]
  static double wrap90(double a) { return std::atan2(std::sin(4 * a), std::cos(4 * a)) / 4; }

  /**
   * The keyframe's direction under rotation R (radians, modulo pi/2) and, optionally, its
   * derivative against the rotation's tangent (right-multiplied, as Rot3 retracts).
   */
  double meanDirection(const Rot3& R, OptionalJacobian<1, 3> H = {}) const {
    const Matrix3 Rm = R.matrix();
    double C = 0.0, S = 0.0;
    Vector3 dC = Vector3::Zero(), dS = Vector3::Zero();
    for (int k = 0; k < normals_.rows(); ++k) {
      const Vector3 n = normals_.row(k).transpose();
      const Vector3 nw = Rm * n;
      const double r2 = nw.x() * nw.x() + nw.y() * nw.y();
      if (r2 < 1e-18) continue;  // a vertical edge: no direction
      const double a = std::atan2(-nw.x(), nw.y());
      const double w = weights_(k);
      C += w * std::cos(4 * a);
      S += w * std::sin(4 * a);
      if (H) {
        // d a / d nw = (-nw_y, nw_x, 0) / r2;  d nw / d omega = -R [n]x  (R exp(omega) n)
        const Vector3 da_dnw(-nw.y() / r2, nw.x() / r2, 0.0);
        const Vector3 da = (-Rm * skewSymmetric(n)).transpose() * da_dnw;
        dC += -4 * w * std::sin(4 * a) * da;
        dS += 4 * w * std::cos(4 * a) * da;
      }
    }
    const double mean = std::atan2(S, C) / 4;
    if (H) {
      const double den = C * C + S * S;
      if (den > 0) {
        const Vector3 g = (C * dS - S * dC) / (4 * den);
        *H = g.transpose();
      } else {
        *H = Matrix13::Zero();
      }
    }
    return mean;
  }

  /// Evaluate the error and, optionally, its derivatives against the pose and the direction
  Vector evaluateError(const Pose3& pose, const double& direction, OptionalMatrixType H1,
                       OptionalMatrixType H2) const override {
    Matrix13 Hr;
    const double mean = meanDirection(pose.rotation(), H1 ? &Hr : nullptr);
    if (H1) {
      Matrix16 J = Matrix16::Zero();
      J.leftCols<3>() = Hr;  // Pose3's tangent: rotation first, then translation
      *H1 = J;
    }
    if (H2) *H2 = -Matrix11::Identity();
    return Vector1(wrap90(mean - direction));
  }

  /// The edges' normals
  const Matrix& normals() const { return normals_; }

  /// The edges' weights
  const Vector& weights() const { return weights_; }

 private:
#if GTSAM_ENABLE_BOOST_SERIALIZATION
  /// Serialization function
  friend class boost::serialization::access;
  template <class ARCHIVE>
  void serialize(ARCHIVE& ar, const unsigned int /*version*/) {
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(Base);
    ar& BOOST_SERIALIZATION_NVP(normals_);
    ar& BOOST_SERIALIZATION_NVP(weights_);
  }
#endif
};

/// traits
template <>
struct traits<DirectionFactor> : public Testable<DirectionFactor> {};

}  // namespace gtsam
