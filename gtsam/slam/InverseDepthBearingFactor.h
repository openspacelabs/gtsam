/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file InverseDepthBearingFactor.h
 * @brief A bearing measurement of a landmark that is parametrized by its inverse depth along a
 * fixed bearing from the pose it was first seen from.
 * @author OpenSpace
 */

#pragma once

#include <gtsam/base/Matrix.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Unit3.h>
#include <gtsam/nonlinear/NonlinearFactor.h>

#include <optional>

namespace gtsam {

/**
 * A landmark is parametrized by the pose it was first seen from (its anchor), the unit bearing
 * along which it was seen in the anchor's frame (fixed, taken from that first observation), and
 * its inverse depth along that bearing, which is zero for a point at infinity and otherwise
 * positive. The factor ties the anchor pose, another observer pose and the inverse depth through
 * the unit bearing the landmark is seen along from the observer, in the observer's frame:
 *
 *   predicted direction  d = R_j^T (R_h b_h + rho (p_h - p_j))
 *   error                the measured bearing's tangent-plane coordinates of d/|d|
 *                        (Unit3::errorVector), zero when the directions coincide.
 *
 * Bearings are directions in the body frame of a spherical (omnidirectional) camera, so there is
 * no image plane for a point to leave and no cheirality to check. Landmarks seen with little or
 * no parallax stay well conditioned at inverse depth near zero, where they constrain rotation
 * only and the Jacobians against the positions vanish; this is the anchored inverse-depth
 * parametrization of Civera, Davison and Montiel (2008), with the bearing held fixed as in
 * VINS-Mono. The inverse depth is not kept non-negative by the factor itself.
 *
 * Each bearing may come from a sensor displaced from its body pose (a lens of a multi-lens
 * rig): the optional sensor poses, in the body frame, are composed onto the anchor and the
 * observer, and the bearings are then in the sensors' frames.
 *
 * An inverse depth seen over a short baseline carries little information. Eliminate the
 * inverse depths before the poses (Ordering::ColamdConstrainedFirst): the Cholesky elimination
 * declares a clique underconstrained when the last two pivots of its frontal block differ by
 * more than 2^12, which such a scalar does next to a pose held by a tight prior.
 * @ingroup slam
 */
class GTSAM_EXPORT InverseDepthBearingFactor
    : public NoiseModelFactorN<Pose3, Pose3, double> {
 protected:
  Unit3 anchorBearing_;  ///< the landmark's bearing in the anchor's (sensor) frame, fixed
  Unit3 measured_;       ///< the bearing measured from the observer, in its (sensor) frame
  std::optional<Pose3> anchorSensor_;    ///< the anchor's sensor pose in its body frame
  std::optional<Pose3> observerSensor_;  ///< the observer's sensor pose in its body frame

 public:
  /// Directions shorter than this (the observer standing on the landmark) carry no bearing.
  static constexpr double kMinDirectionNorm = 1e-12;

  typedef NoiseModelFactorN<Pose3, Pose3, double> Base;
  typedef InverseDepthBearingFactor This;
  typedef std::shared_ptr<This> shared_ptr;

  // Provide access to the Matrix& version of evaluateError:
  using Base::evaluateError;

  /// Default constructor, for serialization
  InverseDepthBearingFactor() : anchorBearing_(0, 0, 1), measured_(0, 0, 1) {}

  /**
   * Constructor
   * @param anchorKey key of the pose the landmark is anchored at (its first observation)
   * @param observerKey key of the pose the bearing was measured from
   * @param inverseDepthKey key of the landmark's inverse depth (a double)
   * @param anchorBearing the landmark's bearing in the anchor's frame
   * @param measured the bearing measured from the observer, in its frame
   * @param model noise model on the two tangent-plane error coordinates
   * @param anchorSensor the anchor's sensor pose in its body frame (default identity)
   * @param observerSensor the observer's sensor pose in its body frame (default identity)
   */
  InverseDepthBearingFactor(Key anchorKey, Key observerKey, Key inverseDepthKey,
                            const Unit3& anchorBearing, const Unit3& measured,
                            const SharedNoiseModel& model,
                            std::optional<Pose3> anchorSensor = {},
                            std::optional<Pose3> observerSensor = {})
      : Base(model, anchorKey, observerKey, inverseDepthKey),
        anchorBearing_(anchorBearing),
        measured_(measured),
        anchorSensor_(anchorSensor),
        observerSensor_(observerSensor) {}

  ~InverseDepthBearingFactor() override {}

  /// @return a deep copy of this factor
  NonlinearFactor::shared_ptr clone() const override {
    return std::static_pointer_cast<NonlinearFactor>(NonlinearFactor::shared_ptr(new This(*this)));
  }

  /// print
  void print(const std::string& s = "",
             const KeyFormatter& keyFormatter = DefaultKeyFormatter) const override {
    std::cout << s << "InverseDepthBearingFactor\n";
    anchorBearing_.print("  anchor bearing: ");
    measured_.print("  measured: ");
    if (anchorSensor_) anchorSensor_->print("  anchor sensor pose in body frame: ");
    if (observerSensor_) observerSensor_->print("  observer sensor pose in body frame: ");
    Base::print("", keyFormatter);
  }

  /// Two optional sensor poses are equal when both are absent or both present and equal
  static bool sensorsEqual(const std::optional<Pose3>& a, const std::optional<Pose3>& b,
                           double tol) {
    return (!a && !b) || (a && b && a->equals(*b, tol));
  }

  /// equals
  bool equals(const NonlinearFactor& p, double tol = 1e-9) const override {
    const This* e = dynamic_cast<const This*>(&p);
    return e && Base::equals(p, tol) && anchorBearing_.equals(e->anchorBearing_, tol) &&
           measured_.equals(e->measured_, tol) &&
           sensorsEqual(anchorSensor_, e->anchorSensor_, tol) &&
           sensorsEqual(observerSensor_, e->observerSensor_, tol);
  }

  /**
   * The landmark as a world point: the anchor's position plus the bearing scaled by the depth.
   * Undefined (infinite) at zero inverse depth.
   */
  static Point3 Landmark(const Pose3& anchor, const Unit3& anchorBearing, double inverseDepth) {
    return anchor.translation() +
           anchor.rotation().rotate(Point3(anchorBearing.unitVector())) / inverseDepth;
  }

  /// The anchor's sensor pose in the world (the anchor pose itself without a sensor pose)
  Pose3 anchorSensorPose(const Pose3& anchor) const {
    return anchorSensor_ ? anchor.compose(*anchorSensor_) : anchor;
  }

  /// The observer's sensor pose in the world (the observer pose itself without a sensor pose)
  Pose3 observerSensorPose(const Pose3& observer) const {
    return observerSensor_ ? observer.compose(*observerSensor_) : observer;
  }

  /// The landmark as a world point, from the body poses (sensor poses applied)
  Point3 landmark(const Pose3& anchor, double inverseDepth) const {
    return Landmark(anchorSensorPose(anchor), anchorBearing_, inverseDepth);
  }

  /// The unnormalized direction to the landmark in the observer's sensor frame, from the
  /// sensor poses in the world
  Vector3 direction(const Pose3& anchorSensorPose, const Pose3& observerSensorPose,
                    double inverseDepth) const {
    const Vector3 inWorld =
        anchorSensorPose.rotation().rotate(Point3(anchorBearing_.unitVector())) +
        inverseDepth * (anchorSensorPose.translation() - observerSensorPose.translation());
    return observerSensorPose.rotation().unrotate(Point3(inWorld));
  }

  /// The bearing, in the observer's sensor frame, the landmark is predicted to be seen along,
  /// from the body poses
  Unit3 predict(const Pose3& anchor, const Pose3& observer, double inverseDepth) const {
    return Unit3(direction(anchorSensorPose(anchor), observerSensorPose(observer), inverseDepth));
  }

  /// Evaluate the error and, optionally, its derivatives against the anchor pose, the
  /// observer pose and the inverse depth
  Vector evaluateError(const Pose3& anchorBody, const Pose3& observerBody,
                       const double& inverseDepth, OptionalMatrixType H1, OptionalMatrixType H2,
                       OptionalMatrixType H3) const override {
    // The sensor poses in the world, and how they move with their body poses
    Matrix66 H_anchor, H_observer;
    const Pose3 anchor =
        anchorSensor_ ? anchorBody.compose(*anchorSensor_, H1 ? &H_anchor : nullptr) : anchorBody;
    const Pose3 observer = observerSensor_
                               ? observerBody.compose(*observerSensor_, H2 ? &H_observer : nullptr)
                               : observerBody;
    const Matrix3 Rh = anchor.rotation().matrix();
    const Matrix3 RjT = observer.rotation().matrix().transpose();
    const Vector3 bh = anchorBearing_.unitVector();
    const Vector3 baseline = anchor.translation() - observer.translation();
    const Vector3 d = RjT * (Rh * bh + inverseDepth * baseline);
    if (d.norm() < kMinDirectionNorm) {
      if (H1) *H1 = Matrix::Zero(2, 6);
      if (H2) *H2 = Matrix::Zero(2, 6);
      if (H3) *H3 = Matrix::Zero(2, 1);
      return Vector2::Zero();
    }
    const bool wantH = H1 || H2 || H3;
    Matrix23 H_u_d;
    const Unit3 u = Unit3::FromPoint3(Point3(d), wantH ? &H_u_d : nullptr);
    Matrix22 H_e_u;
    const Vector2 e = measured_.errorVector(u, {}, wantH ? &H_e_u : nullptr);
    if (wantH) {
      const Matrix23 H_e_d = H_e_u * H_u_d;
      if (H1) {
        // The anchor retracts as (R_h exp(w), p_h + R_h v): d(R_h b_h)/dw = -R_h [b_h]x and
        // dp_h/dv = R_h.
        Matrix36 D;
        D << -RjT * Rh * skewSymmetric(bh), inverseDepth * RjT * Rh;
        *H1 = anchorSensor_ ? Matrix(H_e_d * D * H_anchor) : Matrix(H_e_d * D);
      }
      if (H2) {
        // The observer retracts as (R_j exp(w), p_j + R_j v): R_j^T becomes (I - [w]x) R_j^T, so
        // dd/dw = [d]x, and dd/dv = -rho R_j^T R_j = -rho I.
        Matrix36 D;
        D << skewSymmetric(d), -inverseDepth * Matrix3::Identity();
        *H2 = observerSensor_ ? Matrix(H_e_d * D * H_observer) : Matrix(H_e_d * D);
      }
      if (H3) *H3 = H_e_d * (RjT * baseline);
    }
    return e;
  }

  /// The landmark's bearing in the anchor's frame
  const Unit3& anchorBearing() const { return anchorBearing_; }

  /// The bearing measured from the observer
  const Unit3& measured() const { return measured_; }

  /// The anchor's sensor pose in its body frame, if any
  const std::optional<Pose3>& anchorSensor() const { return anchorSensor_; }

  /// The observer's sensor pose in its body frame, if any
  const std::optional<Pose3>& observerSensor() const { return observerSensor_; }

 private:
#if GTSAM_ENABLE_BOOST_SERIALIZATION
  /// Serialization function
  friend class boost::serialization::access;
  template <class ARCHIVE>
  void serialize(ARCHIVE& ar, const unsigned int /*version*/) {
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(Base);
    ar& BOOST_SERIALIZATION_NVP(anchorBearing_);
    ar& BOOST_SERIALIZATION_NVP(measured_);
    ar& BOOST_SERIALIZATION_NVP(anchorSensor_);
    ar& BOOST_SERIALIZATION_NVP(observerSensor_);
  }
#endif
};

/// traits
template <>
struct traits<InverseDepthBearingFactor> : public Testable<InverseDepthBearingFactor> {};

}  // namespace gtsam
