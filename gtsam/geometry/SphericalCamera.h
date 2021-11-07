/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file SphericalCamera.h
 * @brief Calibrated camera with spherical projection
 * @date Aug 26, 2021
 * @author Luca Carlone
 */

#pragma once

#include <gtsam/geometry/BearingRange.h>
#include <gtsam/geometry/Unit3.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/base/concepts.h>
#include <gtsam/base/Manifold.h>
#include <gtsam/base/ThreadsafeException.h>
#include <gtsam/dllexport.h>
#include <boost/serialization/nvp.hpp>

namespace gtsam {

class GTSAM_EXPORT EmptyCal {
 public:
  enum { dimension = 0 };
  EmptyCal(){}
  virtual ~EmptyCal() = default;
  using shared_ptr = boost::shared_ptr<EmptyCal>;
  void print(const std::string& s) const {
    std::cout << "empty calibration: " <<  s << std::endl;
  }
};

/**
 * A spherical camera class that has a Pose3 and measures bearing vectors
 * @addtogroup geometry
 * \nosubgrouping
 */
class GTSAM_EXPORT SphericalCamera {

 public:

  enum {
    dimension = 6
  };

  typedef Unit3 Measurement;
  typedef std::vector<Unit3> MeasurementVector;
  typedef EmptyCal CalibrationType;

 private:

  Pose3 pose_;  ///< 3D pose of camera

 protected:

  EmptyCal::shared_ptr emptyCal_;

 public:

  /// @}
  /// @name Standard Constructors
  /// @{

  /// Default constructor
  SphericalCamera()
      : pose_(Pose3::identity()),
        emptyCal_(boost::make_shared<EmptyCal>()) {
  }

  /// Constructor with pose
  explicit SphericalCamera(const Pose3& pose)
      : pose_(pose),
        emptyCal_(boost::make_shared<EmptyCal>()) {
  }

  /// Constructor with empty intrinsics (needed for smart factors)
  explicit SphericalCamera(const Pose3& pose,
                           const boost::shared_ptr<EmptyCal>& cal)
      : pose_(pose),
        emptyCal_(boost::make_shared<EmptyCal>()) {
  }

  /// @}
  /// @name Advanced Constructors
  /// @{
  explicit SphericalCamera(const Vector& v)
      : pose_(Pose3::Expmap(v)) {
  }

  /// Default destructor
  virtual ~SphericalCamera() = default;

  /// return shared pointer to calibration
  const boost::shared_ptr<EmptyCal>& sharedCalibration() const {
    return emptyCal_;
  }

  /// return calibration
  const EmptyCal& calibration() const {
    return *emptyCal_;
  }

  /// @}
  /// @name Testable
  /// @{

  /// assert equality up to a tolerance
  bool equals(const SphericalCamera &camera, double tol = 1e-9) const;

  /// print
  virtual void print(const std::string& s = "SphericalCamera") const;

  /// @}
  /// @name Standard Interface
  /// @{

  /// return pose, constant version
  const Pose3& pose() const {
    return pose_;
  }

  /// get rotation
  const Rot3& rotation() const {
    return pose_.rotation();
  }

  /// get translation
  const Point3& translation() const {
    return pose_.translation();
  }

  //  /// return pose, with derivative
  //  const Pose3& getPose(OptionalJacobian<6, 6> H) const;

  /// @}
  /// @name Transformations and measurement functions
  /// @{

  /// Project a point into the image and check depth
  std::pair<Unit3, bool> projectSafe(const Point3& pw) const;

  /** Project point into the image
   * (note: there is no CheiralityException for a spherical camera)
   * @param point 3D point in world coordinates
   * @return the intrinsic coordinates of the projected point
   */
  Unit3 project2(const Point3& pw, OptionalJacobian<2, 6> Dpose = boost::none,
                 OptionalJacobian<2, 3> Dpoint = boost::none) const;

  /** Project point into the image
   * (note: there is no CheiralityException for a spherical camera)
   * @param point 3D direction in world coordinates
   * @return the intrinsic coordinates of the projected point
   */
  Unit3 project2(const Unit3& pwu, OptionalJacobian<2, 6> Dpose = boost::none,
                 OptionalJacobian<2, 2> Dpoint = boost::none) const;

  /// backproject a 2-dimensional point to a 3-dimensional point at given depth
  Point3 backproject(const Unit3& p, const double depth) const;

  /// backproject point at infinity
  Unit3 backprojectPointAtInfinity(const Unit3& p) const;

  /** Project point into the image
   * (note: there is no CheiralityException for a spherical camera)
   * @param point 3D point in world coordinates
   * @return the intrinsic coordinates of the projected point
   */
  Unit3 project(const Point3& point, OptionalJacobian<2, 6> Dpose = boost::none,
                OptionalJacobian<2, 3> Dpoint = boost::none) const;

  /** Compute reprojection error for a given 3D point in world coordinates
   * @param point 3D point in world coordinates
   * @return the tangent space error between the projection and the measurement
   */
  Vector2 reprojectionError(const Point3& point, const Unit3& measured,
                            OptionalJacobian<2, 6> Dpose = boost::none,
                            OptionalJacobian<2, 3> Dpoint = boost::none) const;
  /// @}

  /// move a cameras according to d
  SphericalCamera retract(const Vector6& d) const {
    return SphericalCamera(pose().retract(d));
  }

  /// return canonical coordinate
  Vector6 localCoordinates(const SphericalCamera& p) const {
    return pose().localCoordinates(p.pose());
  }

  /// for Canonical
  static SphericalCamera identity() {
    return SphericalCamera(Pose3::identity());  // assumes that the default constructor is valid
  }

  /// for Linear Triangulation
  Matrix34 getCameraProjectionMatrix() const {
    return Matrix34(pose_.inverse().matrix().block(0, 0, 3, 4));
  }

  /// for Nonlinear Triangulation
  Vector defaultErrorWhenTriangulatingBehindCamera() const {
    return Eigen::Matrix<double,traits<Point2>::dimension,1>::Constant(0.0);
  }

  /// @deprecated
   size_t dim() const {
     return 6;
   }

   /// @deprecated
   static size_t Dim() {
     return 6;
   }

 private:

  /** Serialization function */
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int /*version*/) {
    ar & BOOST_SERIALIZATION_NVP(pose_);
  }
};
// end of class SphericalCamera

template<>
struct traits<SphericalCamera> : public internal::LieGroup<Pose3> {
};

template<>
struct traits<const SphericalCamera> : public internal::LieGroup<Pose3> {
};

}
