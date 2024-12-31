// #ifdef CURVES_WITH_PINOCCHIO_SUPPORT

#define EIGEN_RUNTIME_NO_MALLOC
#define BOOST_TEST_MODULE test_so3_smooth

#include <boost/test/included/unit_test.hpp>
#include <random>

#include "ndcurves/fwd.h"
#include "ndcurves/serialization/curves.hpp"
#include "ndcurves/so3_smooth.h"

using namespace ndcurves;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

double generateRandomNumber(double lower, double upper) {
  // Some random number
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(lower, upper);
  return dist(gen);
}

BOOST_AUTO_TEST_CASE(default_constructor) {
  SO3Smooth_t traj;
  BOOST_CHECK_EQUAL(traj.dim(), 3);
  BOOST_CHECK_EQUAL(traj.min(), 0.0);
  BOOST_CHECK_EQUAL(traj.max(), 1.0);
  BOOST_CHECK_EQUAL(traj.degree(), 5.0);
  BOOST_CHECK(
      ndcurves::matrix3_t::Identity().isApprox(traj.get_init_rotation()));
  BOOST_CHECK(
      ndcurves::matrix3_t::Identity().isApprox(traj.get_end_rotation()));
}

BOOST_AUTO_TEST_CASE(default_generate) {
  SO3Smooth_t traj;

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(false);
  traj.generate();
  BOOST_CHECK_EQUAL(traj.dim(), 3);
  BOOST_CHECK_EQUAL(traj.min(), 0.0);
  BOOST_CHECK_EQUAL(traj.max(), 1.0);
  BOOST_CHECK_EQUAL(traj.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj.get_init_rotation(), Eigen::Matrix3d::Identity());
  BOOST_CHECK_EQUAL(traj.get_end_rotation(), Eigen::Matrix3d::Identity());
  Eigen::internal::set_is_malloc_allowed(true);
}

BOOST_AUTO_TEST_CASE(from_quat_and_time_constructor) {
  // Create some start variables.
  quaternion_t init_quat = quaternion_t::UnitRandom();
  matrix3_t init_rot = init_quat.toRotationMatrix();
  quaternion_t end_quat = quaternion_t::UnitRandom();
  matrix3_t end_rot = end_quat.toRotationMatrix();
  double t_min = generateRandomNumber(0.0, 100.0);
  double t_max = t_min + generateRandomNumber(0.1, 10.0);

  // Create the object.
  SO3Smooth_t traj(init_quat, end_quat, t_min, t_max);

  // Test the object.
  BOOST_CHECK_EQUAL(traj.dim(), 3);
  BOOST_CHECK_EQUAL(traj.min(), t_min);
  BOOST_CHECK_EQUAL(traj.max(), t_max);
  BOOST_CHECK_EQUAL(traj.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj.get_init_rotation(), init_rot);
  BOOST_CHECK_EQUAL(traj.get_end_rotation(), end_rot);
}

BOOST_AUTO_TEST_CASE(from_quat_and_time_generate) {
  // Create some start variables.
  quaternion_t init_quat = quaternion_t::UnitRandom();
  matrix3_t init_rot = init_quat.normalized().toRotationMatrix();
  quaternion_t end_quat = quaternion_t::UnitRandom();
  matrix3_t end_rot = end_quat.toRotationMatrix();
  double t_min = generateRandomNumber(0.0, 100.0);
  double t_max = t_min + generateRandomNumber(0.1, 10.0);

  // Create the object.
  SO3Smooth_t traj;

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(false);

  // Generate the trajectory.
  traj.generate(init_quat, end_quat, t_min, t_max);

  // Test the object.
  BOOST_CHECK_EQUAL(traj.dim(), 3);
  BOOST_CHECK_EQUAL(traj.min(), t_min);
  BOOST_CHECK_EQUAL(traj.max(), t_max);
  BOOST_CHECK_EQUAL(traj.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj.get_init_rotation(), init_rot);
  BOOST_CHECK_EQUAL(traj.get_end_rotation(), end_rot);

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(true);
}

BOOST_AUTO_TEST_CASE(from_rot_and_time_constructor) {
  // Create some start variables.
  quaternion_t init_quat = quaternion_t::UnitRandom();
  matrix3_t init_rot = init_quat.toRotationMatrix();
  quaternion_t end_quat = quaternion_t::UnitRandom();
  matrix3_t end_rot = end_quat.toRotationMatrix();
  double t_min = generateRandomNumber(0.0, 100.0);
  double t_max = t_min + generateRandomNumber(0.1, 10.0);

  // Create the object.
  SO3Smooth_t traj(init_rot, end_rot, t_min, t_max);

  // Test the object.
  BOOST_CHECK_EQUAL(traj.dim(), 3);
  BOOST_CHECK_EQUAL(traj.min(), t_min);
  BOOST_CHECK_EQUAL(traj.max(), t_max);
  BOOST_CHECK_EQUAL(traj.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj.get_init_rotation(), init_rot);
  BOOST_CHECK_EQUAL(traj.get_end_rotation(), end_rot);
}

BOOST_AUTO_TEST_CASE(from_rot_and_time_generate) {
  // Create some start variables.
  quaternion_t init_quat = quaternion_t::UnitRandom();
  matrix3_t init_rot = init_quat.toRotationMatrix();
  quaternion_t end_quat = quaternion_t::UnitRandom();
  matrix3_t end_rot = end_quat.toRotationMatrix();
  double t_min = generateRandomNumber(0.0, 100.0);
  double t_max = t_min + generateRandomNumber(0.1, 10.0);

  // Create the object.
  SO3Smooth_t traj;

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(false);

  // Generate the trajectory.
  traj.generate(init_rot, end_rot, t_min, t_max);

  // Test the object.
  BOOST_CHECK_EQUAL(traj.dim(), 3);
  BOOST_CHECK_EQUAL(traj.min(), t_min);
  BOOST_CHECK_EQUAL(traj.max(), t_max);
  BOOST_CHECK_EQUAL(traj.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj.get_init_rotation(), init_rot);
  BOOST_CHECK_EQUAL(traj.get_end_rotation(), end_rot);

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(true);
}

BOOST_AUTO_TEST_CASE(from_quat_constructor) {
  // Create some start variables.
  quaternion_t init_quat = quaternion_t::UnitRandom();
  matrix3_t init_rot = init_quat.toRotationMatrix();
  quaternion_t end_quat = quaternion_t::UnitRandom();
  matrix3_t end_rot = end_quat.toRotationMatrix();

  // Create the object.
  SO3Smooth_t traj(init_quat, end_quat);

  // Test the object.
  BOOST_CHECK_EQUAL(traj.dim(), 3);
  BOOST_CHECK_EQUAL(traj.min(), 0.0);
  BOOST_CHECK_EQUAL(traj.max(), 1.0);
  BOOST_CHECK_EQUAL(traj.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj.get_init_rotation(), init_rot);
  BOOST_CHECK_EQUAL(traj.get_end_rotation(), end_rot);
}

BOOST_AUTO_TEST_CASE(from_quat_generate) {
  // Create some start variables.
  quaternion_t init_quat = quaternion_t::UnitRandom();
  matrix3_t init_rot = init_quat.toRotationMatrix();
  quaternion_t end_quat = quaternion_t::UnitRandom();
  matrix3_t end_rot = end_quat.toRotationMatrix();

  // Create the object.
  SO3Smooth_t traj;

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(false);

  // Generate the trajectory.
  traj.generate(init_quat, end_quat);

  // Test the object.
  BOOST_CHECK_EQUAL(traj.dim(), 3);
  BOOST_CHECK_EQUAL(traj.min(), 0.0);
  BOOST_CHECK_EQUAL(traj.max(), 1.0);
  BOOST_CHECK_EQUAL(traj.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj.get_init_rotation(), init_rot);
  BOOST_CHECK_EQUAL(traj.get_end_rotation(), end_rot);

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(true);
}

BOOST_AUTO_TEST_CASE(from_rot_constructor) {
  // Create some start variables.
  quaternion_t init_quat = quaternion_t::UnitRandom();
  matrix3_t init_rot = init_quat.toRotationMatrix();
  quaternion_t end_quat = quaternion_t::UnitRandom();
  matrix3_t end_rot = end_quat.toRotationMatrix();

  // Create the object.
  SO3Smooth_t traj(init_rot, end_rot);

  // Test the object.
  BOOST_CHECK_EQUAL(traj.dim(), 3);
  BOOST_CHECK_EQUAL(traj.min(), 0.0);
  BOOST_CHECK_EQUAL(traj.max(), 1.0);
  BOOST_CHECK_EQUAL(traj.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj.get_init_rotation(), init_rot);
  BOOST_CHECK_EQUAL(traj.get_end_rotation(), end_rot);
}

BOOST_AUTO_TEST_CASE(from_rot_generate) {
  // Create some start variables.
  quaternion_t init_quat = quaternion_t::UnitRandom();
  matrix3_t init_rot = init_quat.toRotationMatrix();
  quaternion_t end_quat = quaternion_t::UnitRandom();
  matrix3_t end_rot = end_quat.toRotationMatrix();

  // Create the object.
  SO3Smooth_t traj;

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(false);

  // Generate the trajectory.
  traj.generate(init_rot, end_rot);

  // Test the object.
  BOOST_CHECK_EQUAL(traj.dim(), 3);
  BOOST_CHECK_EQUAL(traj.min(), 0.0);
  BOOST_CHECK_EQUAL(traj.max(), 1.0);
  BOOST_CHECK_EQUAL(traj.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj.get_init_rotation(), init_rot);
  BOOST_CHECK_EQUAL(traj.get_end_rotation(), end_rot);

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(true);
}

BOOST_AUTO_TEST_CASE(copy_real_time) {
  // None real_time critical
  quaternion_t init_quat = quaternion_t::UnitRandom();
  matrix3_t init_rot = init_quat.toRotationMatrix();
  quaternion_t end_quat = quaternion_t::UnitRandom();
  matrix3_t end_rot = end_quat.toRotationMatrix();
  double t_min = generateRandomNumber(0.0, 100.0);
  double t_max = t_min + generateRandomNumber(0.1, 10.0);

  // Create 2 objects and 1 initialized to be copied.
  SO3Smooth_t traj1(init_quat, end_quat, t_min, t_max);

  // Construct and copy the traj1 object.
  SO3Smooth_t traj2(traj1);

  // Test both trajectories independently.
  BOOST_CHECK_EQUAL(traj1.dim(), 3);
  BOOST_CHECK_EQUAL(traj1.min(), t_min);
  BOOST_CHECK_EQUAL(traj1.max(), t_max);
  BOOST_CHECK_EQUAL(traj1.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj1.get_init_rotation(), init_rot);
  BOOST_CHECK_EQUAL(traj1.get_end_rotation(), end_rot);
  BOOST_CHECK_EQUAL(traj2.dim(), 3);
  BOOST_CHECK_EQUAL(traj2.min(), t_min);
  BOOST_CHECK_EQUAL(traj2.max(), t_max);
  BOOST_CHECK_EQUAL(traj2.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj2.get_init_rotation(), init_rot);
  BOOST_CHECK_EQUAL(traj2.get_end_rotation(), end_rot);
}

BOOST_AUTO_TEST_CASE(computation_check) {
  // Create some start variables.
  quaternion_t init_quat = quaternion_t::UnitRandom();
  matrix3_t init_rot = init_quat.toRotationMatrix();
  quaternion_t end_quat = quaternion_t::UnitRandom();
  matrix3_t end_rot = end_quat.toRotationMatrix();
  double t_min = generateRandomNumber(0.0, 100.0);
  double t_max = t_min + generateRandomNumber(0.1, 10.0);
  std::vector<matrix3_t, Eigen::aligned_allocator<matrix3_t>> traj_discr;
  double dt = 1e-4;
  std::size_t nb_sample = static_cast<std::size_t>((t_max - t_min) / dt);
  traj_discr.resize(nb_sample, matrix3_t::Zero());

  // Create the object.
  SO3Smooth_t traj;

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(false);

  // Generate the trajectory.
  traj.generate(init_rot, end_rot, t_min, t_max);

  // Evaluate the trajectory.
  for (std::size_t i = 0; i < traj_discr.size(); ++i) {
    traj_discr[i] = traj(t_min + static_cast<double>(i) * dt);
  }

  // Test the object.
  BOOST_CHECK_EQUAL(traj.dim(), 3);
  BOOST_CHECK_EQUAL(traj.min(), t_min);
  BOOST_CHECK_EQUAL(traj.max(), t_max);
  BOOST_CHECK_EQUAL(traj.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj.get_init_rotation(), init_rot);
  BOOST_CHECK_EQUAL(traj.get_end_rotation(), end_rot);
  BOOST_CHECK(traj(traj.min()).isApprox(init_rot));
  BOOST_CHECK(traj(traj.max()).isApprox(end_rot));

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(true);

  // Continuity test.
  for (std::size_t i = 1; i < traj_discr.size(); ++i) {
    matrix3_t error = traj_discr[i] - traj_discr[i - 1];
    for (Eigen::Index row = 0; row < error.rows(); ++row) {
      for (Eigen::Index col = 0; col < error.cols(); ++col) {
        BOOST_CHECK_LE(std::abs(error(row, col)), 1e-2);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(derivate_computation_check) {
  // Create some start variables.
  quaternion_t init_quat = quaternion_t::UnitRandom();
  matrix3_t init_rot = init_quat.toRotationMatrix();
  quaternion_t end_quat = quaternion_t::UnitRandom();
  matrix3_t end_rot = end_quat.toRotationMatrix();
  double t_min = generateRandomNumber(0.0, 100.0);
  double t_max = t_min + generateRandomNumber(0.1, 10.0);
  Eigen::MatrixXd traj_vel_discr;
  Eigen::MatrixXd traj_acc_discr;
  double dt = 1e-3;
  std::size_t nb_sample = static_cast<std::size_t>((t_max - t_min) / dt);
  traj_vel_discr.resize(3, nb_sample);
  traj_acc_discr.resize(3, nb_sample);

  // Create the object.
  SO3Smooth_t traj;

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(false);

  // Generate the trajectory.
  traj.generate(init_rot, end_rot, t_min, t_max);

  // Evaluate the trajectory.
  for (std::size_t i = 0; i < nb_sample; ++i) {
    double t = t_min + static_cast<double>(i) * dt;
    traj_vel_discr.col(i) = traj.derivate(t, 1);
    traj_acc_discr.col(i) = traj.derivate(t, 2);
  }

  // Test the object.
  BOOST_CHECK_EQUAL(traj.dim(), 3);
  BOOST_CHECK_EQUAL(traj.min(), t_min);
  BOOST_CHECK_EQUAL(traj.max(), t_max);
  BOOST_CHECK_EQUAL(traj.degree(), 5.0);
  BOOST_CHECK_EQUAL(traj.get_init_rotation(), init_rot);
  BOOST_CHECK_EQUAL(traj.get_end_rotation(), end_rot);
  BOOST_CHECK((traj(traj.min()) - init_rot).isMuchSmallerThan(1, -1e-8));
  BOOST_CHECK((traj(traj.max()) - end_rot).isMuchSmallerThan(1, -1e-8));
  BOOST_CHECK_EQUAL(traj.derivate(traj.min(), 1), point3_t::Zero());
  BOOST_CHECK_LE(traj.derivate(traj.max(), 1).norm(), 1e-3);

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(true);

  point3_t vel_max;
  vel_max << traj_vel_discr.row(0).maxCoeff(), traj_vel_discr.row(1).maxCoeff(),
      traj_vel_discr.row(2).maxCoeff();
  point3_t vel_min;
  vel_min << traj_vel_discr.row(0).minCoeff(), traj_vel_discr.row(1).minCoeff(),
      traj_vel_discr.row(2).minCoeff();
  point3_t vel_amplitude = (vel_max - vel_min).array().abs();
  point3_t acc_max;
  acc_max << traj_acc_discr.row(0).maxCoeff(), traj_acc_discr.row(1).maxCoeff(),
      traj_acc_discr.row(2).maxCoeff();
  point3_t acc_min;
  acc_min << traj_acc_discr.row(0).minCoeff(), traj_acc_discr.row(1).minCoeff(),
      traj_acc_discr.row(2).minCoeff();
  point3_t acc_amplitude = (acc_max - acc_min).array().abs();

  // Continuity test.
  for (std::size_t i = 1; i < nb_sample; ++i) {
    point3_t error = traj_vel_discr.col(i) - traj_vel_discr.col(i - 1);
    for (Eigen::Index j = 0; j < error.size(); ++j) {
      BOOST_CHECK_LE(std::abs(error(j)), 1e-0 * vel_amplitude(j));
    }
  }

  for (std::size_t i = 1; i < nb_sample; ++i) {
    point3_t error = traj_acc_discr.col(i) - traj_acc_discr.col(i - 1);
    for (Eigen::Index j = 0; j < error.size(); ++j) {
      BOOST_CHECK_LE(std::abs(error(j)), 1e-0 * acc_amplitude(j));
    }
  }
}

BOOST_AUTO_TEST_CASE(serialization) {
  std::string fileName("fileTest_constant");

  // Create some start variables.
  quaternion_t init_quat = quaternion_t::UnitRandom();
  matrix3_t init_rot = init_quat.toRotationMatrix();
  quaternion_t end_quat = quaternion_t::UnitRandom();
  matrix3_t end_rot = end_quat.toRotationMatrix();
  double t_min = generateRandomNumber(0.0, 100.0);
  double t_max = t_min + generateRandomNumber(0.1, 10.0);

  // Create the object.
  SO3Smooth_t c(init_quat, end_quat, t_min, t_max);

  c.saveAsText<SO3Smooth_t>(fileName + ".txt");
  c.saveAsXML<SO3Smooth_t>(fileName + ".xml", "constant");
  c.saveAsBinary<SO3Smooth_t>(fileName);

  SO3Smooth_t c_txt, c_xml, c_binary;
  c_txt.loadFromText<SO3Smooth_t>(fileName + ".txt");
  c_xml.loadFromXML<SO3Smooth_t>(fileName + ".xml", "constant");
  c_binary.loadFromBinary<SO3Smooth_t>(fileName);

  BOOST_CHECK(c == c_txt);
  BOOST_CHECK(c == c_xml);
  BOOST_CHECK(c == c_binary);
}

BOOST_AUTO_TEST_SUITE_END()

// #endif  // CURVES_WITH_PINOCCHIO_SUPPORT
