import os
import pickle
import unittest
from math import sqrt

import numpy as np
from numpy import array, array_equal, isclose, random, zeros
from numpy.linalg import norm

from ndcurves import (
    CURVES_WITH_PINOCCHIO_SUPPORT,
    Quaternion,
    SE3Curve,
    SO3Linear,
    bezier,
    bezier3,
    convert_to_bezier,
    convert_to_hermite,
    convert_to_polynomial,
    cubic_hermite_spline,
    curve_constraints,
    exact_cubic,
    piecewise,
    piecewise3,
    piecewise_SE3,
    polynomial,
)

if CURVES_WITH_PINOCCHIO_SUPPORT:
    from pinocchio import SE3, Motion


class TestCurves(unittest.TestCase):
    # def print_str(self, inStr):
    #   print inStr
    #   return

    def compareCurves(self, c1, c2):
        t_min = c1.min()
        t_max = c1.max()
        self.assertEqual(
            t_min, c2.min(), "c1 min : " + str(t_min) + " ; c2 min : " + str(c2.min())
        )
        self.assertEqual(
            t_max, c2.max(), "c1 max : " + str(t_max) + " ; c2 max : " + str(c2.max())
        )
        self.assertTrue(
            norm(c1.derivate(t_min, 1) - c2.derivate(t_min, 1)) < 1e-10,
            "dc1(tmin) = "
            + str(c1.derivate(t_min, 1))
            + " ; dc2(tmin) = "
            + str(c2.derivate(t_min, 1)),
        )
        self.assertTrue(
            norm(c1.derivate(t_max, 1) - c2.derivate(t_max, 1)) < 1e-10,
            "dc1(tmax) = "
            + str(c1.derivate(t_max, 1))
            + " ; dc2(tmax) = "
            + str(c2.derivate(t_max, 1)),
        )
        t = t_min
        while t < t_max:
            self.assertTrue(
                norm(c1(t) - c2(t)) < 1e-10,
                " at t = " + str(t) + " c1 = " + str(c1(t)) + " ; c2 = " + str(c2(t)),
            )
            t = t + 0.01

    def test_bezier(self):
        print("test_bezier")
        # To test :
        # - Functions : constructor, min, max, derivate,
        #               compute_derivate, compute_primitive
        # - Variables : degree, nbWayPoints
        __EPS = 1e-6
        waypoints = array([[1.0, 2.0, 3.0]]).T
        a = bezier(waypoints, 0.0, 2.0)
        t = 0.0
        while t < 2.0:
            self.assertTrue(norm(a(t) - array([1.0, 2.0, 3.0]).T) < __EPS)
            t += 0.1
        waypoints = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        # time_waypoints = array([0., 1.]).transpose()
        # Create bezier6 and bezier
        a = bezier(waypoints, 0.0, 3.0)
        # Test waypoints
        self.assertTrue(a.nbWaypoints == 2)
        for i in range(0, a.nbWaypoints):
            if i == 0:
                self.assertTrue((a.waypointAtIndex(0) == array([1.0, 2.0, 3.0])).all())
            elif i == 1:
                self.assertTrue((a.waypointAtIndex(1) == array([4.0, 5.0, 6.0])).all())

        a1 = a.elevate(1)
        for i in range(100):
            dt = float(i) / 100.0 * 3.0
            self.assertTrue(norm(a(dt) - a1(dt)) < __EPS)

        # arithmetic
        b = a + a1
        b += a1
        b = a - a1
        b -= a
        a1 *= 0.1
        a1 /= 0.1
        b = -a1
        c = a.cross(b)
        c(0)
        b += array([1.0, 2.0, 3.0])
        b -= array([1.0, 2.0, 3.0])
        b = a + array([1.0, 2.0, 3.0])
        b = a - array([1.0, 2.0, 3.0])

        # self.assertTrue((a.waypoints == waypoints).all())
        # Test : Degree, min, max, derivate
        # self.print_str(("test 1")
        self.assertEqual(a.degree, a.nbWaypoints - 1)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a(a.min()) == array([1.0, 2.0, 3.0])).all())
        a.derivate(0.4, 2)
        a = a.compute_derivate(100)
        prim = a.compute_primitive(1)
        # Check primitive and derivate - order 1
        for i in range(10):
            t = float(i) / 10.0
            self.assertTrue((a(t) == prim.derivate(t, 1)).all())
        self.assertTrue((prim(0) == array([0.0, 0.0, 0.0])).all())
        # Check primitive and derivate - order 2
        prim = a.compute_primitive(2)
        for i in range(10):
            t = float(i) / 10.0
            self.assertTrue((a(t) == prim.derivate(t, 2)).all())
        self.assertTrue((prim(0) == array([0.0, 0.0, 0.0])).all())
        # Create new bezier curve
        waypoints = array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
            ]
        ).transpose()
        a0 = bezier(waypoints)
        a1 = bezier(waypoints, 0.0, 3.0)
        prim0 = a0.compute_primitive(1)
        prim1 = a1.compute_primitive(1)
        # Check change in argument time_t of bezier
        for i in range(10):
            t = float(i) / 10.0
            self.assertTrue(norm(a0(t) - a1(3 * t)) < __EPS)
            self.assertTrue(
                norm(a0.derivate(t, 1) - a1.derivate(3 * t, 1) * 3.0) < __EPS
            )
            self.assertTrue(
                norm(a0.derivate(t, 2) - a1.derivate(3 * t, 2) * 9.0) < __EPS
            )
            self.assertTrue(norm(prim0(t) - prim1(t * 3) / 3.0) < __EPS)
        self.assertTrue((prim(0) == array([0.0, 0.0, 0.0])).all())
        # testing accessor to waypoints :
        wp_getter = a0.waypoints()
        self.assertEqual(wp_getter.shape[0], waypoints.shape[0])
        self.assertEqual(wp_getter.shape[1], waypoints.shape[1])
        self.assertTrue(array_equal(wp_getter, waypoints))
        # check that it return a copy:
        a0.waypoints()[1, 1] = -15.0
        self.assertEqual(a0.waypoints()[1, 1], waypoints[1, 1])
        # testing bezier with constraints
        c = curve_constraints(3)
        c.init_vel = array([[0.0, 1.0, 1.0]]).transpose()
        c.end_vel = array([[0.0, 1.0, 1.0]]).transpose()
        c.init_acc = array([[0.0, 1.0, -1.0]]).transpose()
        c.end_acc = array([[0.0, 100.0, 1.0]]).transpose()
        # Check derivate with constraints
        waypoints = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        a = bezier(waypoints, c)
        self.assertTrue(norm(a.derivate(0, 1) - c.init_vel) < 1e-10)
        self.assertTrue(norm(a.derivate(1, 2) - c.end_acc) < 1e-10)

        # Test serialization : bezier 3
        a.saveAsText("serialization_curve.test")
        # waypoints = array([[0,0,0,], [0,0,0,]]).transpose()
        b = bezier()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")

        # Bezier dim 4
        waypoints = array([[1.0, 2.0, 3.0, 4.0]]).T
        a = bezier(waypoints, 0.0, 2.0)
        # Test serialization : bezier of dim 4
        a.saveAsText("serialization_curve.test")
        # waypoints = array([[0,0,0,], [0,0,0,]]).transpose()
        b = bezier()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")

        a_pickled = pickle.dumps(a)
        a_from_pickle = pickle.loads(a_pickled)
        self.assertEqual(a_from_pickle, a)
        return

    def test_bezier3(self):
        print("test_bezier3")
        # To test :
        # - Functions : constructor, min, max, derivate,
        #               compute_derivate, compute_primitive
        # - Variables : degree, nbWayPoints
        __EPS = 1e-6
        waypoints = array([[1.0, 2.0, 3.0]]).T
        a = bezier3(waypoints, 0.0, 2.0)
        t = 0.0
        while t < 2.0:
            self.assertTrue(norm(a(t) - array([1.0, 2.0, 3.0]).T) < __EPS)
            t += 0.1
        waypoints = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        # time_waypoints = array([[0., 1.]]).transpose()
        # Create bezier6 and bezier
        a = bezier3(waypoints, 0.0, 3.0)
        a1 = a.elevate(1)
        b = bezier3(waypoints, 0.0, 3.0)
        b.elevateSelf(2)
        assert b.degree == a.degree + 2
        for i in range(100):
            dt = float(i) / 100.0 * 3.0
            self.assertTrue(norm(a(dt) - a1(dt)) < __EPS)
            self.assertTrue(norm(a(dt) - b(dt)) < __EPS)

        # arithmetic
        b = a + a1
        b += a1
        b = a - a1
        b -= a
        a1 *= 0.1
        a1 /= 0.1
        b = -a1
        c = a.cross(b)
        c(0)
        b += array([1.0, 2.0, 3.0])
        b -= array([1.0, 2.0, 3.0])
        b = a + array([1.0, 2.0, 3.0])
        b = a - array([1.0, 2.0, 3.0])

        # Test waypoints
        self.assertTrue(a.nbWaypoints == 2)
        for i in range(0, a.nbWaypoints):
            if i == 0:
                self.assertTrue((a.waypointAtIndex(0) == array([1.0, 2.0, 3.0])).all())
            elif i == 1:
                self.assertTrue((a.waypointAtIndex(1) == array([4.0, 5.0, 6.0])).all())
        # self.assertTrue((a.waypoints == waypoints).all())
        # Test : Degree, min, max, derivate
        # self.print_str(("test 1")
        self.assertEqual(a.degree, a.nbWaypoints - 1)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a(a.min()) == array([1.0, 2.0, 3.0])).all())
        a.derivate(0.4, 2)
        a = a.compute_derivate(100)
        prim = a.compute_primitive(1)
        # Check primitive and derivate - order 1
        for i in range(10):
            t = float(i) / 10.0
            self.assertTrue((a(t) == prim.derivate(t, 1)).all())
        self.assertTrue((prim(0) == array([0.0, 0.0, 0.0])).all())
        # Check primitive and derivate - order 2
        prim = a.compute_primitive(2)
        for i in range(10):
            t = float(i) / 10.0
            self.assertTrue((a(t) == prim.derivate(t, 2)).all())
        self.assertTrue((prim(0) == array([0.0, 0.0, 0.0])).all())
        # Create new bezier curve
        waypoints = array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
            ]
        ).transpose()
        a0 = bezier3(waypoints)
        a1 = bezier3(waypoints, 0.0, 3.0)
        prim0 = a0.compute_primitive(1)
        prim1 = a1.compute_primitive(1)
        # Check change in argument time_t of bezier
        for i in range(10):
            t = float(i) / 10.0
            self.assertTrue(norm(a0(t) - a1(3 * t)) < __EPS)
            self.assertTrue(
                norm(a0.derivate(t, 1) - a1.derivate(3 * t, 1) * 3.0) < __EPS
            )
            self.assertTrue(
                norm(a0.derivate(t, 2) - a1.derivate(3 * t, 2) * 9.0) < __EPS
            )
            self.assertTrue(norm(prim0(t) - prim1(t * 3) / 3.0) < __EPS)
        self.assertTrue((prim(0) == array([0.0, 0.0, 0.0])).all())
        # testing accessor to waypoints :
        wp_getter = a0.waypoints()
        self.assertEqual(wp_getter.shape[0], waypoints.shape[0])
        self.assertEqual(wp_getter.shape[1], waypoints.shape[1])
        self.assertTrue(array_equal(wp_getter, waypoints))
        # check that it return a copy:
        a0.waypoints()[1, 1] = -15.0
        self.assertEqual(a0.waypoints()[1, 1], waypoints[1, 1])
        # testing bezier with constraints
        c = curve_constraints(3)
        c.init_vel = array([[0.0, 1.0, 1.0]]).transpose()
        c.end_vel = array([[0.0, 1.0, 1.0]]).transpose()
        c.init_acc = array([[0.0, 1.0, -1.0]]).transpose()
        c.end_acc = array([[0.0, 100.0, 1.0]]).transpose()
        # Check derivate with constraints
        waypoints = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        a = bezier3(waypoints, c)
        self.assertTrue(norm(a.derivate(0, 1) - c.init_vel) < 1e-10)
        self.assertTrue(norm(a.derivate(1, 2) - c.end_acc) < 1e-10)

        # Test serialization : bezier 3
        a.saveAsText("serialization_curve.test")
        # waypoints = array([[0,0,0,], [0,0,0,]]).transpose()
        b = bezier3()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")

        # Bezier dim 4
        waypoints = array([[1.0, 2.0, 3.0, 4.0]]).T
        a = bezier(waypoints, 0.0, 2.0)
        # Test serialization : bezier of dim 4
        a.saveAsText("serialization_curve.test")
        # waypoints = array([[0,0,0,], [0,0,0,]]).transpose()
        b = bezier()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")
        a_pickled = pickle.dumps(a)
        a_from_pickle = pickle.loads(a_pickled)
        self.assertEqual(a_from_pickle, a)
        return

    def test_polynomial(self):
        print("test_polynomial")
        # To test :
        # - Functions : constructor, min, max, derivate, serialize, deserialize
        waypoints = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        a = polynomial(waypoints)  # Defined on [0.,1.]
        a = polynomial(waypoints, -1.0, 3.0)  # Defined on [-1.,3.]
        a.min()
        a.max()
        a(0.4)

        # arithmetic
        # waypoints2 = array([[1., 2., 3.], [4., 5., 6.], [4., 5., 6.]]).transpose()
        a1 = polynomial(waypoints, -1.0, 3.0)  # Defined on [-1.,3.]
        b = a + a1
        b += a1
        b = a - a1
        b -= a
        a1 *= 0.1
        a1 /= 0.1
        b = -a1
        c = a.cross(array([1.0, 2.0, 3.0]))
        c = a.cross(a)
        c(0)
        b += array([1.0, 2.0, 3.0])
        b -= array([1.0, 2.0, 3.0])
        b = a + array([1.0, 2.0, 3.0])
        b = a - array([1.0, 2.0, 3.0])

        # Test get coefficient at degree
        self.assertTrue((a.coeff() == waypoints).all())
        self.assertTrue((a.coeffAtDegree(0) == array([1.0, 2.0, 3.0])).all())
        self.assertTrue((a.coeffAtDegree(1) == array([4.0, 5.0, 6.0])).all())
        # Other tests
        self.assertTrue((a(a.min()) == array([1.0, 2.0, 3.0])).all())
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        a_derivated = a.compute_derivate(1)
        self.assertTrue((a.derivate(0.4, 1) == a_derivated(0.4)).all())
        # Test serialization
        a.saveAsText("serialization_curve.test")
        b = polynomial()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")
        a_pickled = pickle.dumps(a)
        a_from_pickle = pickle.loads(a_pickled)
        self.assertEqual(a_from_pickle, a)
        return

    def test_polynomial_from_boundary_condition(self):
        p0 = array([1.0, 3.0, -2.0])
        p1 = array([0.6, 2.0, 2.5])
        dp0 = array([-6.0, 2.0, -1.0])
        dp1 = array([10.0, 10.0, 10.0])
        ddp0 = array([1.0, -7.0, 4.5])
        ddp1 = array([6.0, -1.0, -4])
        min = 1.0
        max = 2.5
        polC0 = polynomial(p0, p1, min, max)
        self.assertEqual(polC0.min(), min)
        self.assertEqual(polC0.max(), max)
        # TODO: Why are thoso `.T[0]` needed ?
        self.assertTrue(array_equal(polC0((min + max) / 2.0), 0.5 * p0 + 0.5 * p1))
        polC1 = polynomial(p0, dp0, p1, dp1, min, max)
        self.assertEqual(polC1.min(), min)
        self.assertEqual(polC1.max(), max)
        self.assertTrue(isclose(polC1(min), p0).all())
        self.assertTrue(isclose(polC1(max), p1).all())
        self.assertTrue(isclose(polC1.derivate(min, 1), dp0).all())
        self.assertTrue(isclose(polC1.derivate(max, 1), dp1).all())
        polC2 = polynomial(p0, dp0, ddp0, p1, dp1, ddp1, min, max)
        self.assertEqual(polC2.min(), min)
        self.assertEqual(polC2.max(), max)
        self.assertTrue(isclose(polC2(min), p0).all())
        self.assertTrue(isclose(polC2(max), p1).all())
        self.assertTrue(isclose(polC2.derivate(min, 1), dp0).all())
        self.assertTrue(isclose(polC2.derivate(max, 1), dp1).all())
        self.assertTrue(isclose(polC2.derivate(min, 2), ddp0).all())
        self.assertTrue(isclose(polC2.derivate(max, 2), ddp1).all())
        # check that the exception are correctly raised :
        with self.assertRaises(ValueError):
            polC0 = polynomial(p0, p1, max, min)

        with self.assertRaises(ValueError):
            polC1 = polynomial(p0, dp0, p1, dp1, max, min)

        with self.assertRaises(ValueError):
            polC2 = polynomial(p0, dp0, ddp0, p1, dp1, ddp1, max, min)

    def test_cubic_hermite_spline(self):
        print("test_cubic_hermite_spline")
        points = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        tangents = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        time_points = array([[0.0, 1.0]]).transpose()
        a = cubic_hermite_spline(points, tangents, time_points)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        # Test serialization
        a.saveAsText("serialization_curve.test")
        b = cubic_hermite_spline()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")
        # test dim 4
        points = array([[1.0, 2.0, 3.0, 4.0], [4.0, 5.0, 6.0, 7.0]]).transpose()
        tangents = array([[1.0, 2.0, 3.0, 4.0], [4.0, 5.0, 6.0, 7.0]]).transpose()
        time_points = array([[0.0, 1.0]]).transpose()
        a = cubic_hermite_spline(points, tangents, time_points)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        # Test serialization
        a.saveAsText("serialization_curve.test")
        b = cubic_hermite_spline()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")
        a_pickled = pickle.dumps(a)
        a_from_pickle = pickle.loads(a_pickled)
        self.assertEqual(a_from_pickle, a)
        return

    def test_piecewise_polynomial_curve(self):
        print("test_piecewise_polynomial_curve")
        # To test :
        # - Functions : constructor, min, max, derivate, add_curve,
        #               is_continuous, serialize, deserialize
        waypoints0 = array([[0.0, 0.0, 0.0]]).transpose()
        waypoints1 = array([[1.0, 1.0, 1.0]]).transpose()
        waypoints2 = array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]).transpose()
        polynomial(waypoints0, 0.0, 0.1)
        a = polynomial(waypoints1, 0.0, 1.0)
        b = polynomial(waypoints2, 1.0, 3.0)
        pc = piecewise(a)
        pc.append(b)
        pc.min()
        pc.max()
        pc(0.4)
        self.assertTrue((pc(pc.min()) == array([[1.0, 1.0, 1.0]]).transpose()).all())
        self.assertTrue((pc.derivate(0.4, 0) == pc(0.4)).all())
        pc.derivate(0.4, 2)
        pc.is_continuous(0)
        pc.is_continuous(1)
        a0 = pc.curve_at_index(0)
        self.assertTrue(array_equal(a0(0.5), pc(0.5)))
        self.assertTrue(a0 == a)
        a0 = b  # should not have any effect
        self.assertTrue(array_equal(pc.curve_at_index(0)(0.5), a(0.5)))
        # Test serialization
        pc.saveAsText("serialization_pc.test")
        pc_test = piecewise()
        pc_test.loadFromText("serialization_pc.test")
        self.assertTrue((pc(0.4) == pc_test(0.4)).all())
        os.remove("serialization_pc.test")
        pc_pickled = pickle.dumps(pc)
        pc_from_pickle = pickle.loads(pc_pickled)
        self.assertEqual(pc_from_pickle, pc)

        waypoints3 = array(
            [[1.0, 2.0, 3.0, 0.6, -9.0], [-1.0, 1.6, 1.7, 6.7, 14]]
        ).transpose()
        c = polynomial(waypoints3, 3.0, 5.2)
        with self.assertRaises(ValueError):  # a and c doesn't have the same dimension
            pc.append(c)

        # Test the different append methods :
        pc = piecewise()
        self.assertEqual(pc.num_curves(), 0)
        end_point1 = array([1.0, 3.0, 5.0, 6.5, -2.0])
        max1 = 2.5
        with self.assertRaises(
            RuntimeError
        ):  # cannot add final point in an empty curve
            pc.append(end_point1, max1)
        with self.assertRaises(
            ValueError
        ):  # a and end_point1 doesn't have the same dimension
            pc.append(a)
            pc.append(end_point1, max1)

        pc = piecewise()
        d = polynomial(waypoints3, 0.0, 1.2)
        self.assertEqual(pc.num_curves(), 0)
        pc.append(d)
        self.assertEqual(pc.num_curves(), 1)
        pc.append(end_point1, max1)
        self.assertEqual(pc.num_curves(), 2)
        self.assertEqual(pc.min(), 0.0)
        self.assertEqual(pc.max(), max1)
        self.assertTrue(pc.is_continuous(0))
        self.assertTrue(isclose(pc(0.0), d(0.0)).all())
        self.assertTrue(isclose(pc(max1), end_point1).all())

        return

    def test_piecewise_from_points_list(self):
        N = 7
        rng = random.default_rng()
        points = array(rng.random((3, N)))
        points_derivative = array(rng.random((3, N)))
        points_second_derivative = array(rng.random((3, N)))
        time_points = array(rng.random((1, N))).T
        time_points.sort(0)
        polC0 = piecewise.FromPointsList(points, time_points)
        self.assertEqual(polC0.min(), time_points[0, 0])
        self.assertEqual(polC0.max(), time_points[-1, 0])
        self.assertTrue(polC0.is_continuous(0))
        self.assertTrue(not polC0.is_continuous(1))
        for i in range(N):
            self.assertTrue(isclose(polC0(time_points[i, 0]), points[:, i]).all())

        c0 = polC0.curve_at_index(0)
        self.assertEqual(c0.min(), time_points[0])
        self.assertEqual(c0.max(), time_points[1])
        self.assertEqual(c0.dim(), 3)
        mid_t = (c0.max() + c0.min()) / 2.0
        self.assertTrue(array_equal(polC0(mid_t), c0(mid_t)))

        polC1 = piecewise.FromPointsList(points, points_derivative, time_points)
        self.assertEqual(polC1.min(), time_points[0, 0])
        self.assertEqual(polC1.max(), time_points[-1, 0])
        self.assertTrue(polC1.is_continuous(0))
        self.assertTrue(polC1.is_continuous(1))
        self.assertTrue(not polC1.is_continuous(2))
        for i in range(N):
            self.assertTrue(isclose(polC1(time_points[i, 0]), points[:, i]).all())
            self.assertTrue(
                isclose(
                    polC1.derivate(time_points[i, 0], 1), points_derivative[:, i]
                ).all()
            )

        polC2 = piecewise.FromPointsList(
            points, points_derivative, points_second_derivative, time_points
        )
        self.assertEqual(polC2.min(), time_points[0, 0])
        self.assertEqual(polC2.max(), time_points[-1, 0])
        self.assertTrue(polC2.is_continuous(0))
        self.assertTrue(polC2.is_continuous(1))
        self.assertTrue(polC2.is_continuous(2))
        self.assertTrue(not polC2.is_continuous(3))
        for i in range(N):
            self.assertTrue(isclose(polC2(time_points[i, 0]), points[:, i]).all())
            self.assertTrue(
                isclose(
                    polC2.derivate(time_points[i, 0], 1), points_derivative[:, i]
                ).all()
            )
            self.assertTrue(
                isclose(
                    polC2.derivate(time_points[i, 0], 2), points_second_derivative[:, i]
                ).all()
            )

        # check if exepetion are corectly raised
        # when time_points are not in ascending values
        time_points[0, 0] = 1
        time_points[1, 0] = 0.5
        with self.assertRaises(ValueError):
            polC0 = piecewise.FromPointsList(points, time_points)

        with self.assertRaises(ValueError):
            polC1 = piecewise.FromPointsList(points, points_derivative, time_points)

        with self.assertRaises(ValueError):
            polC2 = piecewise.FromPointsList(
                points, points_derivative, points_second_derivative, time_points
            )

    def test_piecewise_bezier_curve(self):
        # To test :
        # - Functions : constructor, min, max, derivate, add_curve, is_continuous
        waypoints = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        a = bezier(waypoints, 0.0, 1.0)
        b = bezier(waypoints, 1.0, 2.0)
        pc = piecewise(a)
        pc.append(b)
        a.split(array([0.4, 0.8])).curve_at_index(0)
        pc.min()
        pc.max()
        pc(0.4)
        self.assertTrue((pc(pc.min()) == array([1.0, 2.0, 3.0])).all())
        self.assertTrue((pc.derivate(0.4, 0) == pc(0.4)).all())
        pc.derivate(0.4, 2)
        pc.is_continuous(0)
        pc.is_continuous(1)
        # test access to curves :
        self.assertTrue(array_equal(pc.curve_at_index(0)(0.5), a(0.5)))
        waypoints = array([[3.0, 4.0, -3.0], [5.0, 1.0, 2.0]]).transpose()
        c = bezier(waypoints, 1.5, 2.0)
        c0 = pc.curve_at_index(0)
        c0 = c  # should not have any effect
        c0
        self.assertTrue(array_equal(pc.curve_at_index(0)(0.5), a(0.5)))

        # Test serialization
        pc.saveAsText("serialization_pc.test")
        pc_test = piecewise()
        pc_test.loadFromText("serialization_pc.test")
        self.assertTrue((pc(0.4) == pc_test(0.4)).all())
        os.remove("serialization_pc.test")
        pc_pickled = pickle.dumps(pc)
        pc_from_pickle = pickle.loads(pc_pickled)
        self.assertEqual(pc_from_pickle, pc)
        return

    def test_piecewise_cubic_hermite_curve(self):
        print("test_piecewise_cubic_hermite_curve")
        # To test :
        # - Functions : constructor, min, max, derivate, add_curve, is_continuous
        points = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        tangents = array([[2.0, 2.0, 2.0], [4.0, 4.0, 4.0]]).transpose()
        time_points0 = array([[0.0, 1.0]]).transpose()
        time_points1 = array([[1.0, 2.0]]).transpose()
        a = cubic_hermite_spline(points, tangents, time_points0)
        b = cubic_hermite_spline(points, tangents, time_points1)
        pc = piecewise(a)
        pc.append(b)
        pc.min()
        pc.max()
        pc(0.4)
        self.assertTrue((pc(0.0) == array([1.0, 2.0, 3.0])).all())
        self.assertTrue((pc.derivate(0.4, 0) == pc(0.4)).all())
        pc.derivate(0.4, 2)
        pc.is_continuous(0)
        pc.is_continuous(1)
        a0 = pc.curve_at_index(0)
        self.assertTrue(array_equal(a0(0.5), pc(0.5)))
        self.assertTrue(a0 == a)
        a0 = b  # should not have any effect
        self.assertTrue(array_equal(pc.curve_at_index(0)(0.5), a(0.5)))
        # Test serialization
        pc.saveAsText("serialization_pc.test")
        pc_test = piecewise()
        pc_test.loadFromText("serialization_pc.test")
        self.assertTrue((pc(0.4) == pc_test(0.4)).all())
        os.remove("serialization_pc.test")
        pc_pickled = pickle.dumps(pc)
        pc_from_pickle = pickle.loads(pc_pickled)
        self.assertEqual(pc_from_pickle, pc)
        return

    def test_exact_cubic(self):
        print("test_exact_cubic")
        # To test :
        # - Functions : constructor, min, max, derivate, getNumberSplines, getSplineAt
        waypoints = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        time_waypoints = array([[0.0, 1.0]]).transpose()
        a = exact_cubic(waypoints, time_waypoints)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        a.getNumberSplines()
        a.getSplineAt(0)
        # Test serialization
        a.saveAsText("serialization_pc.test")
        b = exact_cubic()
        b.loadFromText("serialization_pc.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_pc.test")
        a_pickled = pickle.dumps(a)
        a_from_pickle = pickle.loads(a_pickled)
        self.assertEqual(a_from_pickle, a)
        return

    def test_exact_cubic_constraint(self):
        print("test_exact_cubic_constraint")
        # To test :
        # - Functions : constructor, min, max, derivate
        waypoints = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        time_waypoints = array([[0.0, 1.0]]).transpose()
        c = curve_constraints(3)
        c.init_vel = array([[0.0, 1.0, 1.0]]).transpose()
        c.end_vel = array([[0.0, 1.0, 1.0]]).transpose()
        c.init_acc = array([[0.0, 1.0, 1.0]]).transpose()
        c.end_acc = array([[0.0, 1.0, 1.0]]).transpose()
        c.init_vel
        c.end_vel
        c.init_acc
        c.end_acc
        exact_cubic(waypoints, time_waypoints)
        exact_cubic(waypoints, time_waypoints, c)
        return

    def test_cubic_hermite_spline_2(self):
        points = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        tangents = array([[2.0, 2.0, 2.0], [4.0, 4.0, 4.0]]).transpose()
        time_points = array([[0.0, 1.0]]).transpose()
        a = cubic_hermite_spline(points, tangents, time_points)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a(0.0) == array([1.0, 2.0, 3.0])).all())
        self.assertTrue(
            isclose(a.derivate(0.0, 1), array([[2.0, 2.0, 2.0]]).transpose()).all()
        )
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        return

    def test_conversion_curves(self):
        print("test_conversion_curves")
        __EPS = 1e-6
        waypoints = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        a = bezier(waypoints)
        # converting bezier to polynomial
        a_pol = convert_to_polynomial(a)
        self.assertTrue(norm(a(0.3) - a_pol(0.3)) < __EPS)
        # converting polynomial to hermite
        a_chs = convert_to_hermite(a_pol)
        self.assertTrue(norm(a_chs(0.3) - a_pol(0.3)) < __EPS)
        # converting hermite to bezier
        a_bc = convert_to_bezier(a_chs)
        self.assertTrue(norm(a_chs(0.3) - a_bc(0.3)) < __EPS)
        self.assertTrue(norm(a(0.3) - a_bc(0.3)) < __EPS)
        # converting bezier to hermite
        a_chs = convert_to_hermite(a)
        self.assertTrue(norm(a(0.3) - a_chs(0.3)) < __EPS)
        # converting hermite to polynomial
        a_pol = convert_to_polynomial(a_chs)
        self.assertTrue(norm(a_pol(0.3) - a_chs(0.3)) < __EPS)
        # converting polynomial to bezier
        a_bc = convert_to_bezier(a_pol)
        self.assertTrue(norm(a_bc(0.3) - a_pol(0.3)) < __EPS)
        self.assertTrue(norm(a(0.3) - a_bc(0.3)) < __EPS)
        return

    def test_piecewise_se3_curve(self):
        init_quat = Quaternion.Identity()
        end_quat = Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0, 0)
        init_rot = init_quat.matrix()
        end_rot = end_quat.matrix()
        waypoints = array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
            ]
        ).transpose()
        min = 0.2
        max = 1.5
        translation = bezier3(waypoints, min, max)
        # test with bezier
        se3 = SE3Curve(translation, init_rot, end_rot)
        pc = piecewise_SE3(se3)
        self.assertEqual(pc.num_curves(), 1)
        self.assertEqual(pc.min(), min)
        self.assertEqual(pc.max(), max)
        pmin = pc(min)
        pmax = pc(max)
        self.assertTrue(isclose(pmin[:3, :3], init_rot).all())
        self.assertTrue(isclose(pmax[:3, :3], end_rot).all())
        self.assertTrue(isclose(pmin[0:3, 3], translation(min)).all())
        self.assertTrue(isclose(pmax[0:3, 3], translation(max)).all())
        # add another curve :
        end_pos2 = array([-2, 0.2, 1.6])
        max2 = 2.7
        se3_2 = SE3Curve(translation(max), end_pos2, end_rot, end_rot, max, max2)
        pc.append(se3_2)
        self.assertEqual(pc.num_curves(), 2)
        pmin2 = pc(max)
        pmax2 = pc(max2)
        self.assertTrue(isclose(pmin2[:3, :3], end_rot).all())
        self.assertTrue(isclose(pmax2[:3, :3], end_rot).all())
        self.assertTrue(isclose(pmin2[0:3, 3], se3_2.translation(max)).all())
        self.assertTrue(isclose(pmax2[0:3, 3], se3_2.translation(max2)).all())
        self.assertTrue(pc.is_continuous(0))
        self.assertFalse(pc.is_continuous(1))

        # check if error are correctly raised :
        with self.assertRaises(ValueError):  # time intervals do not match
            se3_3 = SE3Curve(se3_2(max2), se3_2(max2 - 0.5), max2 - 0.5, max2 + 1.5)
            pc.append(se3_3)
        with self.assertRaises(ValueError):
            se3_3 = SE3Curve(se3_2(max2), se3_2(max2 - 0.5), max2 + 0.1, max2 + 1.5)
            pc.append(se3_3)

        pc.saveAsText("serialization_curve.txt")
        pc_txt = piecewise_SE3()
        pc_txt.loadFromText("serialization_curve.txt")
        self.compareCurves(pc, pc_txt)

        pc.saveAsXML("serialization_curve.xml", "pc")
        pc_xml = piecewise_SE3()
        pc_xml.loadFromXML("serialization_curve.xml", "pc")
        self.compareCurves(pc, pc_xml)

        pc.saveAsBinary("serialization_curve")
        pc_bin = piecewise_SE3()
        pc_bin.loadFromBinary("serialization_curve")
        self.compareCurves(pc, pc_bin)

        pc_pickled = pickle.dumps(pc)
        pc_from_pickle = pickle.loads(pc_pickled)
        self.assertEqual(pc_from_pickle, pc)

        se3_3 = SE3Curve(se3(max), se3_2(max2 - 0.5), max2, max2 + 1.5)
        pc.append(se3_3)
        self.assertFalse(pc.is_continuous(0))

        #  test the different append methods :
        init_quat = Quaternion.Identity()
        end_quat = Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0, 0)
        init_rot = init_quat.matrix()
        end_rot = end_quat.matrix()
        waypoints = array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
            ]
        ).transpose()
        min = 0.2
        max = 1.5
        translation = bezier3(waypoints, min, max)
        # test with bezier
        se3 = SE3Curve(translation, init_rot, end_rot)
        pc = piecewise_SE3()
        self.assertEqual(pc.num_curves(), 0)
        pc.append(se3)
        self.assertEqual(pc.num_curves(), 1)
        self.assertEqual(pc.min(), min)
        self.assertEqual(pc.max(), max)
        pmin = pc(min)
        pmax = pc(max)
        self.assertTrue(isclose(pmin[:3, :3], init_rot).all())
        self.assertTrue(isclose(pmax[:3, :3], end_rot).all())
        self.assertTrue(isclose(pmin[0:3, 3], translation(min)).all())
        self.assertTrue(isclose(pmax[0:3, 3], translation(max)).all())
        # append a final tranform :
        end_quat = Quaternion(sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0)
        end_rot = end_quat.matrix()
        end_translation = array([1.7, -0.8, 3.0]).T
        end_pose = array(np.identity(4))
        end_pose[:3, :3] = end_rot
        end_pose[:3, 3] = end_translation
        max2 = 3.8
        pc.append(end_pose, max2)
        self.assertEqual(pc.num_curves(), 2)
        self.assertEqual(pc.min(), min)
        self.assertEqual(pc.max(), max2)
        pmin = pc(min)
        pmax = pc(max2)
        self.assertTrue(isclose(pmin[:3, :3], init_rot).all())
        self.assertTrue(isclose(pmax[:3, :3], end_rot).all())
        self.assertTrue(isclose(pmin[0:3, 3], translation(min)).all())
        self.assertTrue(isclose(pmax[0:3, 3], end_translation).all())
        self.assertTrue(pc.is_continuous(0))
        if CURVES_WITH_PINOCCHIO_SUPPORT:
            end_quat = Quaternion(sqrt(2.0) / 2.0, 0.0, 0, sqrt(2.0) / 2.0)
            end_rot = end_quat.matrix()
            end_translation = array([-17.0, 3.7, 1.0])
            end_pose = SE3.Identity()
            end_pose.rotation = end_rot
            end_pose.translation = end_translation
            max3 = 6.5
            pc.append(end_pose, max3)
            self.assertEqual(pc.num_curves(), 3)
            self.assertEqual(pc.min(), min)
            self.assertEqual(pc.max(), max3)
            pmin = pc(min)
            pmax = pc(max3)
            self.assertTrue(isclose(pmin[:3, :3], init_rot).all())
            self.assertTrue(isclose(pmax[:3, :3], end_rot).all())
            self.assertTrue(isclose(pmin[0:3, 3], translation(min)).all())
            self.assertTrue(isclose(pmax[0:3, 3], end_translation).all())
            self.assertTrue(pc.is_continuous(0))
        pc = piecewise_SE3()
        with self.assertRaises(RuntimeError):
            pc.append(end_pose, max)

    if CURVES_WITH_PINOCCHIO_SUPPORT:

        def test_piecewise_se3_curve_linear_pinocchio(self):
            print("test piecewise SE3 pinocchio")
            init_quat = Quaternion.Identity()
            end_quat = Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0, 0)
            init_rot = init_quat.matrix()
            end_rot = end_quat.matrix()
            init_translation = array([0.2, -0.7, 0.6])
            end_translation = array([3.6, -2.2, -0.9])
            init_pose = SE3.Identity()
            end_pose = SE3.Identity()
            init_pose.rotation = init_rot
            end_pose.rotation = end_rot
            init_pose.translation = init_translation
            end_pose.translation = end_translation
            min = 0.7
            max = 12.0
            se3 = SE3Curve(init_pose, end_pose, min, max)
            pc = piecewise_SE3(se3)
            self.assertEqual(pc.num_curves(), 1)
            p = pc.evaluateAsSE3(min)
            self.assertTrue(isinstance(p, SE3))
            self.assertTrue(pc.evaluateAsSE3(min).isApprox(init_pose, 1e-6))
            self.assertTrue(pc.evaluateAsSE3(max).isApprox(end_pose, 1e-6))
            self.assertEqual(pc.min(), min)
            self.assertEqual(pc.max(), max)
            self.assertTrue(isclose(pc.rotation(min), init_rot).all())
            self.assertTrue(isclose(pc.translation(min), init_translation).all())
            self.assertTrue(isclose(pc.rotation(max), end_rot).all())
            self.assertTrue(isclose(pc.translation(max), end_translation).all())
            # add another curve :
            end_translation2 = array([-2.0, 1.6, -14.0])
            end_pose2 = SE3.Identity()
            end_pose2.rotation = end_rot
            end_pose2.translation = end_translation2
            max2 = 23.9
            se3_2 = SE3Curve(end_pose, end_pose2, max, max2)
            pc.append(se3_2)
            self.assertEqual(pc.num_curves(), 2)
            p = pc.evaluateAsSE3(max2)
            self.assertTrue(isinstance(p, SE3))
            self.assertTrue(pc.evaluateAsSE3(max2).isApprox(end_pose2, 1e-6))
            self.assertEqual(pc.min(), min)
            self.assertEqual(pc.max(), max2)
            self.assertTrue(isclose(pc.rotation(max2), end_rot).all())
            self.assertTrue(isclose(pc.translation(max2), end_translation2).all())

    def test_conversion_piecewise_curves(self):
        waypoints = array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        a = bezier(waypoints, 0.0, 1.0)
        b = bezier(waypoints, 1.0, 2.0)
        pc = piecewise(a)
        pc.append(b)
        # Convert to piecewise polynomial
        pc_pol = pc.convert_piecewise_curve_to_polynomial()
        self.compareCurves(pc_pol, pc)
        # Convert to piecewise cubic hermite
        pc_chs = pc.convert_piecewise_curve_to_cubic_hermite()
        self.compareCurves(pc_chs, pc)
        # Convert to piecewise bezier
        pc_bc = pc_chs.convert_piecewise_curve_to_bezier()
        self.compareCurves(pc_bc, pc)
        return

    def test_so3_linear(self):
        print("test SO3 Linear")
        init_quat = Quaternion.Identity()
        end_quat = Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0, 0)
        init_rot = init_quat.matrix()
        end_rot = end_quat.matrix()
        min = 0.2
        max = 1.5

        so3Rot = SO3Linear(init_rot, end_rot, min, max)
        so3Quat = SO3Linear(init_quat, end_quat, min, max)
        self.assertEqual(so3Rot.min(), min)
        self.assertEqual(so3Rot.max(), max)
        self.assertEqual(so3Quat.min(), min)
        self.assertEqual(so3Quat.max(), max)
        self.assertTrue(isclose(so3Rot(min), init_rot).all())
        self.assertTrue(isclose(so3Rot(max), end_rot).all())
        self.assertTrue(isclose(so3Quat(min), init_rot).all())
        self.assertTrue(isclose(so3Quat(max), end_rot).all())
        self.assertTrue(so3Rot.computeAsQuaternion(min).isApprox(init_quat))
        self.assertTrue(so3Rot.computeAsQuaternion(max).isApprox(end_quat))
        self.assertTrue(so3Quat.computeAsQuaternion(min).isApprox(init_quat))
        self.assertTrue(so3Quat.computeAsQuaternion(max).isApprox(end_quat))
        t = min
        while t < max:
            self.assertTrue(isclose(so3Quat(t), so3Rot(t)).all())
            t += 0.01
        # check the derivatives :
        vel = array([1.20830487, 0.0, 0.0])
        zeros3 = zeros((3, 1))
        t = min
        while t < max:
            self.assertTrue(isclose(so3Quat.derivate(t, 1), vel).all())
            self.assertTrue(isclose(so3Rot.derivate(t, 1), vel).all())
            t += 0.01
        for i in range(2, 5):
            t = min
            while t < max:
                self.assertTrue(isclose(so3Quat.derivate(t, i), zeros3).all())
                self.assertTrue(isclose(so3Rot.derivate(t, i), zeros3).all())
                t += 0.01

        # check that errors are correctly raised when necessary :
        with self.assertRaises(ValueError):
            so3Rot(0.0)
        with self.assertRaises(ValueError):
            so3Rot(-0.1)
        with self.assertRaises(ValueError):
            so3Rot(3)
        # TODO: these are not passing
        with self.assertRaises(ValueError):
            so3Rot.derivate(0, 1)
        with self.assertRaises(ValueError):
            so3Rot.derivate(3.0, 1)
        with self.assertRaises(ValueError):
            so3Rot.derivate(1.0, 0)
        with self.assertRaises(ValueError):
            SO3Linear(init_rot, end_rot, max, min)

    def test_so3_linear_serialization(self):
        print("test SO3 Linear")
        init_quat = Quaternion(0.590, -0.002, -0.766, 0.255)
        end_quat = Quaternion(-0.820, 0.162, 0.381, 0.396)
        init_quat.normalize()
        end_quat.normalize()
        init_rot = init_quat.matrix()
        end_rot = end_quat.matrix()
        min = 0.2
        max = 1.5
        so3Rot = SO3Linear(init_rot, end_rot, min, max)
        so3Rot.saveAsText("serialization_curve.txt")
        so3_txt = SO3Linear()
        so3_txt.loadFromText("serialization_curve.txt")
        self.compareCurves(so3Rot, so3_txt)

        so3Rot.saveAsXML("serialization_curve.xml", "so3")
        so3_xml = SO3Linear()
        so3_xml.loadFromXML("serialization_curve.xml", "so3")
        self.compareCurves(so3Rot, so3_xml)

        so3Rot.saveAsBinary("serialization_curve")
        so3_bin = SO3Linear()
        so3_bin.loadFromBinary("serialization_curve")
        self.compareCurves(so3Rot, so3_bin)

        so3Rot_pickled = pickle.dumps(so3Rot)
        so3Rot_from_pickle = pickle.loads(so3Rot_pickled)
        self.assertEqual(so3Rot_from_pickle, so3Rot)

    def test_se3_curve_linear(self):
        print("test SE3 Linear")
        init_quat = Quaternion.Identity()
        end_quat = Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0, 0)
        init_rot = init_quat.matrix()
        end_rot = end_quat.matrix()
        init_translation = array([1, 1.2, -0.6]).T
        end_translation = array([2.3, 0, 0.9]).T
        init_pose = array(np.identity(4))
        end_pose = array(np.identity(4))
        init_pose[:3, :3] = init_rot
        end_pose[:3, :3] = end_rot
        init_pose[:3, 3] = init_translation
        end_pose[:3, 3] = end_translation
        min = 0.2
        max = 1.5
        se3 = SE3Curve(init_pose, end_pose, min, max)
        p = se3(min)
        if CURVES_WITH_PINOCCHIO_SUPPORT:
            self.assertTrue(isinstance(se3.evaluateAsSE3(min), SE3))
            self.assertTrue(se3.evaluateAsSE3(min).isApprox(SE3(init_pose), 1e-6))
            self.assertTrue(se3.evaluateAsSE3(max).isApprox(SE3(end_pose), 1e-6))
        self.assertEqual(p.shape[0], 4)
        self.assertEqual(p.shape[1], 4)
        self.assertTrue(isclose(se3(min), init_pose).all())
        self.assertTrue(isclose(se3(max), end_pose).all())
        self.assertEqual(se3.min(), min)
        self.assertEqual(se3.max(), max)
        self.assertTrue(isclose(se3.rotation(min), init_rot).all())
        self.assertTrue(isclose(se3.translation(min), init_translation).all())
        self.assertTrue(isclose(se3.rotation(max), end_rot).all())
        self.assertTrue(isclose(se3.translation(max), end_translation).all())
        # check value of derivative (should be constant here)
        d = se3.derivate(min, 1)
        if CURVES_WITH_PINOCCHIO_SUPPORT:
            motion = se3.derivateAsMotion(min, 1)
            self.assertTrue(isinstance(motion, Motion))
            self.assertTrue(
                isclose(
                    motion.linear, ((end_translation - init_translation) / (max - min))
                ).all()
            )
            self.assertTrue(isclose(motion.angular[0], 1.20830487))
            self.assertTrue(isclose(motion.angular[1:3], array([0, 0]).T).all())
            self.assertTrue(motion.isApprox(se3.derivateAsMotion(0.5, 1), 1e-6))
            self.assertTrue(motion.isApprox(se3.derivateAsMotion(max, 1), 1e-6))
            self.assertTrue(se3.derivateAsMotion(min, 2).isApprox(Motion.Zero(), 1e-6))
            self.assertTrue(se3.derivateAsMotion(min, 3).isApprox(Motion.Zero(), 1e-6))
        self.assertEqual(d.shape[0], 6)
        # self.assertEqual(d.shape[1], 1)
        self.assertTrue(
            isclose(d[0:3], ((end_translation - init_translation) / (max - min))).all()
        )
        self.assertTrue(isclose(d[3], 1.20830487))
        self.assertTrue(isclose(d[4:6], array([0, 0]).T).all())
        self.assertTrue(isclose(d, se3.derivate(0.5, 1)).all())
        self.assertTrue(isclose(d, se3.derivate(max, 1)).all())
        self.assertTrue(isclose(se3.derivate(min, 2), zeros((6, 1))).all())
        self.assertTrue(isclose(se3.derivate(min, 3), zeros((6, 1))).all())

        # test accessor to translation_curve :
        tr_se3 = se3.translation_curve()
        self.assertTrue(
            array_equal(tr_se3((max + min) / 2.0), se3.translation((max + min) / 2.0))
        )
        # test accessor to rotation :
        rot_se3 = se3.rotation_curve()
        rot_se3((max + min) / 2.0)
        self.assertTrue(
            isclose(rot_se3((max + min) / 2.0), (se3.rotation((max + min) / 2.0))).all()
        )
        # check that it return a CONST reference :
        waypoints2 = array(
            [[1.0, -2.0, 3.5], [5.6, 5.0, -6.0], [4.0, 1.2, 0.5]]
        ).transpose()
        tr_se3 = bezier3(waypoints2, min, max)  # should not have any effect
        self.assertFalse(
            array_equal(tr_se3((max + min) / 2.0), se3.translation((max + min) / 2.0))
        )
        # check that errors are correctly raised when necessary :
        with self.assertRaises(ValueError):
            se3(0.0)
        with self.assertRaises(ValueError):
            se3(-0.1)
        with self.assertRaises(ValueError):
            se3(3)
        with self.assertRaises(ValueError):
            se3.derivate(0, 1)
        with self.assertRaises(ValueError):
            se3.derivate(3.0, 1)
        with self.assertRaises(ValueError):
            se3.derivate(1.0, 0)
        with self.assertRaises(ValueError):
            SE3Curve(init_pose, end_pose, max, min)

    def test_se3_from_translation_curve(self):
        print("test SE3 From translation curves")
        init_quat = Quaternion.Identity()
        end_quat = Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0, 0)
        init_rot = init_quat.matrix()
        end_rot = end_quat.matrix()
        waypoints = array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
            ]
        ).transpose()
        min = 0.2
        max = 1.5
        translation = bezier3(waypoints, min, max)
        # test with bezier
        se3 = SE3Curve(translation, init_rot, end_rot)
        self.assertEqual(se3.min(), min)
        self.assertEqual(se3.max(), max)
        pmin = se3(min)
        pmax = se3(max)
        self.assertTrue(isclose(pmin[:3, :3], init_rot).all())
        self.assertTrue(isclose(pmax[:3, :3], end_rot).all())
        self.assertTrue(isclose(pmin[0:3, 3], translation(min)).all())
        self.assertTrue(isclose(pmax[0:3, 3], translation(max)).all())
        if CURVES_WITH_PINOCCHIO_SUPPORT:
            pminSE3 = se3.evaluateAsSE3(min)
            pmaxSE3 = se3.evaluateAsSE3(max)
            self.assertTrue(isclose(pminSE3.rotation, init_rot).all())
            self.assertTrue(isclose(pmaxSE3.rotation, end_rot).all())
            self.assertTrue(isclose(pminSE3.translation, translation(min)).all())
            self.assertTrue(isclose(pmaxSE3.translation, translation(max)).all())
        t = min
        while t < max:
            if CURVES_WITH_PINOCCHIO_SUPPORT:
                self.assertTrue(
                    isclose(se3.evaluateAsSE3(t).translation, translation(t)).all()
                )
                self.assertTrue(
                    isclose(
                        se3.derivateAsMotion(t, 1).linear, translation.derivate(t, 1)
                    ).all()
                )
            self.assertTrue(isclose(se3(t)[0:3, 3], translation(t)).all())
            self.assertTrue(
                isclose(se3.derivate(t, 1)[0:3], translation.derivate(t, 1)).all()
            )
            t += 0.02

        # test accessor to translation_curve :
        tr_se3 = se3.translation_curve()
        self.assertTrue(tr_se3 == translation)
        self.assertTrue(
            array_equal(tr_se3((max + min) / 2.0), se3.translation((max + min) / 2.0))
        )
        # test accessor to rotation :
        rot_se3 = se3.rotation_curve()
        rot_se3((max + min) / 2.0)
        se3.rotation((max + min) / 2.0)
        self.assertTrue(
            isclose(rot_se3((max + min) / 2.0), se3.rotation((max + min) / 2.0)).all()
        )
        # check that it return a CONST reference :
        waypoints2 = array(
            [[1.0, -2.0, 3.5], [5.6, 5.0, -6.0], [4.0, 1.2, 0.5]]
        ).transpose()
        tr_se3 = bezier3(waypoints2, min, max)  # should not have any effect
        self.assertFalse(
            array_equal(tr_se3((max + min) / 2.0), se3.translation((max + min) / 2.0))
        )

        # test with bezier3
        translation = bezier3(waypoints, min, max)
        se3 = SE3Curve(translation, init_rot, end_rot)
        self.assertEqual(se3.min(), min)
        self.assertEqual(se3.max(), max)
        pmin = se3(min)
        pmax = se3(max)
        self.assertTrue(isclose(pmin[:3, :3], init_rot).all())
        self.assertTrue(isclose(pmax[:3, :3], end_rot).all())
        self.assertTrue(isclose(pmin[0:3, 3], translation(min)).all())
        self.assertTrue(isclose(pmax[0:3, 3], translation(max)).all())
        if CURVES_WITH_PINOCCHIO_SUPPORT:
            pminSE3 = se3.evaluateAsSE3(min)
            pmaxSE3 = se3.evaluateAsSE3(max)
            self.assertTrue(isclose(pminSE3.rotation, init_rot).all())
            self.assertTrue(isclose(pmaxSE3.rotation, end_rot).all())
            self.assertTrue(isclose(pminSE3.translation, translation(min)).all())
            self.assertTrue(isclose(pmaxSE3.translation, translation(max)).all())
        t = min
        while t < max:
            if CURVES_WITH_PINOCCHIO_SUPPORT:
                self.assertTrue(
                    isclose(se3.evaluateAsSE3(t).translation, translation(t)).all()
                )
                self.assertTrue(
                    isclose(
                        se3.derivateAsMotion(t, 1).linear, translation.derivate(t, 1)
                    ).all()
                )
            self.assertTrue(isclose(se3(t)[0:3, 3], translation(t)).all())
            self.assertTrue(
                isclose(se3.derivate(t, 1)[0:3], translation.derivate(t, 1)).all()
            )
            t += 0.02

        # test with piecewise polynomial
        N = 7
        rng = random.default_rng()
        points = array(rng.random((3, N)))
        # points_derivative = array(random.rand(3, N))
        # points_second_derivative = array(random.rand(3, N))
        time_points = array(rng.random((1, N))).T
        time_points.sort(0)
        translation = piecewise3.FromPointsList(points, time_points)
        min = translation.min()
        max = translation.max()
        se3 = SE3Curve(translation, init_rot, end_rot)
        self.assertEqual(se3.min(), min)
        self.assertEqual(se3.max(), max)
        pmin = se3(min)
        pmax = se3(max)
        self.assertTrue(isclose(pmin[:3, :3], init_rot).all())
        self.assertTrue(isclose(pmax[:3, :3], end_rot).all())
        self.assertTrue(isclose(pmin[0:3, 3], translation(min)).all())
        self.assertTrue(isclose(pmax[0:3, 3], translation(max)).all())
        if CURVES_WITH_PINOCCHIO_SUPPORT:
            pminSE3 = se3.evaluateAsSE3(min)
            pmaxSE3 = se3.evaluateAsSE3(max)
            self.assertTrue(isclose(pminSE3.rotation, init_rot).all())
            self.assertTrue(isclose(pmaxSE3.rotation, end_rot).all())
            self.assertTrue(isclose(pminSE3.translation, translation(min)).all())
            self.assertTrue(isclose(pmaxSE3.translation, translation(max)).all())
        t = min
        while t < max:
            if CURVES_WITH_PINOCCHIO_SUPPORT:
                self.assertTrue(
                    isclose(se3.evaluateAsSE3(t).translation, translation(t)).all()
                )
                self.assertTrue(
                    isclose(
                        se3.derivateAsMotion(t, 1).linear, translation.derivate(t, 1)
                    ).all()
                )
            self.assertTrue(isclose(se3(t)[0:3, 3], translation(t)).all())
            self.assertTrue(
                isclose(se3.derivate(t, 1)[0:3], translation.derivate(t, 1)).all()
            )
            t += 0.02

    def test_se3_from_curves(self):
        print("test SE3 from curves")
        init_quat = Quaternion.Identity()
        end_quat = Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0, 0)
        init_rot = init_quat.matrix()
        end_rot = end_quat.matrix()
        waypoints = array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
            ]
        ).transpose()
        min = 0.2
        max = 1.5
        translation = bezier3(waypoints, min, max)
        rotation = SO3Linear(init_rot, end_rot, min, max)
        se3 = SE3Curve(translation, rotation)
        self.assertEqual(se3.min(), min)
        self.assertEqual(se3.max(), max)
        pmin = se3(min)
        pmax = se3(max)
        self.assertTrue(isclose(pmin[:3, :3], init_rot).all())
        self.assertTrue(isclose(pmax[:3, :3], end_rot).all())
        self.assertTrue(isclose(pmin[0:3, 3], translation(min)).all())
        self.assertTrue(isclose(pmax[0:3, 3], translation(max)).all())
        if CURVES_WITH_PINOCCHIO_SUPPORT:
            pminSE3 = se3.evaluateAsSE3(min)
            pmaxSE3 = se3.evaluateAsSE3(max)
            self.assertTrue(isclose(pminSE3.rotation, init_rot).all())
            self.assertTrue(isclose(pmaxSE3.rotation, end_rot).all())
            self.assertTrue(isclose(pminSE3.translation, translation(min)).all())
            self.assertTrue(isclose(pmaxSE3.translation, translation(max)).all())
        t = min
        while t < max:
            if CURVES_WITH_PINOCCHIO_SUPPORT:
                self.assertTrue(
                    isclose(se3.evaluateAsSE3(t).translation, translation(t)).all()
                )
                self.assertTrue(
                    isclose(se3.evaluateAsSE3(t).rotation, rotation(t)).all()
                )
                self.assertTrue(
                    isclose(
                        se3.derivateAsMotion(t, 1).linear, translation.derivate(t, 1)
                    ).all()
                )
                self.assertTrue(
                    isclose(
                        se3.derivateAsMotion(t, 1).angular, rotation.derivate(t, 1)
                    ).all()
                )
            self.assertTrue(isclose(se3(t)[0:3, 3], translation(t)).all())
            self.assertTrue(isclose(se3(t)[0:3, 0:3], rotation(t)).all())
            self.assertTrue(
                isclose(se3.derivate(t, 1)[0:3], translation.derivate(t, 1)).all()
            )
            self.assertTrue(
                isclose(se3.derivate(t, 1)[3:6], rotation.derivate(t, 1)).all()
            )
            t += 0.02

        # test accessor to translation_curve :
        tr_se3 = se3.translation_curve()
        self.assertTrue(tr_se3 == translation)
        self.assertTrue(
            array_equal(tr_se3((max + min) / 2.0), se3.translation((max + min) / 2.0))
        )
        # test accessor to rotation :
        rot_se3 = se3.rotation_curve()
        rot_se3((max + min) / 2.0)
        se3.rotation((max + min) / 2.0)
        self.assertTrue(
            array_equal(rot_se3((max + min) / 2.0), se3.rotation((max + min) / 2.0))
        )
        # check that it return a CONST reference :
        waypoints2 = array(
            [[1.0, -2.0, 3.5], [5.6, 5.0, -6.0], [4.0, 1.2, 0.5]]
        ).transpose()
        tr_se3 = bezier3(waypoints2, min, max)  # should not have any effect
        self.assertFalse(
            array_equal(tr_se3((max + min) / 2.0), se3.translation((max + min) / 2.0))
        )

        # check if errors are correctly raised :
        rotation = SO3Linear(init_rot, end_rot, min + 0.2, max)
        with self.assertRaises(ValueError):
            se3 = SE3Curve(translation, rotation)
        rotation = SO3Linear(init_rot, end_rot, min - 0.1, max)
        with self.assertRaises(ValueError):
            se3 = SE3Curve(translation, rotation)
        rotation = SO3Linear(init_rot, end_rot, min, max + 0.5)
        with self.assertRaises(ValueError):
            se3 = SE3Curve(translation, rotation)
        rotation = SO3Linear(init_rot, end_rot, min, max - 0.1)
        with self.assertRaises(ValueError):
            se3 = SE3Curve(translation, rotation)

    if CURVES_WITH_PINOCCHIO_SUPPORT:

        def test_se3_curve_linear_pinocchio(self):
            print("test SE3 Linear pinocchio")
            init_quat = Quaternion.Identity()
            end_quat = Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0, 0)
            init_rot = init_quat.matrix()
            end_rot = end_quat.matrix()
            init_translation = array([0.2, -0.7, 0.6])
            end_translation = array([3.6, -2.2, -0.9])
            init_pose = SE3.Identity()
            end_pose = SE3.Identity()
            init_pose.rotation = init_rot
            end_pose.rotation = end_rot
            init_pose.translation = init_translation
            end_pose.translation = end_translation
            min = 0.7
            max = 12.0
            se3 = SE3Curve(init_pose, end_pose, min, max)
            p = se3.evaluateAsSE3(min)
            self.assertTrue(isinstance(p, SE3))
            self.assertTrue(se3.evaluateAsSE3(min).isApprox(init_pose, 1e-6))
            self.assertTrue(se3.evaluateAsSE3(max).isApprox(end_pose, 1e-6))
            self.assertEqual(se3.min(), min)
            self.assertEqual(se3.max(), max)
            self.assertTrue(isclose(se3.rotation(min), init_rot).all())
            self.assertTrue(isclose(se3.translation(min), init_translation).all())
            self.assertTrue(isclose(se3.rotation(max), end_rot).all())
            self.assertTrue(isclose(se3.translation(max), end_translation).all())
            # check value of derivative (should be constant here)
            d = se3.derivateAsMotion(min, 1)
            self.assertTrue(isinstance(d, Motion))
            self.assertTrue(
                isclose(
                    d.linear, ((end_translation - init_translation) / (max - min))
                ).all()
            )
            self.assertTrue(isclose(d.angular[0], 0.139009))
            self.assertTrue(isclose(d.angular[1:3], array([0, 0]).T).all())
            self.assertTrue(
                d.isApprox(se3.derivateAsMotion((min + max) / 2.0, 1), 1e-6)
            )
            self.assertTrue(d.isApprox(se3.derivateAsMotion(max, 1), 1e-6))
            self.assertTrue(se3.derivateAsMotion(min, 2).isApprox(Motion.Zero(), 1e-6))
            self.assertTrue(se3.derivateAsMotion(min, 3).isApprox(Motion.Zero(), 1e-6))
            # check that errors are correctly raised when necessary :
            with self.assertRaises(ValueError):
                se3(0.0)
            with self.assertRaises(ValueError):
                se3(-0.1)
            with self.assertRaises(ValueError):
                se3(13.0)
            with self.assertRaises(ValueError):
                se3.derivate(0, 1)
            with self.assertRaises(ValueError):
                se3.derivate(13.0, 1)
            with self.assertRaises(ValueError):
                se3.derivate(1.0, 0)

    def test_se3_serialization(self):
        print("test serialization SE3")
        init_quat = Quaternion(0.590, -0.002, -0.766, 0.255)
        end_quat = Quaternion(-0.820, 0.162, 0.381, 0.396)
        init_quat.normalize()
        end_quat.normalize()
        init_rot = init_quat.matrix()
        end_rot = end_quat.matrix()
        init_translation = array([1, 1.2, -0.6]).T
        end_translation = array([2.3, 0, 0.9]).T
        init_pose = array(np.identity(4))
        end_pose = array(np.identity(4))
        init_pose[:3, :3] = init_rot
        end_pose[:3, :3] = end_rot
        init_pose[:3, 3] = init_translation
        end_pose[:3, 3] = end_translation
        min = 0.2
        max = 1.5
        se3_linear = SE3Curve(init_pose, end_pose, min, max)

        se3_linear.saveAsText("serialization_curve.txt")
        se3_txt = SE3Curve()
        se3_txt.loadFromText("serialization_curve.txt")
        self.compareCurves(se3_linear, se3_txt)

        se3_linear.saveAsXML("serialization_curve.xml", "se3")
        se3_xml = SE3Curve()
        se3_xml.loadFromXML("serialization_curve.xml", "se3")
        self.compareCurves(se3_linear, se3_xml)

        se3_linear.saveAsBinary("serialization_curve")
        se3_bin = SE3Curve()
        se3_bin.loadFromBinary("serialization_curve")
        self.compareCurves(se3_linear, se3_bin)

        se3_pickled = pickle.dumps(se3_linear)
        se3_from_pickle = pickle.loads(se3_pickled)
        self.assertEqual(se3_from_pickle, se3_linear)

        # test from two curves :
        init_quat = Quaternion.Identity()
        end_quat = Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0, 0)
        init_rot = init_quat.matrix()
        end_rot = end_quat.matrix()
        waypoints = array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
            ]
        ).transpose()
        min = 0.2
        max = 1.5
        translation = bezier3(waypoints, min, max)
        rotation = SO3Linear(init_rot, end_rot, min, max)
        se3_curves = SE3Curve(translation, rotation)

        se3_curves.saveAsText("serialization_curve.txt")
        se3_txt = SE3Curve()
        se3_txt.loadFromText("serialization_curve.txt")
        self.compareCurves(se3_curves, se3_txt)

        se3_curves.saveAsXML("serialization_curve.xml", "se3")
        se3_xml = SE3Curve()
        se3_xml.loadFromXML("serialization_curve.xml", "se3")
        self.compareCurves(se3_curves, se3_xml)

        se3_curves.saveAsBinary("serialization_curve")
        se3_bin = SE3Curve()
        se3_bin.loadFromBinary("serialization_curve")
        self.compareCurves(se3_curves, se3_bin)

        se3_pickled = pickle.dumps(se3_curves)
        se3_from_pickle = pickle.loads(se3_pickled)
        self.assertEqual(se3_from_pickle, se3_curves)

    def test_operatorEqual(self):
        # test with bezier
        waypoints = array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
            ]
        ).transpose()
        a0 = bezier(waypoints)
        a1 = bezier(waypoints, 0.0, 3.0)
        a2 = bezier(waypoints, 0.0, 3.0)
        self.assertTrue(a0 != a1)
        self.assertTrue(a1 == a2)

        # test with polynomials of degree 5
        p0 = array([1.0, 3.0, -2.0])
        p1 = array([0.6, 2.0, 2.5])
        dp0 = array([-6.0, 2.0, -1.0])
        dp1 = array([10.0, 10.0, 10.0])
        ddp0 = array([1.0, -7.0, 4.5])
        ddp1 = array([6.0, -1.0, -4])
        min = 1.0
        max = 2.5
        pol_1 = polynomial(p0, dp0, ddp0, p1, dp1, ddp1, min, max)
        pol_2 = polynomial(p0, dp0, ddp0, p1, dp1, ddp1, min, max)
        pol_3 = polynomial(p1, dp0, ddp0, p0, dp1, ddp1, min, max)
        self.assertTrue(pol_1 == pol_2)
        self.assertTrue(pol_1 != pol_3)

        # test with polynomial/bezier
        pol_4 = polynomial(p0, dp0, p1, dp1, min, max)
        b_4 = convert_to_bezier(pol_4)
        self.assertTrue(pol_4.isEquivalent(b_4))
        self.assertTrue(pol_4.isEquivalent(b_4, 1e-6))
        self.assertTrue(pol_4.isEquivalent(b_4, 1e-6, 2))
        self.assertFalse(pol_4.isEquivalent(a1))

        # test with SE3:
        init_quat = Quaternion.Identity()
        end_quat = Quaternion(sqrt(2.0) / 2.0, sqrt(2.0) / 2.0, 0, 0)
        init_rot = init_quat.matrix()
        end_rot = end_quat.matrix()
        waypoints = array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
                [4.0, 5.0, 6.0],
            ]
        ).transpose()
        min = 0.2
        max = 1.5
        translation = bezier3(waypoints, min, max)
        se3_1 = SE3Curve(translation, init_rot, end_rot)
        se3_2 = SE3Curve(translation, init_rot, end_rot)
        waypoints2 = array(
            [
                [1.0, 2.0, 3.5],
                [4.0, 5.0, 6.0],
                [-4, 5.0, 6.0],
                [4.0, 8.0, 6.0],
                [4.0, 5.0, 6.0],
            ]
        ).transpose()
        translation3 = bezier3(waypoints2, min, max)
        se3_3 = SE3Curve(translation3, init_rot, end_rot)
        end_quat2 = Quaternion(sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0)
        end_rot2 = end_quat2.matrix()
        se3_4 = SE3Curve(translation, init_rot, end_rot2)
        self.assertTrue(se3_1 == se3_2)
        self.assertTrue(se3_1 != se3_3)
        self.assertTrue(se3_1 != se3_4)


if __name__ == "__main__":
    unittest.main()
