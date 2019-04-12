#!/usr/bin/env python

import unittest

from numpy import matrix
from numpy.linalg import norm
from numpy.testing import assert_allclose

from curves import bezier3, bezier6, curve_constraints, exact_cubic, from_bezier, polynom, spline_deriv_constraint


class TestCurve(unittest.TestCase):
    def test_curve(self):
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).T
        waypoints6 = matrix([[1., 2., 3., 7., 5., 5.], [4., 5., 6., 4., 5., 6.]]).T
        time_waypoints = matrix([0., 1.]).T

        # Test : Bezier3/Bezier6, polynom, exact_cubic, curve_constraints, spline_deriv_constraints, bernstein
        # Done : Bezier3/Bezier6, polynom, exact_cubic, curve_constraints, spline_deriv_constraints
        # To do ? Bernstein

        # TESTING BEZIER CURVE 
        # - Functions : constructor, min, max, derivate,compute_derivate, compute_primitive
        # - Variables : degree, nbWayPoints
        # Create bezier6 and bezier3
        a = bezier6(waypoints6)
        a = bezier3(waypoints, 3.)
        # Test : Degree, min, max, derivate
        self.assertEqual(a.degree, a.nbWaypoints - 1)
        a.min()
        a.max()
        a(0.4)
        assert_allclose(a.derivate(0.4, 0), a(0.4))
        a.derivate(0.4, 2)
        a = a.compute_derivate(100)
        prim = a.compute_primitive(1)
        # Check primitive and derivate - order 1
        for i in range(10):
            t = float(i) / 10.
            assert_allclose(a(t), prim.derivate(t, 1))
        assert_allclose(prim(0), matrix([0., 0., 0.]).T)
        # Check primitive and derivate - order 2
        prim = a.compute_primitive(2)
        for i in range(10):
            t = float(i) / 10.
            assert_allclose(a(t), prim.derivate(t, 2), atol=1e-20)
        assert_allclose(prim(0), matrix([0., 0., 0.]).T)
        # Create new bezier3 curve
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.], [4., 5., 6.], [4., 5., 6.], [4., 5., 6.]]).T
        a0 = bezier3(waypoints)
        a1 = bezier3(waypoints, 3.)
        prim0 = a0.compute_primitive(1)
        prim1 = a1.compute_primitive(1)
        # Check change in argument time_t of bezier3
        for i in range(10):
            t = float(i) / 10.
            assert_allclose(a0(t), a1(3 * t))
            assert_allclose(a0.derivate(t, 1), a1.derivate(3 * t, 1) * 3.)
            assert_allclose(a0.derivate(t, 2), a1.derivate(3 * t, 2) * 9.)
            assert_allclose(prim0(t), prim1(t * 3) / 3.)
        assert_allclose(prim(0), matrix([0., 0., 0.]).T)
        with self.assertRaises(AssertionError):
            assert_allclose(prim(0), matrix([0., 0., 0.]))
        # testing bezier with constraints
        c = curve_constraints()
        c.init_vel = matrix([0., 1., 1.]).T
        c.end_vel = matrix([0., 1., 1.]).T
        c.init_acc = matrix([0., 1., -1.]).T
        c.end_acc = matrix([0., 100., 1.]).T
        #Check derivate with constraints
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).T
        a = bezier3(waypoints, c)
        assert_allclose(a.derivate(0, 1), c.init_vel)
        assert_allclose(a.derivate(1, 2), c.end_acc)

        # TESTING POLYNOM FUNCTION
        # - Functions : constructor, min, max, derivate
        a = polynom(waypoints)
        a = polynom(waypoints, -1., 3.)
        a.min()
        a.max()
        a(0.4)
        assert_allclose(a.derivate(0.4, 0), a(0.4))
        a.derivate(0.4, 2)

        # TESTING EXACT_CUBIC FUNCTION
        # - Functions : constructor, min, max, derivate
        a = exact_cubic(waypoints, time_waypoints)
        a.min()
        a.max()
        a(0.4)
        assert_allclose(a.derivate(0.4, 0), a(0.4))
        a.derivate(0.4, 2)

        # TESTING SPLINE_DERIV_CONSTRAINTS
        # - Functions : constructor, min, max, derivate
        c = curve_constraints()
        c.init_vel
        c.end_vel
        c.init_acc
        c.end_acc
        c.init_vel = matrix([0., 1., 1.]).T
        c.end_vel = matrix([0., 1., 1.]).T
        c.init_acc = matrix([0., 1., 1.]).T
        c.end_acc = matrix([0., 1., 1.]).T
        a = spline_deriv_constraint(waypoints, time_waypoints)
        a = spline_deriv_constraint(waypoints, time_waypoints, c)

        # CONVERTING BEZIER TO POLYNOM
        a = bezier3(waypoints)
        a_pol = from_bezier(a)
        assert_allclose(a(0.3), a_pol(0.3))


if __name__ == '__main__':
    unittest.main()
