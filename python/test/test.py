#!/usr/bin/env python

import unittest

from numpy import matrix
from numpy.linalg import norm
from numpy.testing import assert_allclose

from hpp_spline import bezier, bezier6, curve_constraints, exact_cubic, from_bezier, polynom, spline_deriv_constraint


class TestSpline(unittest.TestCase):
    def test_spline(self):
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).T
        waypoints6 = matrix([[1., 2., 3., 7., 5., 5.], [4., 5., 6., 4., 5., 6.]]).T
        time_waypoints = matrix([0., 1.]).T

        # testing bezier curve
        a = bezier6(waypoints6)
        a = bezier(waypoints, 3.)

        self.assertEqual(a.degree, a.nbWaypoints - 1)
        a.min()
        a.max()
        a(0.4)
        assert_allclose(a.derivate(0.4, 0), a(0.4))
        a.derivate(0.4, 2)
        a = a.compute_derivate(100)

        prim = a.compute_primitive(1)

        for i in range(10):
            t = float(i) / 10.
            assert_allclose(a(t), prim.derivate(t, 1))
        assert_allclose(prim(0), matrix([0., 0., 0.]).T)

        prim = a.compute_primitive(2)
        for i in range(10):
            t = float(i) / 10.
            assert_allclose(a(t), prim.derivate(t, 2), atol=1e-20)
        assert_allclose(prim(0), matrix([0., 0., 0.]).T)

        waypoints = matrix([[1., 2., 3.], [4., 5., 6.], [4., 5., 6.], [4., 5., 6.], [4., 5., 6.]]).T
        a0 = bezier(waypoints)
        a1 = bezier(waypoints, 3.)
        prim0 = a0.compute_primitive(1)
        prim1 = a1.compute_primitive(1)

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

        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).T
        a = bezier(waypoints, c)
        assert_allclose(a.derivate(0, 1), c.init_vel)
        assert_allclose(a.derivate(1, 2), c.end_acc)

        # testing polynom function
        a = polynom(waypoints)
        a = polynom(waypoints, -1., 3.)
        a.min()
        a.max()
        a(0.4)
        assert_allclose(a.derivate(0.4, 0), a(0.4))
        a.derivate(0.4, 2)

        # testing exact_cubic function
        a = exact_cubic(waypoints, time_waypoints)
        a.min()
        a.max()
        a(0.4)
        assert_allclose(a.derivate(0.4, 0), a(0.4))
        a.derivate(0.4, 2)

        # testing spline_deriv_constraints
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

        # converting bezier to polynom

        a = bezier(waypoints)
        a_pol = from_bezier(a)
        assert_allclose(a(0.3), a_pol(0.3))


if __name__ == '__main__':
    unittest.main()
