import unittest
import os

from numpy import matrix
from numpy.linalg import norm

#from curves import ( serialize_polynomial, deserialize_polynomial, serialize_piecewise_polynomial_curve, deserialize_piecewise_polynomial_curve )

from curves import (bezier3, bezier6, bezier_from_hermite, bezier_from_polynomial, cubic_hermite_spline,
                    curve_constraints, exact_cubic, hermite_from_bezier, hermite_from_polynomial,
                    piecewise_bezier3_curve, piecewise_bezier6_curve, piecewise_cubic_hermite_curve,
                    piecewise_polynomial_curve, polynomial, polynomial_from_bezier, polynomial_from_hermite)

#import curves

class TestCurves(unittest.TestCase):
    #def print_str(self, inStr):
    #   print inStr
    #   return

    def test_bezier(self):
        # To test :
        # - Functions : constructor, min, max, derivate,compute_derivate, compute_primitive
        # - Variables : degree, nbWayPoints
        __EPS = 1e-6
        waypoints = matrix([[1., 2., 3.]]).T
        a = bezier3(waypoints, 0., 2.)
        t = 0.
        while t < 2.:
            self.assertTrue(norm(a(t) - matrix([1., 2., 3.]).T) < __EPS)
            t += 0.1
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        waypoints6 = matrix([[1., 2., 3., 7., 5., 5.], [4., 5., 6., 4., 5., 6.]]).transpose()
        time_waypoints = matrix([0., 1.]).transpose()
        # Create bezier6 and bezier3
        a6 = bezier6(waypoints6)
        a = bezier3(waypoints, 0., 3.)
        # Test : Degree, min, max, derivate
        #self.print_str(("test 1")
        self.assertEqual(a.degree, a.nbWaypoints - 1)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        a = a.compute_derivate(100)
        prim = a.compute_primitive(1)
        # Check primitive and derivate - order 1
        for i in range(10):
            t = float(i) / 10.
            self.assertTrue((a(t) == prim.derivate(t, 1)).all())
        self.assertTrue((prim(0) == matrix([0., 0., 0.])).all())
        # Check primitive and derivate - order 2
        prim = a.compute_primitive(2)
        for i in range(10):
            t = float(i) / 10.
            self.assertTrue((a(t) == prim.derivate(t, 2)).all())
        self.assertTrue((prim(0) == matrix([0., 0., 0.])).all())
        # Create new bezier3 curve
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.], [4., 5., 6.], [4., 5., 6.], [4., 5., 6.]]).transpose()
        a0 = bezier3(waypoints)
        a1 = bezier3(waypoints, 0., 3.)
        prim0 = a0.compute_primitive(1)
        prim1 = a1.compute_primitive(1)
        # Check change in argument time_t of bezier3
        for i in range(10):
            t = float(i) / 10.
            self.assertTrue(norm(a0(t) - a1(3 * t)) < __EPS)
            self.assertTrue(norm(a0.derivate(t, 1) - a1.derivate(3 * t, 1) * 3.) < __EPS)
            self.assertTrue(norm(a0.derivate(t, 2) - a1.derivate(3 * t, 2) * 9.) < __EPS)
            self.assertTrue(norm(prim0(t) - prim1(t * 3) / 3.) < __EPS)
        self.assertTrue((prim(0) == matrix([0., 0., 0.])).all())
        # testing bezier with constraints
        c = curve_constraints()
        c.init_vel = matrix([0., 1., 1.]).transpose()
        c.end_vel = matrix([0., 1., 1.]).transpose()
        c.init_acc = matrix([0., 1., -1.]).transpose()
        c.end_acc = matrix([0., 100., 1.]).transpose()
        #Check derivate with constraints
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        a = bezier3(waypoints, c)
        self.assertTrue(norm(a.derivate(0, 1) - c.init_vel) < 1e-10)
        self.assertTrue(norm(a.derivate(1, 2) - c.end_acc) < 1e-10)

        # Test serialization : bezier 3
        a.saveAsText("serialization_curve.test")
        #waypoints = matrix([[0,0,0,], [0,0,0,]]).transpose()
        b = bezier3()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())

        # Test serialization : Bezier 6
        a6.saveAsText("serialization_curve.test")
        b6 = bezier6()
        b6.loadFromText("serialization_curve.test")
        self.assertTrue((a6(0.4) == b6(0.4)).all())
        os.remove("serialization_curve.test")
        return

    def test_polynomial(self):
        # To test :
        # - Functions : constructor, min, max, derivate, serialize, deserialize
        waypoints_0 = matrix([[0.,0.,0.], [0.,0.,0.]]).transpose()
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        a = polynomial(waypoints)  # Defined on [0.,1.]
        a = polynomial(waypoints, -1., 3.)  # Defined on [-1.,3.]
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        # Test serialization
        a.saveAsText("serialization_curve.test")
        b = polynomial()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")
        return

    def test_cubic_hermite_spline(self):
        points = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        tangents = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        time_points = matrix([0., 1.]).transpose()
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
        return

    def test_piecewise_polynomial_curve(self):
        # To test :
        # - Functions : constructor, min, max, derivate, add_curve, is_continuous, serialize, deserialize
        waypoints0 = matrix([[0., 0., 0.]]).transpose()
        waypoints1 = matrix([[1., 1., 1.]]).transpose()
        waypoints2 = matrix([[1., 1., 1.], [1., 1., 1.]]).transpose()
        pol0 = polynomial(waypoints0, 0., 0.1)
        a = polynomial(waypoints1, 0., 1.)
        b = polynomial(waypoints2, 1., 3.)
        pc = piecewise_polynomial_curve(a)
        pc.add_curve(b)
        pc.min()
        pc.max()
        pc(0.4)
        self.assertTrue((pc.derivate(0.4, 0) == pc(0.4)).all())
        pc.derivate(0.4, 2)
        pc.is_continuous(0)
        pc.is_continuous(1)
        # Test serialization
        pc.saveAsText("serialization_pc.test")
        pc_test = piecewise_polynomial_curve()
        pc_test.loadFromText("serialization_pc.test")
        self.assertTrue((pc(0.4) == pc_test(0.4)).all())
        os.remove("serialization_pc.test")
        return

    def test_piecewise_bezier3_curve(self):
        # To test :
        # - Functions : constructor, min, max, derivate, add_curve, is_continuous
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        a = bezier3(waypoints, 0., 1.)
        b = bezier3(waypoints, 1., 2.)
        pc = piecewise_bezier3_curve(a)
        pc.add_curve(b)
        pc.min()
        pc.max()
        pc(0.4)
        self.assertTrue((pc.derivate(0.4, 0) == pc(0.4)).all())
        pc.derivate(0.4, 2)
        pc.is_continuous(0)
        pc.is_continuous(1)
        # Test serialization
        pc.saveAsText("serialization_pc.test")
        pc_test = piecewise_bezier3_curve()
        pc_test.loadFromText("serialization_pc.test")
        self.assertTrue((pc(0.4) == pc_test(0.4)).all())
        os.remove("serialization_pc.test")
        return

    def test_piecewise_bezier6_curve(self):
        # To test :
        # - Functions : constructor, min, max, derivate, add_curve, is_continuous
        waypoints = matrix([[1., 2., 3., 7., 5., 5.], [4., 5., 6., 4., 5., 6.]]).transpose()
        a = bezier6(waypoints, 0., 1.)
        b = bezier6(waypoints, 1., 2.)
        pc = piecewise_bezier6_curve(a)
        pc.add_curve(b)
        pc.min()
        pc.max()
        pc(0.4)
        self.assertTrue((pc.derivate(0.4, 0) == pc(0.4)).all())
        pc.derivate(0.4, 2)
        pc.is_continuous(0)
        pc.is_continuous(1)
        # Test serialization
        pc.saveAsText("serialization_pc.test")
        pc_test = piecewise_bezier6_curve()
        pc_test.loadFromText("serialization_pc.test")
        self.assertTrue((pc(0.4) == pc_test(0.4)).all())
        os.remove("serialization_pc.test")
        return

    def test_piecewise_cubic_hermite_curve(self):
        # To test :
        # - Functions : constructor, min, max, derivate, add_curve, is_continuous
        points = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        tangents = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        time_points0 = matrix([0., 1.]).transpose()
        time_points1 = matrix([1., 2.]).transpose()
        a = cubic_hermite_spline(points, tangents, time_points0)
        b = cubic_hermite_spline(points, tangents, time_points1)
        pc = piecewise_cubic_hermite_curve(a)
        pc.add_curve(b)
        pc.min()
        pc.max()
        pc(0.4)
        self.assertTrue((pc.derivate(0.4, 0) == pc(0.4)).all())
        pc.derivate(0.4, 2)
        pc.is_continuous(0)
        pc.is_continuous(1)
        # Test serialization
        pc.saveAsText("serialization_pc.test")
        pc_test = piecewise_cubic_hermite_curve()
        pc_test.loadFromText("serialization_pc.test")
        self.assertTrue((pc(0.4) == pc_test(0.4)).all())
        os.remove("serialization_pc.test")
        return

    def test_exact_cubic(self):
        # To test :
        # - Functions : constructor, min, max, derivate, getNumberSplines, getSplineAt
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        time_waypoints = matrix([0., 1.]).transpose()
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
        return

    def test_exact_cubic_constraint(self):
        # To test :
        # - Functions : constructor, min, max, derivate
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        time_waypoints = matrix([0., 1.]).transpose()
        c = curve_constraints()
        c.init_vel
        c.end_vel
        c.init_acc
        c.end_acc
        c.init_vel = matrix([0., 1., 1.]).transpose()
        c.end_vel = matrix([0., 1., 1.]).transpose()
        c.init_acc = matrix([0., 1., 1.]).transpose()
        c.end_acc = matrix([0., 1., 1.]).transpose()
        a = exact_cubic(waypoints, time_waypoints)
        a = exact_cubic(waypoints, time_waypoints, c)
        return

    def test_conversion_curves(self):
        __EPS = 1e-6
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        a = bezier3(waypoints)
        # converting bezier to polynomial
        a_pol = polynomial_from_bezier(a)
        self.assertTrue(norm(a(0.3) - a_pol(0.3)) < __EPS)
        # converting polynomial to hermite
        a_chs = hermite_from_polynomial(a_pol)
        self.assertTrue(norm(a_chs(0.3) - a_pol(0.3)) < __EPS)
        # converting hermite to bezier
        a_bc = bezier_from_hermite(a_chs)
        self.assertTrue(norm(a_chs(0.3) - a_bc(0.3)) < __EPS)
        self.assertTrue(norm(a(0.3) - a_bc(0.3)) < __EPS)
        # converting bezier to hermite
        a_chs = hermite_from_bezier(a)
        self.assertTrue(norm(a(0.3) - a_chs(0.3)) < __EPS)
        # converting hermite to polynomial
        a_pol = polynomial_from_hermite(a_chs)
        self.assertTrue(norm(a_pol(0.3) - a_chs(0.3)) < __EPS)
        # converting polynomial to bezier
        a_bc = bezier_from_polynomial(a_pol)
        self.assertTrue(norm(a_bc(0.3) - a_pol(0.3)) < __EPS)
        self.assertTrue(norm(a(0.3) - a_bc(0.3)) < __EPS)
        return


if __name__ == '__main__':
    unittest.main()
