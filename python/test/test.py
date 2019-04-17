# FOR PRINT, TO REMOVE
import StringIO
import sys
# END FOR PRINT


from numpy import matrix
from numpy.linalg import norm

import unittest

from curves import bezier3, bezier6, curve_constraints, exact_cubic, from_bezier, polynom, spline_deriv_constraint


class TestCurves(unittest.TestCase):
	#def print_str(self, inStr):
	#	print inStr
	#	return

	def test_bezier(self):
		# To test :
		# - Functions : constructor, min, max, derivate,compute_derivate, compute_primitive
		# - Variables : degree, nbWayPoints
		__EPS = 1e-6
		waypoints = matrix([[1., 2., 3.]]).T
		a = bezier3(waypoints,2.)
		t = 0.
		while t < 2.:
			self.assertTrue (norm(a(t) - matrix([1., 2., 3.]).T) < __EPS)
			t += 0.1
		waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
		waypoints6 = matrix([[1., 2., 3., 7., 5., 5.], [4., 5., 6., 4., 5., 6.]]).transpose()
		time_waypoints = matrix([0., 1.]).transpose()
		# Create bezier6 and bezier3
		a = bezier6(waypoints6)
		a = bezier3(waypoints, 3.)
		# Test : Degree, min, max, derivate
		#self.print_str(("test 1")
		self.assertEqual (a.degree, a.nbWaypoints - 1)
		a.min()
		a.max()
		a(0.4)
		self.assertTrue ((a.derivate(0.4, 0) == a(0.4)).all())
		a.derivate(0.4, 2)
		a = a.compute_derivate(100)
		prim = a.compute_primitive(1)
		# Check primitive and derivate - order 1
		for i in range(10):
		    t = float(i) / 10.
		    self.assertTrue ((a(t) == prim.derivate(t, 1)).all())
		self.assertTrue ((prim(0) == matrix([0., 0., 0.])).all())
		# Check primitive and derivate - order 2
		prim = a.compute_primitive(2)
		for i in range(10):
		    t = float(i) / 10.
		    self.assertTrue ((a(t) == prim.derivate(t, 2)).all())
		self.assertTrue ((prim(0) == matrix([0., 0., 0.])).all())
		# Create new bezier3 curve
		waypoints = matrix([[1., 2., 3.], [4., 5., 6.], [4., 5., 6.], [4., 5., 6.], [4., 5., 6.]]).transpose()
		a0 = bezier3(waypoints)
		a1 = bezier3(waypoints, 3.)
		prim0 = a0.compute_primitive(1)
		prim1 = a1.compute_primitive(1)
		# Check change in argument time_t of bezier3
		for i in range(10):
		    t = float(i) / 10.
		    self.assertTrue (norm(a0(t) - a1(3 * t)) < __EPS)
		    self.assertTrue (norm(a0.derivate(t, 1) - a1.derivate(3 * t, 1) * 3.) < __EPS)
		    self.assertTrue (norm(a0.derivate(t, 2) - a1.derivate(3 * t, 2) * 9.) < __EPS)
		    self.assertTrue (norm(prim0(t) - prim1(t * 3) / 3.) < __EPS)
		self.assertTrue ((prim(0) == matrix([0., 0., 0.])).all())
		# testing bezier with constraints
		c = curve_constraints()
		c.init_vel = matrix([0., 1., 1.]).transpose()
		c.end_vel = matrix([0., 1., 1.]).transpose()
		c.init_acc = matrix([0., 1., -1.]).transpose()
		c.end_acc = matrix([0., 100., 1.]).transpose()
		#Check derivate with constraints
		waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
		a = bezier3(waypoints, c)
		self.assertTrue (norm(a.derivate(0, 1) - c.init_vel) < 1e-10)
		self.assertTrue (norm(a.derivate(1, 2) - c.end_acc) < 1e-10)
		return

	def test_polynom(self):
		# To test :
		# - Functions : constructor, min, max, derivate
		waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
		a = polynom(waypoints)
		a = polynom(waypoints, -1., 3.)
		a.min()
		a.max()
		a(0.4)
		self.assertTrue ((a.derivate(0.4, 0) == a(0.4)).all())
		a.derivate(0.4, 2)
		return

	def test_exact_cubic(self):
		# To test :
		# - Functions : constructor, min, max, derivate
		waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
		time_waypoints = matrix([0., 1.]).transpose()
		a = exact_cubic(waypoints, time_waypoints)
		a.min()
		a.max()
		a(0.4)
		self.assertTrue ((a.derivate(0.4, 0) == a(0.4)).all())
		a.derivate(0.4, 2)
		return

	def test_spline_deriv_constraint(self):
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
		a = spline_deriv_constraint(waypoints, time_waypoints)
		a = spline_deriv_constraint(waypoints, time_waypoints, c)
		return

	def test_from_bezier(self):
		# converting bezier to polynom
		__EPS = 1e-6
		waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
		a = bezier3(waypoints)
		a_pol = from_bezier(a)
		self.assertTrue (norm(a(0.3) - a_pol(0.3)) < __EPS)
		return

if __name__ == '__main__':
    unittest.main()
