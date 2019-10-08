import os
import unittest
import eigenpy
eigenpy.switchToNumpyArray()

from numpy import array_equal, isclose, matrix, random, array, zeros
from numpy.linalg import norm

from curves import (bezier, bezier_from_hermite, bezier_from_polynomial, cubic_hermite_spline, curve_constraints,
                    exact_cubic, hermite_from_bezier, hermite_from_polynomial, piecewise_bezier_curve,
                    piecewise_cubic_hermite_curve, piecewise_polynomial_curve, polynomial, polynomial_from_bezier,
                    polynomial_from_hermite)

from curves.optimization import(constraint_flag, integral_cost_flag, quadratic_problem, setup_control_points, generate_problem, generate_integral_problem, problemData, problemDefinition)


#generate problem data

c = curve_constraints()
c.init_vel = matrix([0., 1., 1.]).transpose()
c.end_vel = matrix([0., 1., 1.]).transpose()
c.init_acc = matrix([0., 1., -1.]).transpose()
c.end_acc = matrix([0., 0., 0]).transpose()

pD = problemDefinition()
pD.start
pD.curveConstraints = c
pD.start = array([[0.,0.,0.]]).T
pD.end = array([[1.,1.,1.]]).T
pD.flag = constraint_flag.INIT_VEL | constraint_flag.INIT_POS

#~ generate_integral_problem(pD,integral_cost_flag.ACCELERATION)

#generate problem
problem = setup_control_points(pD)
bezier = problem.bezier()
res = bezier.evaluate(array([zeros(12)]).T)(0.)
