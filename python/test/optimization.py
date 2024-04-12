import unittest

from numpy import array, zeros
from numpy.linalg import norm

from ndcurves.optimization import (
    constraint_flag,
    generate_integral_problem,
    integral_cost_flag,
    problem_definition,
    setup_control_points,
)


class TestProblemDefinition(unittest.TestCase):
    # generate problem data
    def test_problem_definition(self):
        pD = problem_definition(3)
        pD.init_pos
        pD.init_vel = array([0.0, 1.0, 1.0]).T
        pD.end_vel = array([0.0, 1.0, 1.0]).T
        pD.init_acc = array([0.0, 1.0, -1.0]).T
        pD.end_acc = array([0.0, 0.0, 0]).T
        pD.init_pos = array([[0.0, 0.0, 0.0]]).T
        pD.end_pos = array([[1.0, 1.0, 1.0]]).T
        pD.flag = constraint_flag.INIT_VEL | constraint_flag.INIT_POS

        generate_integral_problem(pD, integral_cost_flag.ACCELERATION)

        # generate problem
        problem = setup_control_points(pD)
        bezierLinear = problem.bezier()
        bezierFixed = bezierLinear.evaluate(array([zeros(12)]).T)
        self.assertTrue(bezierFixed.nbWaypoints == pD.degree + 1)
        self.assertTrue(norm(bezierFixed(0.0) - pD.init_pos) <= 0.001)
        self.assertTrue(norm(bezierFixed.derivate(0.0, 1) - pD.init_vel) <= 0.001)


if __name__ == "__main__":
    unittest.main()
