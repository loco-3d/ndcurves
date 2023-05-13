import unittest

from numpy import array, dot, identity, zeros

# importing the bezier curve class
from ndcurves import bezier


# dummy methods
def plot(*karrgs):
    pass


class TestNotebook(unittest.TestCase):
    # def print_str(self, inStr):
    #   print inStr
    #   return

    def test_notebook(self):
        print("test_notebook")

        # We describe a degree 3 curve as a Bezier curve with 4 control points
        waypoints = array(
            [[1.0, 2.0, 3.0], [-4.0, -5.0, -6.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
        ).transpose()
        ref = bezier(waypoints)

        numSamples = 10
        fNumSamples = float(numSamples)
        ptsTime = [
            (ref(float(t) / fNumSamples), float(t) / fNumSamples)
            for t in range(numSamples + 1)
        ]

        from ndcurves.optimization import problem_definition, setup_control_points

        # dimension of our problem (here 3 as our curve is 3D)
        dim = 3
        refDegree = 3
        pD = problem_definition(dim)
        # we want to fit a curve of the same degree
        # as the reference curve for the sanity check
        pD.degree = refDegree

        # generates the variable bezier curve with the parameters of problemDefinition
        problem = setup_control_points(pD)
        # for now we only care about the curve itself
        variableBezier = problem.bezier()

        variableBezier(0.0)

        # least square form of ||Ax-b||**2
        def to_least_square(A, b):
            return dot(A.T, A), -dot(A.T, b)

        def genCost(variableBezier, ptsTime):
            # first evaluate variableBezier for each time sampled
            allsEvals = [(variableBezier(time), pt) for (pt, time) in ptsTime]
            # then compute the least square form of the cost for each points
            allLeastSquares = [
                to_least_square(el.B(), el.c() + pt) for (el, pt) in allsEvals
            ]
            # and finally sum the costs
            Ab = [sum(x) for x in zip(*allLeastSquares)]
            return Ab[0], Ab[1]

        A, b = genCost(variableBezier, ptsTime)

        def quadprog_solve_qp(P, q, G=None, h=None, C=None, d=None, verbose=False):
            return zeros(P.shape[0])

        res = quadprog_solve_qp(A, b)

        def evalAndPlot(variableBezier, res):
            fitBezier = variableBezier.evaluate(res.reshape((-1, 1)))
            return fitBezier

        fitBezier = evalAndPlot(variableBezier, res)

        pD.degree = refDegree - 1

        problem = setup_control_points(pD)
        variableBezier = problem.bezier()

        A, b = genCost(variableBezier, ptsTime)
        res = quadprog_solve_qp(A, b)
        fitBezier = evalAndPlot(variableBezier, res)

        from ndcurves.optimization import constraint_flag

        pD.flag = constraint_flag.INIT_POS | constraint_flag.END_POS
        # set initial position
        pD.init_pos = array([ptsTime[0][0]]).T
        # set end position
        pD.end_pos = array([ptsTime[-1][0]]).T
        problem = setup_control_points(pD)
        variableBezier = problem.bezier()

        prob = setup_control_points(pD)
        variableBezier = prob.bezier()
        A, b = genCost(variableBezier, ptsTime)
        res = quadprog_solve_qp(A, b)
        evalAndPlot(variableBezier, res)

        # values are 0 by default, so if the constraint is zero this can be skipped
        pD.init_vel = array([[0.0, 0.0, 0.0]]).T
        pD.init_acc = array([[0.0, 0.0, 0.0]]).T
        pD.end_vel = array([[0.0, 0.0, 0.0]]).T
        pD.end_acc = array([[0.0, 0.0, 0.0]]).T
        pD.flag = (
            constraint_flag.END_POS
            | constraint_flag.INIT_POS
            | constraint_flag.INIT_VEL
            | constraint_flag.END_VEL
            | constraint_flag.INIT_ACC
            | constraint_flag.END_ACC
        )

        err = False
        try:
            prob = setup_control_points(pD)
        except RuntimeError:
            err = True
        assert err

        pD.degree = refDegree + 4
        prob = setup_control_points(pD)
        variableBezier = prob.bezier()
        A, b = genCost(variableBezier, ptsTime)
        res = quadprog_solve_qp(A, b)
        fitBezier = evalAndPlot(variableBezier, res)

        pD.degree = refDegree + 60
        prob = setup_control_points(pD)
        variableBezier = prob.bezier()
        A, b = genCost(variableBezier, ptsTime)
        # regularization matrix
        reg = identity(A.shape[1]) * 0.001
        res = quadprog_solve_qp(A + reg, b)
        fitBezier = evalAndPlot(variableBezier, res)

        # set initial / terminal constraints
        pD.flag = constraint_flag.END_POS | constraint_flag.INIT_POS
        pD.degree = refDegree
        prob = setup_control_points(pD)
        variableBezier = prob.bezier()

        # get value of the curve first order derivative at t = 0.8
        t08Constraint = variableBezier.derivate(0.8, 1)
        target = zeros(3)

        A, b = genCost(variableBezier, ptsTime)
        # solve optimization problem with quadprog

        res = quadprog_solve_qp(A, b, C=t08Constraint.B(), d=target - t08Constraint.c())
        fitBezier = evalAndPlot(variableBezier, res)

        # returns a curve composed of the split curves, 2 in our case
        piecewiseCurve = ref.split(array([[0.6]]).T)

        # displaying the obtained curves

        # first, split the variable curve
        piecewiseCurve = variableBezier.split(array([[0.4, 0.8]]).T)

        constrainedCurve = piecewiseCurve.curve_at_index(1)

        # find the number of variables
        problemSize = prob.numVariables * dim
        # find the number of constraints, as many as waypoints
        nConstraints = constrainedCurve.nbWaypoints

        waypoints = constrainedCurve.waypoints()

        ineqMatrix = zeros((nConstraints, problemSize))
        ineqVector = zeros(nConstraints)

        # finding the z equation of each control point
        for i in range(nConstraints):
            wayPoint = constrainedCurve.waypointAtIndex(i)
            ineqMatrix[i, :] = wayPoint.B()[2, :]
            ineqVector[i] = -wayPoint.c()[2]

        res = quadprog_solve_qp(A, b, G=ineqMatrix, h=ineqVector)
        fitBezier = variableBezier.evaluate(res.reshape((-1, 1)))
        fitBezier


if __name__ == "__main__":
    unittest.main()
