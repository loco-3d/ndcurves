import numpy as np
from curves import bezier, curve_constraints
from curves.optimization import (
    constraint_flag,
    problem_definition,
    setup_control_points,
)
from numpy import array, identity
from qp import quadprog_solve_qp, to_least_square

from .plot_bezier import plotBezier, plt

np.set_printoptions(formatter={"float": lambda x: f"{x:0.1f}"})


# points given as pairs (pt, t)
# does not try to fit anything
def default_fit(
    degree,
    ptsTime,
    dim=3,
    totalTime=1.0,
    constraintFlag=constraint_flag.NONE,
    curveConstraints=None,
):
    pD = None
    if curveConstraints is None:
        pD = problem_definition(dim)
    else:
        pD = problem_definition(curveConstraints)
    # set initial and goal positions to the ones of the list
    # (might not be used if constraints do not correspond)
    pD.init_pos = array([ptsTime[0][0]]).T
    pD.end_pos = array([ptsTime[-1][0]]).T
    # assign curve constraints
    pD.totalTime = totalTime
    pD.degree = degree
    pD.flag = constraintFlag
    problem = setup_control_points(pD)
    bezierLinear = problem.bezier()
    pD.totalTime = totalTime
    pD.degree = degree
    problem = setup_control_points(pD)
    bezierLinear = problem.bezier()

    allsEvals = [(bezierLinear(time), pt) for (pt, time) in ptsTime]
    allLeastSquares = [to_least_square(el.B(), el.c() + pt) for (el, pt) in allsEvals]
    Ab = [sum(x) for x in zip(*allLeastSquares)]

    res = quadprog_solve_qp(
        Ab[0] + identity(problem.numVariables * dim) * 0.0001, Ab[1]
    )
    return bezierLinear.evaluate(res.reshape((-1, 1)))


# return A, b

if __name__ == "__main__":
    waypoints = array(
        [
            [1.0, 2.0, 3.0, 0.0],
            [-4.0, -5.0, -6.0, 0.0],
            [4.0, 5.0, 6.0, 0.0],
            [7.0, 8.0, 9.0, 0.0],
        ]
    ).transpose()
    a = bezier(waypoints)
    totalTime = 1.0
    degree = a.nbWaypoints - 1
    ptsTime = [(a(float(i) / 10.0), float(i) / 10.0) for i in range(11)]
    dim = waypoints.shape[0]

    c = curve_constraints(dim)
    c.init_vel = array([[0.0, 0.0, 0.0, 0.0]]).T

    def plot_fit_and_original(ax, plotControlPoints=False):
        plotBezier(a, ax=ax)
        plotBezier(bFit, ax=ax, color="g")
        if plotControlPoints:
            plotControlPoints(a, ax=ax)
            plotControlPoints(bFit, ax=ax, color="g")

    fig = plt.figure()

    bFit = default_fit(degree, ptsTime, dim=dim)
    ax = fig.add_subplot(221, projection="3d")
    plot_fit_and_original(ax)

    bFit = default_fit(degree - 1, ptsTime, dim=dim)
    ax = fig.add_subplot(222, projection="3d")
    plot_fit_and_original(ax)

    bFit = default_fit(
        degree - 1,
        ptsTime,
        dim=dim,
        constraintFlag=constraint_flag.END_POS | constraint_flag.INIT_POS,
        curveConstraints=c,
    )
    ax = fig.add_subplot(223, projection="3d")
    plot_fit_and_original(ax)

    bFit = default_fit(
        degree + 30,
        ptsTime,
        dim=dim,
        constraintFlag=constraint_flag.END_POS
        | constraint_flag.INIT_VEL
        | constraint_flag.INIT_POS,
        curveConstraints=c,
    )
    ax = fig.add_subplot(224, projection="3d")
    plot_fit_and_original(ax)

    plt.show(block=False)
