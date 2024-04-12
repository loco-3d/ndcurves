import uuid

import matplotlib.pyplot as plt
import numpy as np
from numpy import array, cross, identity, vstack, zeros
from numpy.linalg import norm

from .convex_hull import genFromLine
from .plot_cord import plotBezier, plotControlPoints, plotPoly
from .qp import quadprog_solve_qp
from .qp_cord import accelerationcost, concat, concatvec, lineConstraint
from .varBezier import varBezier

__EPS = 1e-6

idxFile = 0
rng = np.random.default_rng()
uuid.uuid4().hex.upper()[0:6]


# ####################### generate problems ########################
def genBezierInput(numvars=3):
    valDep = array([rng.uniform(0.0, 1.0), rng.uniform(0.0, 5.0), 0.0])
    valEnd = array([rng.uniform(5.0, 10.0), rng.uniform(0.0, 5.0), 0.0])
    return varBezier([valDep] + ["" for _ in range(numvars)] + [valEnd], 1.0)


def genSplit(numCurves):
    splits = []
    lastval = rng.uniform(0.0, 1.0)
    for i in range(numCurves):
        splits += [lastval]
        lastval += rng.uniform(0.0, 1.0)
    return [el / lastval for el in splits[:-1]]


def getRightMostLine(ptList):
    pt1 = array([0.0, 0.0, 0.0])
    pt2 = array([0.0, 0.0, 0.0])
    for pt in ptList:
        if pt[0] > pt1[0]:
            pt1 = pt
        elif pt[0] > pt2[0]:
            pt2 = pt
    if pt1[1] < pt2[1]:
        return [pt2, pt1]
    else:
        return [pt1, pt2]


# ####################### generate problems ########################

# ####################### solve a given problem ########################


# inequality constraint from line
def getLineFromSegment(line):
    a = line[0]
    b = line[1]
    c = a.copy()
    c[2] = 1.0
    normal = cross((b - a), (c - a))
    normal /= norm(normal)
    # get inequality
    coeff = normal
    rhs = a.dot(normal)
    return (coeff, array([rhs]))


def computeTrajectory(bezVar, splits, save, filename=uuid.uuid4().hex.upper()[0:6]):
    global idxFile
    colors = ["b", "g", "r", "c", "m", "y", "k", "w"]
    subs = bezVar.split(splits)
    # generate random constraints for each line
    line_current = [array([1.0, 1.0, 0.0]), array([0.0, 1.0, 0.0])]
    # qp vars
    P = accelerationcost(bezVar)[0]
    P = P + identity(P.shape[0]) * 0.0001
    q = zeros(3 * bezVar.bezier.nbWaypoints)
    q[-1] = -1
    G = zeros([2, q.shape[0]])
    h = zeros(2)
    G[0, -1] = 1
    h[0] = 1.0
    G[1, -1] = -1
    h[1] = 0.0
    dimExtra = 0
    for i, bez in enumerate(subs):
        color = colors[i]
        init_points = []
        if i == 0:
            init_points = [bezVar.waypoints()[1][:3, 0][:2]]
        if i == len(subs) - 1:
            init_points = [bezVar.waypoints()[1][-3:, -1][:2]]
        lines, ptList = genFromLine(line_current, 5, [[0, 5], [0, 5]], init_points)
        matineq0 = None
        vecineq0 = None
        for line in lines:
            (mat, vec) = getLineFromSegment(line)
            (mat, vec) = lineConstraint(bez, mat, vec, dimExtra)
            (matineq0, vecineq0) = (concat(matineq0, mat), concatvec(vecineq0, vec))
        ineq = (matineq0, vecineq0)
        G = vstack([G, ineq[0]])
        h = concatvec(h, ineq[1])
        line_current = getRightMostLine(ptList)
        plotPoly(lines, color)
    C = None
    d = None
    try:
        res = quadprog_solve_qp(P, q, G=G, h=h, C=C, d=d)
        # plot bezier
        for i, bez in enumerate(subs):
            color = colors[i]
            test = bez.toBezier3(res[:])
            plotBezier(test, color)
        if save:
            plt.savefig(filename + str(idxFile))
        # plot subs control points
        for i, bez in enumerate(subs):
            color = colors[i]
            test = bez.toBezier3(res[:])
            plotBezier(test, color)
            plotControlPoints(test, color)
        if save:
            plt.savefig("subcp" + filename + str(idxFile))
        final = bezVar.toBezier3(res[:])
        plotControlPoints(final, "black")
        if save:
            plt.savefig("cp" + filename + str(idxFile))
        idxFile += 1
        if save:
            plt.close()
        else:
            plt.show()
        return final
    except ValueError:
        plt.close()


# ####################### solve a given problem ########################


# solve and gen problem
def gen(save=False):
    testConstant = genBezierInput(20)
    splits = genSplit(4)
    return computeTrajectory(testConstant, splits, save), testConstant


res = None
for i in range(1):
    res = gen(False)
    # if res[0] != None:
    # break
