
from numpy import array_equal, isclose, matrix, random, zeros, array
from numpy.linalg import norm
import numpy as np
np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)})

import eigenpy
eigenpy.switchToNumpyArray()

from curves import (bezier, curve_constraints)
from curves.optimization import(constraint_flag, integral_cost_flag, setup_control_points, problemData, problemDefinition)

from qp import to_least_square, quadprog_solve_qp


  

#points given as pairs (pt, t)
#does not try to fit anything
def default_fit_3d(degree, ptsTime, totalTime =1., zeroInitVel = False):
    pD = problemDefinition()
    if zeroInitVel:
        pD.flag = constraint_flag.INIT_VEL | constraint_flag.INIT_POS
        c = curve_constraints()
        pD.start = array([ptsTime[0][0]]).T
        c.init_vel = array([[0., 0., 0.]]).T
        pD.curveConstraints = c
    pD.totalTime = totalTime
    pD.degree = degree
    problem = setup_control_points(pD)
    bezierLinear = problem.bezier()

    nVar = problem.numVariables * 3 #dimension 3
    A, b = zeros((nVar,nVar)), zeros(nVar)
    for (pt, time) in ptsTime:
        assert time <= totalTime, "total time inferior to sampling"
        exprAtT = bezierLinear(time)
        nA, nb = to_least_square(exprAtT.A, exprAtT.b + pt)
        A += nA
        b += nb
    
    res = quadprog_solve_qp(A, b)
    return bezierLinear.evaluate(res.reshape((-1,1)) ) 


#~ return A, b

if __name__ == '__main__':

    from plot_bezier import *
    
    waypoints = array([[1., 2., 3.], [-4., -5., -6.], [4., 5., 6.], [7., 8., 9.]]).transpose()
    a = bezier(waypoints)
    totalTime = 1.
    degree = a.nbWaypoints-1
    ptsTime = [ (a(float(i) / 10.), float(i) / 10.) for i in range(11)]
    
    bFit = default_fit_3d(degree, ptsTime)    

    fig = plt.figure()
    ax = fig.add_subplot(221, projection="3d")    
    plotBezier(a, ax = ax)
    #~ plotControlPoints(a, ax = ax)
    plotBezier(bFit, ax = ax, color = "g")
    #~ plotControlPoints(bFit, ax = ax, color = "g")
    
    
    bFit = default_fit_3d(degree-1, ptsTime)    
    ax = fig.add_subplot(222, projection="3d")    
    plotBezier(a, ax = ax)
    #~ plotControlPoints(a, ax = ax)
    plotBezier(bFit, ax = ax, color = "g")
    #~ plotControlPoints(bFit, ax = ax, color = "g")
    
    
    bFit = default_fit_3d(degree, ptsTime, zeroInitVel = True)    
    ax = fig.add_subplot(223, projection="3d")    
    plotBezier(a, ax = ax)
    #~ plotControlPoints(a, ax = ax)
    plotBezier(bFit, ax = ax, color = "g")
    #~ plotControlPoints(bFit, ax = ax, color = "g")
    
    bFit = default_fit_3d(degree + 3, ptsTime, zeroInitVel = True)    
    ax = fig.add_subplot(224, projection="3d")    
    plotBezier(a, ax = ax)
    #~ plotControlPoints(a, ax = ax)
    plotBezier(bFit, ax = ax, color = "g")
    #~ plotControlPoints(bFit, ax = ax, color = "g")
    
    plt.show(block=False)
    
