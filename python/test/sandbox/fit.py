
from numpy import array_equal, isclose, matrix, random, zeros, array, identity
from numpy.linalg import norm
import numpy as np
np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)})

import eigenpy
eigenpy.switchToNumpyArray()

from curves import (bezier, curve_constraints)
from curves.optimization import(constraint_flag, integral_cost_flag, setup_control_points, problem_data, problem_definition)

from qp import to_least_square, quadprog_solve_qp


  

#points given as pairs (pt, t)
#does not try to fit anything
def default_fit(degree, ptsTime, dim = 3, totalTime =1., constraintFlag = constraint_flag.NONE, curveConstraints = None):
    pD = problem_definition(dim)
    #set initial and goal positions to the ones of the list (might not be used if constraints do not correspond)
    pD.start = array([ptsTime[0][0]]).T
    pD.end   = array([ptsTime[-1][0]]).T
    #assign curve constraints
    pD.totalTime = totalTime
    pD.degree = degree
    pD.flag = constraintFlag
    problem = setup_control_points(pD)
    bezierLinear = problem.bezier()
    if curveConstraints is not None:
        pD.curveConstraints = curveConstraints
    pD.totalTime = totalTime
    pD.degree = degree
    problem = setup_control_points(pD)
    bezierLinear = problem.bezier()

    allsEvals = [(bezierLinear(time), pt) for (pt,time) in ptsTime]
    allLeastSquares = [to_least_square(el.A, el.b + pt) for (el, pt) in  allsEvals]
    Ab = [sum(x) for x in zip(*allLeastSquares)]
    
    res = quadprog_solve_qp(Ab[0] + identity(problem.numVariables * dim) * 0.0001, Ab[1])
    return bezierLinear.evaluate(res.reshape((-1,1)) ) 


#~ return A, b

if __name__ == '__main__':

    from plot_bezier import *
    
    waypoints = array([[1., 2., 3., 0.], [-4., -5., -6., 0.], [4., 5., 6., 0.], [7., 8., 9., 0.]]).transpose()
    a = bezier(waypoints)
    totalTime = 1.
    degree = a.nbWaypoints-1
    ptsTime = [ (a(float(i) / 10.), float(i) / 10.) for i in range(11)]
    dim = waypoints.shape[0]
    
    c = curve_constraints(dim)
    c.init_vel = array([[0., 0., 0.,0.]]).T
    
    def plot_fit_and_original(ax, plotControlPoints = False):
        plotBezier(a, ax = ax)
        plotBezier(bFit, ax = ax, color = "g")
        if plotControlPoints:
            plotControlPoints(a, ax = ax)
            plotControlPoints(bFit, ax = ax, color = "g")
    
    fig = plt.figure()
    
    bFit = default_fit(degree, ptsTime, dim = dim)    
    ax = fig.add_subplot(221, projection="3d")    
    plot_fit_and_original(ax)
    
    bFit = default_fit(degree-1, ptsTime, dim = dim)    
    ax = fig.add_subplot(222, projection="3d")     
    plot_fit_and_original(ax)
            
    
    bFit = default_fit(degree-1, ptsTime, dim = dim, constraintFlag = constraint_flag.INIT_VEL | constraint_flag.INIT_POS, curveConstraints = c)    
    ax = fig.add_subplot(223, projection="3d")    
    plot_fit_and_original(ax)
    
    
    bFit = default_fit(degree + 30, ptsTime, dim = dim, constraintFlag = constraint_flag.END_POS |  constraint_flag.INIT_VEL | constraint_flag.INIT_POS, curveConstraints = c)    
    ax = fig.add_subplot(224, projection="3d")     
    plot_fit_and_original(ax)
    
    plt.show(block=False)
    
