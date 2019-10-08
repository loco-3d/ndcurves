
from numpy import array_equal, isclose, matrix, random, zeros, array
from numpy.linalg import norm
import numpy as np

import eigenpy
eigenpy.switchToNumpyArray()

from curves import (bezier)
from curves.optimization import(constraint_flag, integral_cost_flag, setup_control_points, problemData, problemDefinition)

from qp import to_least_square, quadprog_solve_qp

import ipdb; 


np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)})

waypoints = array([[1., 2., 3.], [-4., -5., -6.], [4., 5., 6.], [7., 8., 9.]]).transpose()
a = bezier(waypoints)
totalTime = 1.
degree = 3
ptsTime = [ (a(float(i) / 10.), float(i) / 10.) for i in range(11)]
    

#points given as pairs (pt, t)
#does not try to fit anything
#~ def default_fit_3d(degree, ptsTime, totalTime =1.):
pD = problemDefinition()
pD.flag = constraint_flag.NONE
pD.totalTime = totalTime
pD.degree = degree - 1
problem = setup_control_points(pD)
bezierLinear = problem.bezier()

nVar = problem.numVariables * 3 #dimension 3
#~ nVar = 0
A, b = zeros((nVar,nVar)), zeros(nVar)
#~ ipdb.set_trace()
for (pt, time) in ptsTime:
    assert time <= totalTime, "total time inferior to sampling"
    exprAtT = bezierLinear(time)
    nA, nb = to_least_square(exprAtT.A, exprAtT.b + pt)
    A += nA
    b += nb
    
res = quadprog_solve_qp(A, b)

bFit = bezierLinear.evaluate(res.reshape((-1,1)) ) 

from plot_bezier import *

fig = plt.figure()
#~ ax = fig.add_subplot(131, projection="3d")
ax = fig.add_subplot(111, projection="3d")    
plotBezier(a, ax = ax)
plotControlPoints(a, ax = ax)
plotBezier(bFit, ax = ax, color = "g")
plotControlPoints(bFit, ax = ax, color = "g")
plt.show(block=False)

#~ return A, b
    
