from spline import bezier

from numpy import matrix

a = bezier(matrix([[1.,2.,3.],[4.,5.,6.]]))
a = bezier(matrix([[1.,2.,3.],[4.,5.,6.]]), -1., 3.)
