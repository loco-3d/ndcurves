from spline import bezier

from numpy import matrix

#testing bezier curve
a = bezier(matrix([[1.,2.,3.],[4.,5.,6.]]))
a = bezier(matrix([[1.,2.,3.],[4.,5.,6.]]), -1., 3.)
a.min()
a.max()
a(0.4)
