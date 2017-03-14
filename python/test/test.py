from spline import bezier, spline

from numpy import matrix

waypoints = matrix([[1.,2.,3.],[4.,5.,6.]]).transpose()

#testing bezier curve
a = bezier(waypoints)
a = bezier(waypoints, -1., 3.)
a.min()
a.max()
a(0.4)


#testing spline function
a = spline(waypoints)
a = spline(waypoints, -1., 3.)
a.min()
a.max()
a(0.4)
assert((a.derivate(0.4,0) == a(0.4)).all())
a.derivate(0.4,2)
