from spline import bezier, spline, exact_cubic, spline_constraints, spline_constraints, spline_deriv_constraint

from numpy import matrix

waypoints = matrix([[1.,2.,3.],[4.,5.,6.]]).transpose()
time_waypoints = matrix([0.,1.])

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

#testing exact_cubic function
a = exact_cubic(waypoints, time_waypoints)
a.min()
a.max()
a(0.4)
assert((a.derivate(0.4,0) == a(0.4)).all())
a.derivate(0.4,2)

#testing spline_deriv_constraints
c = spline_constraints();
#~ c.init_vel; // TODO: error in reading DATA at the time ...
#~ c.end_vel;
c.init_acc = matrix([0.,1.,1.]);
#~ c.end_acc;

a = spline_deriv_constraint (waypoints, time_waypoints)
a = spline_deriv_constraint (waypoints, time_waypoints, c)
