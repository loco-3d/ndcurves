from spline import bezier, spline, exact_cubic, curve_constraints, spline_deriv_constraint

from numpy import matrix

waypoints = matrix([[1.,2.,3.],[4.,5.,6.]]).transpose()
time_waypoints = matrix([0.,1.])

#testing bezier curve
a = bezier(waypoints, -1., 3.)

a = bezier(waypoints)
assert(a.degree == a.nbWaypoints -1)
a.min()
a.max()
a(0.4)
assert((a.derivate(0.4,0) == a(0.4)).all())
a.derivate(0.4,2)
a = a.compute_derivate(100)

prim = a.compute_primitive(1)

for i in range(10):
	t = float(i) / 10.
	assert(a(t) == prim.derivate(t,1)).all()
assert(prim(0) == matrix([0.,0.,0.])).all()

prim = a.compute_primitive(2)
for i in range(10):
	t = float(i) / 10.
	assert(a(t) == prim.derivate(t,2)).all()
assert(prim(0) == matrix([0.,0.,0.])).all()



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
c = curve_constraints();
c.init_vel; 
c.end_vel;
c.init_acc;
c.end_acc;


c.init_vel = matrix([0.,1.,1.]);
c.end_vel  = matrix([0.,1.,1.]);
c.init_acc = matrix([0.,1.,1.]);
c.end_acc  = matrix([0.,1.,1.]);

a = spline_deriv_constraint (waypoints, time_waypoints)
a = spline_deriv_constraint (waypoints, time_waypoints, c)
