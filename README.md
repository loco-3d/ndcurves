Spline
===================


A template-based Library for creating curves of arbitrary order and dimension, eventually subject to derivative constraints. The main use of the library is the creation of end-effector trajectories for legged robots.

To do so, tools are provided to:
> - create **exact** splines of arbitrary order (that pass exactly by an arbitrary number waypoints)
> - constrain initial / end velocities and acceleration for the spline.
> - constrain take-off and landing phases to follow a straight line along a given normal (to avoid undesired collisions between the effector and the contact surface)
> - automatically handle 3d rotation of the effector.

The library is template-based, thus generic:  the curves can be of any dimension, and can be implemented in double, float  ...

While a Bezier curve implementation is provided (limited to degree 3), the main interest 
of this library is to create spline curves of arbitrary order

----------
Example of use for and end-effector trajectory
-------------
The library comes with an helper class to automatically generate end-effector trajectories.
For instance, to create a 2 second long trajectory from the point (0,0,0) to (1,1,0), with a waypoint
at (0.5,0.5,0.5), one can use the following code:

```
typedef std::pair<double, Eigen::Vector3d> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;

// loading helper class namespace
using namespace spline::helpers;

// Create waypoints
waypoints.push_back(std::make_pair(0., Eigen::Vector3d(0,0,0)));
waypoints.push_back(std::make_pair(1., Eigen::Vector3d(0.5,0.5,0.5)));
waypoints.push_back(std::make_pair(2., Eigen::Vector3d(1,1,0)));

exact_cubic_t* eff_traj = effector_spline(waypoints.begin(),waypoints.end());

effector_spline_rotation traj(waypoints.begin(), waypoints.end());

// evaluate spline
(*eff_traj)(0.); // (0,0,0)
(*eff_traj)(2.); // (1,1,0)
```
If rotation of the effector must be considered, the code is almost the same:

```
// initial rotation is 0, end rotation is a rotation by Pi around x axis
quat_t init_rot(0,0,0,1), end_rot(1,0,0,0);
effector_spline_rotation eff_traj_rot(waypoints.begin(),waypoints.end(), init_quat, end_quat);

// evaluate spline
eff_traj_rot(0.); // (0,0,0,0,0,0,1)
eff_traj_rot(1.); // (0.5,0.5,0.5,0.707107,0,0,0.707107) // Pi/2 around x axis
eff_traj_rot(2.); // (0,0,0,1,0,0,0)
```

Additional parameters for the same methods an be used to specify parameters for the take off and
landing phases: height and duration of the phase, and along which normal.
Please refer to the Main.cpp files to see all the unit tests and possibilities offered by the library

Installation
-------------
The library is header only, so you just need to copy the include folder where you need.
Eigen is required for the library to work.

To run the tests, there is a CMakeLists.txt:
```
	cd $SPLINE_DIR && mkdir build && cd build
	cmake .. && make 
	../bin/tests
```

If everything went fine you should obtain the following output:
```
performing tests... 
no errors found 
```
