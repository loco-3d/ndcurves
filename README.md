Curves
===================

[![Pipeline status](https://gepgitlab.laas.fr/loco-3d/curves/badges/master/pipeline.svg)](https://gepgitlab.laas.fr/loco-3d/curves/commits/master)
[![Coverage report](https://gepgitlab.laas.fr/loco-3d/curves/badges/master/coverage.svg?job=doc-coverage)](http://projects.laas.fr/gepetto/doc/loco-3d/curves/master/coverage/)


A template-based Library for creating curves of arbitrary order and dimension, eventually subject to derivative constraints. The main use of the library is the creation of end-effector trajectories for legged robots.

To do so, tools are provided to:
 - create **exact** splines of arbitrary order (that pass exactly by an arbitrary number waypoints)
 - constrain initial / end velocities and acceleration for the spline.
 - constrain take-off and landing phases to follow a straight line along a given normal (to avoid undesired collisions between the effector and the contact surface)
 - automatically handle 3d rotation of the effector.

Several type of formulation are provided:
 - Polynomials
 - Bezier
 - Hermite (only cubic hermite for now)

The library is template-based, thus generic:  the curves can be of any dimension, and can be implemented in double or float and can work with kind variables like Vector, Transform, Matrix, ...

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

This package is available as binary in [robotpkg/wip](http://robotpkg.openrobots.org/robotpkg-wip.html)

## Dependencies
* [Eigen (version >= 3.2.2)](http://eigen.tuxfamily.org/index.php?title=Main_Page)

## Additional dependencies for python bindings
* [Boost.Python](http://www.boost.org/doc/libs/1_63_0/libs/python/doc/html/index.html)
* [eigenpy](https://github.com/stack-of-tasks/eigenpy)

To handle this with cmake, use the recursive option to clone the repository.
For instance, using http:
```
git clone --recursive https://github.com/loco-3d/curves.git
```
The library is header only, so the build only serves to build the tests and python bindings:

```
	cd curves && mkdir build && cd build
	cmake .. && make
	make test
```

If everything went fine you should obtain the following output:
```
performing tests...
no errors found
```
### Optional: Python bindings installation
To install the Python bindings first enable the BUILD_PYTHON_INTERFACE option:
```
cmake -DBUILD_PYTHON_INTERFACE=ON ..
```

Then rebuild the library:
```
cd curves/build
make install
```
To see example of use, you can refer to the [test file](https://github.com/loco-3d/curves/blob/master/python/test/test.py)
which is rather self explanatory:

In spite of an exhaustive documentation, please refer to the C++ documentation, which mostly applies
to python.
