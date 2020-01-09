Curves
===================

[![Pipeline status](https://gepgitlab.laas.fr/loco-3d/curves/badges/master/pipeline.svg)](https://gepgitlab.laas.fr/loco-3d/curves/commits/master)
[![Coverage report](https://gepgitlab.laas.fr/loco-3d/curves/badges/master/coverage.svg?job=doc-coverage)](http://projects.laas.fr/gepetto/doc/loco-3d/curves/master/coverage/)


A template-based Library for creating curves of arbitrary order and dimension, eventually subject to derivative constraints. The main use of the library is the creation of end-effector trajectories for legged robots.

To do so, tools are provided to:
> - create **exact** splines of arbitrary order (that pass exactly by an arbitrary number waypoints)
> - constrain initial / end velocities and acceleration for the spline.
> - constrain take-off and landing phases to follow a straight line along a given normal (to avoid undesired collisions between the effector and the contact surface)
> - create curves in SO3
> - support partial symbolic differentiation of curves. You can represent control points as linear variables, and integrate / differentiate those variable curves. You can also compute the cross product of two curves, which is relevant for centroidal dynamics.


The library is template-based, thus generic:  the curves can be of any dimension, and can be implemented in double, float  ...


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
git clone --recursive https://github.com/loco-3d/curves $CURVES_DIR
```
The library is header only, so the build only serves to build the tests and python bindings:

```sh
cd $CURVES_DIR && mkdir build && cd build
cmake .. && make && make test
```

If everything went fine you should obtain the following output:
```sh
100% tests passed, 0 tests failed out of 3
```
### Optional: Python bindings installation
To install the Python bindings, in the CMakeLists.txt file, first enable the BUILD_PYTHON_INTERFACE option:
```cmake
OPTION (BUILD_PYTHON_INTERFACE "Build the python binding" ON)
```

Then rebuild the library:
```sh
cd $CURVES_DIR/build
cmake -DCMAKE_INSTALL_PREFIX=${DEVEL_DIR}/install ..
make install
```
The python bindings should then be accessible through the package centroidal_dynamics.


Documentation and tutorial
-------------

For a python tutorial, you can refer to the [jupyter notebook](https://github.com/loco-3d/curves/blob/devel/python/test/sandbox/test.ipynb) . The [test file](https://github.com/loco-3d/curves/blob/master/python/test/test.py)
is more exhaustive and rather self explanatory.

Please refer to the C++ manual, which mostly applies
to python.
