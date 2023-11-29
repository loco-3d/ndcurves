# NDCurves

[![Pipeline status](https://gitlab.laas.fr/loco-3d/ndcurves/badges/master/pipeline.svg)](https://gitlab.laas.fr/loco-3d/ndcurves/commits/master)
[![Coverage report](https://gitlab.laas.fr/loco-3d/ndcurves/badges/master/coverage.svg?job=doc-coverage)](https://gepettoweb.laas.fr/doc/loco-3d/ndcurves/master/coverage/)
[![PyPI version](https://badge.fury.io/py/ndcurves.svg)](https://pypi.org/project/ndcurves)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/loco-3d/ndcurves/master.svg)](https://results.pre-commit.ci/latest/github/loco-3d/ndcurves)


A template-based Library for creating curves of arbitrary order and dimension, eventually subject to derivative constraints. The main use of the library is the creation of end-effector trajectories for legged robots.

To do so, tools are provided to:
 - create **exact** splines of arbitrary order (that pass exactly by an arbitrary number waypoints)
 - constrain initial / end velocities and acceleration for the spline.
 - constrain take-off and landing phases to follow a straight line along a given normal (to avoid undesired collisions between the effector and the contact surface)
 - automatically handle 3d rotation of the effector.
 - create curves in SO3
 - support partial symbolic differentiation of curves. You can represent control points as linear variables, and integrate / differentiate those variable curves. You can also compute the cross product of two curves, which is relevant for centroidal dynamics.

Several type of formulation are provided:
 - Polynomials
 - Bezier
 - Hermite (only cubic hermite for now)

The library is template-based, thus generic:  the curves can be of any dimension, and can be implemented in double or float and can work with kind variables like Vector, Transform, Matrix, ...

## Installation

### Installation through robotpkg

You can install this package through robotpkg. robotpkg is a package manager tailored for robotics softwares.
It greatly simplifies the release of new versions along with the management of their dependencies.
You just need to add the robotpkg apt repository to your sources.list and then use `sudo apt install robotpkg-py3\*-ndcurves`.

If you have never added robotpkg as a softwares repository, please follow first the instructions from 1 to 3; otherwise, go directly to instruction 4.
Those instructions are similar to the installation procedures presented in [http://robotpkg.openrobots.org/debian.html](http://robotpkg.openrobots.org/debian.html).

1. Add robotpkg as source repository to apt:

```bash
sudo tee /etc/apt/sources.list.d/robotpkg.list <<EOF
deb [arch=amd64] http://robotpkg.openrobots.org/wip/packages/debian/pub $(lsb_release -sc) robotpkg
deb [arch=amd64] http://robotpkg.openrobots.org/packages/debian/pub $(lsb_release -sc) robotpkg
EOF
```

2. Register the authentication certificate of robotpkg:

```bash
curl http://robotpkg.openrobots.org/packages/debian/robotpkg.key | sudo apt-key add -
```

3. You need to run at least once apt update to fetch the package descriptions:

```bash
sudo apt-get update
```

4. The installation of nd-curves:

```bash
sudo apt install robotpkg-py3\*-ndcurves # for Python 3

sudo apt install robotpkg-py27-ndcurves # for Python 2
```

Finally you will need to configure your environment variables, e.g.:

```bash
export PATH=/opt/openrobots/bin:$PATH
export PKG_CONFIG_PATH=/opt/openrobots/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/opt/openrobots/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/opt/openrobots/lib/python2.7/site-packages:$PYTHONPATH
```

### Installation from source

#### Dependencies
* [Eigen (version >= 3.2.2)](http://eigen.tuxfamily.org/index.php?title=Main_Page)

#### Additional dependencies for python bindings
* [Boost.Python](http://www.boost.org/doc/libs/1_63_0/libs/python/doc/html/index.html)
* [eigenpy](https://github.com/stack-of-tasks/eigenpy)

To handle this with cmake, use the recursive option to clone the repository.
For instance, using http:
```
git clone --recursive https://github.com/loco-3d/ndcurves $NDCURVES_DIR
```
Where $NDCURVES_DIR is to be replaced to your selected source folder.
The library is header only, so the build only serves to build the tests and python bindings:

```sh
cd $NDCURVES_DIR && mkdir build && cd build
cmake -DBUILD_PYTHON_INTERFACE=ON .. && make && make test
```

Switch the BUILD_PYTHON argument to OFF if you don't want to use the python bindings of the package.

If everything went fine you should obtain the following output:
```sh
100% tests passed, 0 tests failed out of 3
```
#### Optional: Python bindings installation
To install the Python bindings first enable the `BUILD_PYTHON_INTERFACE` option:
```
cmake -DBUILD_PYTHON_INTERFACE=ON ..
```

Then rebuild the library:
```
cd ${NDCURVES_DIR}/build
make && make test
```
To see example of use, you can refer to the [test file](https://github.com/loco-3d/ndcurves/blob/master/python/test/test.py)
which is rather self explanatory:

In spite of an exhaustive documentation, please refer to the C++ documentation, which mostly applies to python.

## Documentation and tutorial

For the C++ / Python detailed documentation, you can consult [this pdf](https://github.com/loco-3d/ndcurves/blob/devel/doc/curves.pdf) available in the doc folder

For a python tutorial, you can refer to the [jupyter notebook](https://github.com/loco-3d/ndcurves/blob/master/python/test/sandbox/test.ipynb).
The [test file](https://github.com/loco-3d/ndcurves/blob/master/python/test/test.py) is more exhaustive and rather self explanatory.

## Citation

If you wish to cite this work please use the bibtex below:

```bib
@software{ndcurves,
  author = {Tonneau, Steve and Chemin, Jason and Fernbach, Pierre and Saurel, Guilhem},
  title = {ndcurves},
  url = {https://github.com/loco-3d/ndcurves},
  year = {2013}
}
```
