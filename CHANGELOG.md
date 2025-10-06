# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.1.0] - 1980-01-01

- tooling & packaging updates

## [2.0.0] - 2024-12-05

- :warning: BREAKING: switch from boost smart pointers to std ones
- setup nix
- setup mergify

## [1.4.1] - 2024-04-12

- tests: fix use of np.random in tests
- CMake: enable python bindings by default

## [1.4.0] - 2024-04-12

- Add a SO3 curves which is C2. It's called S03Smooth
- Add some typedef on using 3D curves.
- Fix typos around the package
- Add some tests on polynomials and the new structs.
- update tooling
- update packaging
- add pip CI

## [1.3.1] - 2023-11-29

### Changed

- ⚠️ require CMake 3.10 ⚠️

### Added

- python: define CopyableVisitor, fix SerializableVisitor and use them

### Fixed

- cmake generation on macosx
- Supporting pinocchio installed with cppad (#108)
- fix E721

## [1.3.0] - 2023-07-19

- add stubs support
- fix RPATH for OSX

## [1.2.0] - 2023-05-13

- [python] Enabled support to copy and deepcopy
- [python] Removed deprecated macros in eigenpy
- update for eigenpy v3
- pre-commit update
- sync submodule
- CMake: fetch submodule, set default build type, bump standard

## [1.1.6] - 2023-01-24

- don't build python bindings by default, to be more gentle with PAL robotics buildfarm.

## [1.1.5] - 2022-08-29

- fix allocator
- modernize cmake

## [1.1.4] - 2022-06-29

- remove unary function, deprecated since C++11 and removed in C++17

## [1.1.3] - 2022-05-31

## [1.1.2] - 2022-02-09

- [cmake] Improve Eigen detection/usage
- fix format

## [1.1.1] - 2022-01-06

- fix virtual functions
- fix optional dependency to pinocchio
- primitive now accept initial value\n

## [1.1.0] - 2021-05-03

- make dependency on pinocchio not mandatory
- fix link to Boost::serialization
- use eigen matrix serialization from pinocchio >= 2.6.0 if available
- export -DCURVES_WITH_PINOCCHIO_SUPPORT


## [1.0.0] - 2021-03-18

- ⚠️ library renamed from curves to ndcurves ⚠️
- added arithmetic operations on curves

## [0.5.2] - 2020-09-24

- fix CMake for fedora

## [0.5.1] - 2020-07-27

- fix serialization versioning

## [0.5.0] - 2020-07-25

- Add piecewise::load_from_text_file
- New curves type required by sot-talos-balance
- Use Bezier formulation for cubic hermite
- Check input degree before converting curve to bezier or hermite
- remove cubic and quintic classes


## [0.4.1] - 2020-04-02

Changes since v0.4.0:
- [CMake] fix pinocchio detection
- [CMake] export CURVES_WITH_PINOCCHIO_SUPPORT

## [0.4.0] - 2020-03-30

Changes since v0.3.3:
- [Python] Add pickle support
- Add serialization/curves header

## [0.3.3] - 2020-03-11

Changes since v0.3.2:
- CMake Exports


## [0.3.2] - 2020-02-13

Changes since v0.3.1:
- [Python] Fix binding of Piecewise.curve_at_index
- Replace several critical asserts with exceptions
- Install the header python_definitions.h in include/python
- Add python API to retrieve bezier waypoints as 2D array
- Add C++ and python API to retrieve the translation or rotation curve contained in a SE3 curve
- Correctly specify the corresponding shared_pointer to all python class
- Fix SO3 constructors when t_min == t_max
- Correctly check and raise error when trying to use polynomial constructors from boundary condition when t_min == t_max
- Reworking of the exposition of the abstract class in Python (fix  https://gepgitlab.laas.fr/loco-3d/curves/issues/32)
- Correctly register the shared_pointer of the base abstract classes in boost::Python


## [0.3.1] - 2020-02-13

Changes since v0.3.0:
- [CMake] add INSTALL_PYTHON_INTERFACE_ONLY option
- Update README

## [0.3.0] - 2020-01-10

Changes since v0.2.0:
- [CMake] update minimal eigenpy version
- Add operator == and != for all curves
- Add methods isApprox() and isEquivalent() to all curves
- Fix optional dependency to pinocchio for python bindings
- [Python][Tests] remove all unecessary reshape(-1,1) in python
- Rework the class piecewise_curve to make it generic and remove the need to specify the type of curve used as a template argument, can now mix any kind of curves inside
- Rework the methods convert_to_X_from_Y to convert_to_X and remove the need to specify the input type as template argument
- Add compute_derivate_ptr() method to curve_abc and implement it in all child classes
- Add degree() method to curve_abc and implement it in all child class
- Factorize the most commonly used typedef with template argument in fwd.h
- [Python] correctly define eigenpy matrix type for point3 and point6
- [CMake] fix hardcoded path
- [CMake] fix install path of optimization files
- [Tests][Python] use Quaternion.isApprox to test equality instead of ==
- Add SE3 with pinocchio
- Add conversion functions from piecewise curve to python bindings
- Optimization
- Fix all compiler warnings
- Export plot

## [0.2.0] - 2019-10-04

- Initial release

[Unreleased]: https://github.com/loco-3d/ndcurves/compare/v2.1.0...HEAD
[2.1.0]: https://github.com/loco-3d/ndcurves/compare/v2.0.0...v2.1.0
[2.0.0]: https://github.com/loco-3d/ndcurves/compare/v1.4.1...v2.0.0
[1.4.1]: https://github.com/loco-3d/ndcurves/compare/v1.4.0...v1.4.1
[1.4.0]: https://github.com/loco-3d/ndcurves/compare/v1.3.1...v1.4.0
[1.3.1]: https://github.com/loco-3d/ndcurves/compare/v1.3.0...v1.3.1
[1.3.0]: https://github.com/loco-3d/ndcurves/compare/v1.2.0...v1.3.0
[1.2.0]: https://github.com/loco-3d/ndcurves/compare/v1.1.6...v1.2.0
[1.1.6]: https://github.com/loco-3d/ndcurves/compare/v1.1.5...v1.1.6
[1.1.5]: https://github.com/loco-3d/ndcurves/compare/v1.1.4...v1.1.5
[1.1.4]: https://github.com/loco-3d/ndcurves/compare/v1.1.3...v1.1.4
[1.1.3]: https://github.com/loco-3d/ndcurves/compare/v1.1.2...v1.1.3
[1.1.2]: https://github.com/loco-3d/ndcurves/compare/v1.1.1...v1.1.2
[1.1.1]: https://github.com/loco-3d/ndcurves/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/loco-3d/ndcurves/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/loco-3d/ndcurves/compare/v0.5.2...v1.0.0
[0.5.2]: https://github.com/loco-3d/ndcurves/compare/v0.5.1...v0.5.2
[0.5.1]: https://github.com/loco-3d/ndcurves/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/loco-3d/ndcurves/compare/v0.4.1...v0.5.0
[0.4.1]: https://github.com/loco-3d/ndcurves/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/loco-3d/ndcurves/compare/v0.3.3...v0.4.0
[0.3.3]: https://github.com/loco-3d/ndcurves/compare/v0.3.2...v0.3.3
[0.3.2]: https://github.com/loco-3d/ndcurves/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/loco-3d/ndcurves/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/loco-3d/ndcurves/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/loco-3d/ndcurves/releases/tag/v0.2.0
