from numpy import array

import eigenpy
from curves import bezierVar

__EPS = 1e-6

eigenpy.switchToNumpyArray()

waypointsA = array([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.], [1., 2., 3.], [4., 5., 6.], [7., 8., 9.]]).transpose()
waypointsb = array([[1., 2., 3.], [1., 2., 3.]]).transpose()

# testing bezier curve
a = bezierVar(waypointsA, waypointsb, 3.)

subBeziers = a.split(array([[0.2, 0.4]]).transpose())
assert (subBeziers.size == 3)
assert (subBeziers.at(0).max() - 0.2 <= __EPS)
assert (subBeziers.at(1).max() - 0.2 <= __EPS)
assert (subBeziers.at(2).max() - 2.6 <= __EPS)
