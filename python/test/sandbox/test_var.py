from curves import bezierVar
from numpy import array

__EPS = 1e-6


waypointsA = array(
    [
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
    ]
).transpose()
waypointsb = array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]).transpose()

# testing bezier curve
a = bezierVar(waypointsA, waypointsb, 0.0, 3.0)

subBeziers = a.split(array([[0.2, 0.4]]).transpose())
assert subBeziers.size == 3
assert subBeziers.at(0).max() - 0.2 <= __EPS
assert subBeziers.at(1).max() - subBeziers.at(1).min() - 0.2 <= __EPS
assert subBeziers.at(2).max() - subBeziers.at(2).min() - 2.6 <= __EPS
