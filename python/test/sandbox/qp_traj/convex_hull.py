import numpy as np
from numpy import array, cross
from numpy.linalg import norm
from scipy.spatial import ConvexHull


def genConvexHullLines(points):
    hull = ConvexHull(points)
    lineList = [points[el] for el in hull.vertices] + [points[hull.vertices[0]]]
    lineList = [array([*el[:2].tolist(), 0.0]) for el in lineList]
    return [
        [lineList[i], lineList[i + 1]] for i in range(len(hull.vertices))
    ], lineList[:-1]
    # now compute lines


def getLineFromSegment(line):
    a = line[0]
    b = line[1]
    c = a.copy()
    c[2] = 1.0
    normal = cross((b - a), (c - a))
    normal /= norm(normal)
    # get inequality
    coeff = normal
    rhs = a.dot(normal)
    return (coeff, array([rhs]))


def genFromLine(line, num_points, ranges, existing_points=[]):
    """generate right of the line"""
    coeff, rhs = getLineFromSegment(line)
    gen = [*existing_points, line[0][:2], line[1][:2]]
    rng = np.random.default_rng()
    while len(gen) < num_points:
        pt = array(
            [
                rng.uniform(ranges[0][0], ranges[0][1]),
                rng.uniform(ranges[1][0], ranges[1][1]),
            ]
        )
        if coeff[:2].dot(pt) <= rhs:
            gen += [pt]
    return genConvexHullLines(gen)


if __name__ == "__main__":
    genFromLine([array([0.5, 0.0, 0.0]), array([0.5, 0.5, 0.0])], 5)
    genFromLine([array([0.5, 0.0, 0.0]), array([0.5, -0.5, 0.0])], 5)
