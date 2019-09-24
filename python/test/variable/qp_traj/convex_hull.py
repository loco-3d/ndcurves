import quadprog
from numpy import array, dot, vstack, hstack, asmatrix, identity, cross
from numpy.linalg import norm
import numpy as np

from scipy.spatial import ConvexHull

import matplotlib.pyplot as plt

def genConvexHullLines(points):
    hull = ConvexHull(points)
    lineList = [points[el] for el in hull.vertices] + [points[hull.vertices[0]]]
    lineList = [array(el[:2].tolist() + [0.]) for el in lineList]
    return [[lineList[i],lineList[i+1]] for i in range(len(hull.vertices))], lineList[:-1]
    #now compute lines


def getLineFromSegment(line):
    a = line[0]; b = line[1]; c = a.copy() ; c[2] = 1.
    normal = cross((b-a),(c-a))
    normal /= norm(normal)
    # get inequality
    dire = b - a
    coeff = normal
    rhs = a.dot(normal)
    return (coeff, array([rhs]))




#generate right of the line
def genFromLine(line, num_points, ranges, existing_points = []):
    coeff, rhs = getLineFromSegment(line)
    num_gen = 0
    gen = existing_points + [line[0][:2], line[1][:2]]
    while(len(gen) < num_points):
        pt = array([np.random.uniform(ranges[0][0], ranges[0][1]), np.random.uniform(ranges[1][0], ranges[1][1])])
        if coeff[:2].dot(pt) <= rhs :
            gen += [pt]
    return genConvexHullLines(gen)
        
if __name__ == '__main__':
    genFromLine([array([0.5,0.,0.]),array([0.5,0.5,0.])],5)
    genFromLine([array([0.5,0.,0.]),array([0.5,-0.5,0.])],5)
        
