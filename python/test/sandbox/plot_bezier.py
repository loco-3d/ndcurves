from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


def plotControlPoints2D(bez, axes = [0,1], color = "r", ax = None):
    wps = [bez.waypointAtIndex(i) for i in range(bez.nbWaypoints)]
    x = np.array([wp[axes[0]] for wp in wps])
    y = np.array([wp[axes[1]] for wp in wps])
    if ax is not None:
        ax.scatter (x, y, color = color)
    else:        
        plt.scatter(x, y, color = color)
        
def plotBezier2D(bez, axes = [0,1], step = 100., color = "b", showControlPoints = False, ax = None):
    points1 = np.array([(bez(i / step * (bez.max() - bez.min()) + bez.min())[axes[0]], bez(i / step * (bez.max() - bez.min()) + bez.min())[axes[1]]) for i in range(int(step)+1)])
    x = points1[:, 0]
    y = points1[:, 1]
    if ax is None:        
        fig = plt.figure()
        ax = fig.add_subplot(111)
    ax.plot (x, y, color, linewidth=2.0)
    if showControlPoints:
        plotControlPoints2D(bez, color = color, ax = ax)


        

def plotControlPoints(bez,color = "r", ax = None):
    wps = [bez.waypointAtIndex(i) for i in range(bez.nbWaypoints)]
    x = np.array([wp[0] for wp in wps])
    y = np.array([wp[1] for wp in wps])
    z = np.array([wp[2] for wp in wps])
    if ax is None:        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    ax.scatter (x, y, z, color = color)
        
def plotBezier(bez, step = 100., color = "b", linewidth = 2., showControlPoints = False, ax = None):
    points1 = np.array([(bez(i / step * (bez.max() - bez.min()) + bez.min())[0], bez(i / step * (bez.max() - bez.min()) + bez.min())[1], bez(i / step * (bez.max() - bez.min()) + bez.min())[2]) for i in range(int(step)+1)])
    x = points1[:, 0]
    y = points1[:, 1]
    z = points1[:, 2]
    
    if ax is None:        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    ax.plot (x, y, z, color, linewidth=linewidth)
    if showControlPoints:
        plotControlPoints(bez, color = color, ax = ax)




if __name__ == '__main__':
    import eigenpy    
    eigenpy.switchToNumpyArray()
    
    from numpy import array_equal, isclose, matrix, random, array, zeros
    from numpy.linalg import norm

    from curves import (bezier)
    
    waypoints = array([[1., 2., 3.], [-4., -5., -6.], [4., 5., 6.], [7., 8., 9.]]).transpose()
    a = bezier(waypoints)
    
    fig = plt.figure()
    #~ ax = fig.add_subplot(131, projection="3d")
    ax = fig.add_subplot(121)    
    plotBezier2D(a, axes = [1,0], ax = ax)
    plotControlPoints2D(a, axes = [1,0], ax = ax)
    ax = fig.add_subplot(122, projection="3d")    
    plotBezier(a, ax = ax)
    plotControlPoints(a, ax = ax)
    plt.show(block=False)
