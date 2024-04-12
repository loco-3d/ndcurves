from curves import bezier, bezierVar
from numpy import array, zeros

__EPS = 1e-6

_zeroMat = array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]).transpose()
_I3 = array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]).transpose()
_zeroVec = array([[0.0, 0.0, 0.0]]).transpose()


def createControlPoint(val):
    if isinstance(val, str):
        return (_I3, _zeroVec)
    else:
        return (_zeroMat, val.reshape([3, 1]))


def createWaypointList(waypoints):
    mat = zeros([3, len(waypoints) * 3])
    vec = zeros([3, len(waypoints)])
    for i, val in enumerate(waypoints):
        mvar, vvar = createControlPoint(val)
        mat[:, i * 3 : i * 3 + 3] = mvar
        vec[:, i : i + 1] = vvar
    return mat, vec


class varBezier:
    # waypoints is a list that contains either
    # 3d arrays (constants), or a string "variable"
    def __init__(self, waypoints=None, time=1.0):
        if waypoints is not None:
            mat, vec = createWaypointList(waypoints)
            self.bezier = bezierVar(mat, vec, time)

    def fromBezier(self, bez):
        var = varBezier()
        var.bezier = bez
        return var

    # splits into n+1 continuous curves, with n the number of time values
    def split(self, times):
        timearray = array(times).reshape([-1, 1])
        subBeziers = self.bezier.split(timearray)
        dim = subBeziers.size
        return [self.fromBezier(subBeziers.at(i)) for i in range(dim)]

    # for each control point of the curve, gives the linear depency
    # of a given variable
    def waypoints(self, varId=-1):
        if varId < 0:
            return self.bezier.waypoints().A, self.bezier.waypoints().b
        assert self.bezier.nbWaypoints > varId
        mat = self.bezier.waypoints().A[:, varId * 3 : varId * 3 + 3]
        vec = self.bezier.waypoints().b[:, varId]
        return mat, vec

    def matrixFromWaypoints(self, varId):
        assert varId >= 0
        mat, vec = self.waypoints(varId)
        resvec = zeros(3)
        for i in range(0, mat.shape[0] / 3, 1):
            resvec += vec[i * 3 : i * 3 + 3]
        return mat.transpose(), resvec

    def toBezier3(self, x):
        wps = []
        for i in range(self.bezier.nbWaypoints):
            mat, vec = self.matrixFromWaypoints(i)
            wps += [mat.dot(x) + vec]
        return bezier(array(wps).transpose(), self.bezier.max())
