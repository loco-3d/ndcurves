from numpy import array, vstack, zeros

from .varBezier import varBezier

__EPS = 1e-6


# ### helpers for stacking matrices ####
def concat(m1, m2):
    if m1 is None:
        return m2
    return vstack([m1, m2]).reshape([-1, m2.shape[-1]])


def concatvec(m1, m2):
    if m1 is None:
        return m2
    return array(m1.tolist() + m2.tolist())


# ### helpers for stacking matrices ####


# constraint to lie between 2 extremities of a segment
def segmentConstraint(
    varBez, a, b, wpIndex, constraintVarIndex, totalAddVarConstraints
):
    mat, vec = varBez.matrixFromWaypoints(wpIndex)
    vec = b - vec
    resmat = zeros([mat.shape[0], mat.shape[1] + totalAddVarConstraints])
    resmat[:, : mat.shape[1]] = mat
    resmat[:, mat.shape[1] + constraintVarIndex - 1] = b - a
    return (resmat, vec)


def lineConstraint(varBez, C, d, totalAddVarConstraints):
    """constraint to lie one side of a line"""
    resmat = None
    resvec = None
    for wpIndex in range(varBez.bezier.nbWaypoints):
        mat, vec = varBez.matrixFromWaypoints(wpIndex)
        mat = C.dot(mat)
        vec = d - C.dot(vec)
        resmat = concat(resmat, mat)
        resvec = concatvec(resvec, vec)
    augmented = zeros([resmat.shape[0], resmat.shape[1] + totalAddVarConstraints])
    augmented[:, : resmat.shape[1]] = resmat
    return (augmented, resvec)


# #### cost function integrals #####
def compute_primitive(wps):
    """compute bezier primitive of variable waypoints given these waypoints"""
    new_degree = len(wps) - 1 + 1
    current_sum = [zeros(wps[0][0].shape), 0]
    n_wp = [(current_sum[0], 0.0)]
    for wp in wps:
        current_sum[0] = current_sum[0] + wp[0]
        current_sum[1] = current_sum[1] + wp[1]
        n_wp += [(current_sum[0] / new_degree, current_sum[1] / new_degree)]
    return n_wp


def accelerationcost(bezVar):
    """cost function for analytical integral cost of squared acceleration"""
    acc = bezVar.bezier.compute_derivate(2)
    hacc = varBezier()
    hacc.bezier = acc
    wps_i = [[], [], []]
    for i in range(3):  # integrate for each axis
        for j in range(acc.nbWaypoints):
            A_i = hacc.matrixFromWaypoints(j)[0][i, :].reshape([1, -1])
            b_i = hacc.matrixFromWaypoints(j)[1][i]
            wps_i[i] += [(A_i.transpose().dot(A_i), b_i * b_i)]
    # now integrate each bezier curve
    for i, wps in enumerate(wps_i):
        wps_i[i] = compute_primitive(wps)

    resmat = wps_i[0][-1][0]
    resvec = wps_i[0][-1][1]
    for i in range(1, 3):
        resmat = resmat + wps_i[i][-1][0]
        resvec = resvec + wps_i[i][-1][1]
    return (resmat, resvec)


# #### cost function integrals #####
