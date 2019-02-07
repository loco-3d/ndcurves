import quadprog
from numpy import array, dot, vstack, hstack, asmatrix, identity

#min (1/2)x' P x + q' x  
#subject to  G x <= h
#subject to  C x  = d
def quadprog_solve_qp(P, q, G=None, h=None, C=None, d=None):
    qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_a = -q
    if C is not None:
        if G is not None:
                qp_C = -vstack([C, G]).T
                qp_b = -hstack([d, h])   
        else:
                qp_C = -C.transpose()
                qp_b = -d 
        meq = C.shape[0]
    else:  # no equality constraint 
        qp_C = -G.T
        qp_b = -h
        meq = 0
    return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]


#min ||Ax-b||**2 
#subject to  G x <= h
#subject to  C x  = d
def solve_least_square(A,b,G=None, h=None, C=None, d=None):
        P = dot(A.T, A)
        #~ q = 2*dot(b, A).reshape(b.shape[0])
        q = 2*dot(b, A)
        return quadprog_solve_qp(P, q, G, h)


#min q' x  
#subject to  G x <= h
#subject to  C x  = d
def solve_lp(q, G=None, h=None, C=None, d=None):
    qp_G = identity(q.shape[0]) * 0.000001
    qp_a = -q
    if C is not None:
        if G is not None:
                qp_C = -vstack([C, G]).T
                qp_b = -hstack([d, h])  
        else:
                qp_C = -C.transpose()
                qp_b = -d 
        meq = C.shape[0]
    else:  # no equality constraint 
        qp_C = -G.T
        qp_b = -h
        meq = 0
    return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]
        

if __name__ == '__main__':
        
        from numpy.linalg import norm
        
        A = array([[1., 2., 0.], [-8., 3., 2.], [0., 1., 1.]])
        b = array([3., 2., 3.])
        P = dot(A.T, A)
        q = 2*dot(b, A).reshape((3,))
        G = array([[1., 2., 1.], [2., 0., 1.], [-1., 2., -1.]])
        h = array([3., 2., -2.]).reshape((3,))

        res2 = solve_least_square(A, b, G, h)
        res1 =  quadprog_solve_qp(P, q, G, h)
        print res1
        print res2
