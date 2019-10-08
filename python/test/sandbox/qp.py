import quadprog
from numpy import array, dot, vstack, hstack, asmatrix, identity

from scipy.optimize import linprog


#min (1/2)x' P x + q' x  
#subject to  G x <= h
#subject to  C x  = d
def quadprog_solve_qp(P, q, G=None, h=None, C=None, d=None, verbose = False):
    #~ qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_a = -q
    qp_C = None
    qp_b = None
    meq = 0
    if C is not None:
        if G is not None:
                qp_C = -vstack([C, G]).T
                qp_b = -hstack([d, h])   
        else:
                qp_C = -C.transpose()
                qp_b = -d 
        meq = C.shape[0]
    elif G is not None:  # no equality constraint 
        qp_C = -G.T
        qp_b = -h
    res = quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)
    if verbose:
            return res
    #~ print 'qp status ', res
    return res[0]

#least square form of ||Ax-b||**2 
def to_least_square(A, b):
    print "b ", b
    print "nb ", - dot(A.T, b)
    return dot(A.T, A), - dot(A.T, b)


#min ||Ax-b||**2 
#subject to  G x <= h
#subject to  C x  = d
def solve_least_square(A,b,G=None, h=None, C=None, d=None):
        P, q = to_least_square(A, b)
        #~ P = dot(A.T, A)
        #~ q = -dot(A.T, b).reshape(b.shape[0])
        return quadprog_solve_qp(P, q, G, h, C, d)


#min q' x  
#subject to  G x <= h
#subject to  C x  = d
def solve_lp(q, G=None, h=None, C=None, d=None): 
    res = linprog(q, A_ub=G, b_ub=h, A_eq=C, b_eq=d, bounds=[(-100000.,10000.) for _ in range(q.shape[0])], method='interior-point', callback=None, options={'presolve': True})
    return res    
        

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
