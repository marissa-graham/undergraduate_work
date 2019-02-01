import numpy as np
from scipy import linalg as la
from matplotlib import pyplot as plt

def startingPoint(A, b, c):
    '''
    Calculate an initial guess to the solution of the
    linear program min c^T x, Ax = b, x>=0.
    Inputs:
        A -- array of shape (m,n) with linearly independent rows
        b -- array of length m
        c -- array of length n
    Returns:
        x -- array of length n
        lam -- array of length m
        s -- array of length n
    Ref: Nocedal and Wright, p. 410
    '''
    # first calculate x, lam, s of minimal norm satisfying the primal and dual constraints
    B = la.inv(A.dot(A.T))
    x = A.T.dot(B.dot(b))
    lam = B.dot(A.dot(c))
    s = c - A.T.dot(lam)

    # perturb x and s so they are nonnegative
    dx = max((-3./2)*x.min(), 0)
    ds = max((-3./2)*s.min(), 0)
    x += dx*np.ones(x.shape)
    s += ds*np.ones(s.shape)

    # perturb x and s so they are not too close to zero, not too dissimilar
    dx = .5*(x*s).sum()/s.sum()
    ds = .5*(x*s).sum()/x.sum()
    x += dx*np.ones(x.shape)
    s += ds*np.ones(s.shape)

    return x, lam, s

def interiorPoint(A, b, c, niter=20, verbose=False, starting_point=None, pts=False):
    '''
    Solve the linear programming problem min c^T x, Ax = b, x>=0
    using an Interior Point method. This code is not optimized, but
    forms the basis for a common practical approach known as the
    Predictor-Corrector Algorithm.
    Inputs:
        A -- array of shape (m,n) with linearly independent rows
        b -- array of length m
        c -- array of length n
        niter -- positive integer giving the number of iterations
        starting_point -- tuple of arrays giving the initial values for x, l, and s.
                          if unspecified, the function startingPoint is used.
    Returns:
        x -- the optimal point
        val -- the minimum value of the objective function
        (pts -- list of points traced by the algorithm, returned if pts=True)
    Ref: Nocedal and Wright, p. 411
    '''
    pts = []
    # initialize variables
    m,n = A.shape
    if starting_point:
        x, l, s = starting_point
    else:
        x,l,s = startingPoint(A,b,c)
    pts.append(x)
    N = np.zeros((n+m+n, n+m+n))
    N[:n, n:n+m] = A.T
    N[:n, n+m:] = np.eye(n)
    N[n:n+m, :n] = A
    sol = np.empty(n+m+n)
    newsol = np.empty(n+m+n)
    for k in xrange(niter):
        # finish initializing parts of the step equation
        N[n+m:, :n] = np.diag(s)
        N[n+m:, n+m:] = np.diag(x)
        r_c = (A.T).dot(l)+s-c
        r_b = A.dot(x)-b
        rhs = np.hstack((-r_c.ravel(), -r_b.ravel(), -x*s))

        # solve dx_aff, dl_aff, ds_aff using LU decomposition
        lu_piv = la.lu_factor(N)
        sol[:] = la.lu_solve(lu_piv, rhs)
        
        dx_aff = sol[:n]
        dl_aff = sol[n:n+m]
        ds_aff = sol[n+m:]

        mu = np.dot(x,s)/n
        
        # calculate a_p, a_d, mu_aff, sig
        mask1 = dx_aff < 0
        if np.any(mask1):
            a_p = min(1, (-x/dx_aff)[mask1].min())
        mask2 = ds_aff < 0
        if np.any(mask1):
            a_d = min(1, (-s/ds_aff)[mask2].min())
        mu_aff = 1/n*(x+a_p*dx_aff).T*(s+a_d*ds_aff)
        sig = (mu_aff/mu)**3
        
        # calculate dx, dl, ds
        newrhs = np.hstack((-r_c.ravel(), -r_b.ravel(), -x*s-dx_aff*ds_aff+mu*sig))
        newsol[:] = la.lu_solve(lu_piv,newrhs)
        dx = newsol[:n]
        dl = newsol[n:n+m]
        ds = newsol[n+m:]

        # calculate ap, ad
        mask3 = dx < 0
        if np.any(mask3):
            bp = min(1, (-x/dx)[mask3].min())
        mask4 = ds < 0
        if np.any(mask4):
            bd = min(1, (-s/ds)[mask4].min())
        ap = min(1,0.95*bp)
        ad = min(1,0.95*bd)

        # step to new point
        x = x + ap*dx
        pts.append(x)
        l, s = l+ad*dl, s+ad*ds
        if verbose:
            print 'x: ', x
            print 'mu: ', mu

    if pts:
        return pts
    else:
        return x, (c*x).sum()


def leastAbsoluteDeviations():
    '''
    For this problem, write an explanation on paper of the following:
        Explicitly write and explain: c, y, x, A
        Write a small explanation of Least Absolute Deviations and how interior point provides a solution.
    '''
    
    data = np.loadtxt('simdata.txt')
    m = data.shape[0]
    n = data.shape[1] - 1
    c = np.zeros(3*m + 2*(n + 1))
    c[:m] = 1
    y = np.empty(2*m)
    y[::2] = -data[:, 0]
    y[1::2] = data[:, 0]
    x = data[:, 1:]
    A = np.ones((2*m, 3*m + 2*(n + 1)))
    A[::2, :m] = np.eye(m)
    A[1::2, :m] = np.eye(m)
    A[::2, m:m+n] = -x
    A[1::2, m:m+n] = x
    A[::2, m+n:m+2*n] = x
    A[1::2, m+n:m+2*n] = -x
    A[::2, m+2*n] = -1
    A[1::2, m+2*n+1] = -1
    A[:, m+2*n+2:] = -np.eye(2*m, 2*m)
    
    sol = interiorPoint(A, y, c, niter=10)[-1]
    beta = (sol[m:m+n] - sol[m+n:m+2*n])[0]
    b = sol[m+2*n] - sol[m+2*n+1]
    
    dom = np.linspace(0,10,2)
    plt.scatter(data[:,1], data[:,0])
    plt.plot(dom, beta*dom+b)
    plt.show()
    print 'Beta:', beta
    print 'b:', b