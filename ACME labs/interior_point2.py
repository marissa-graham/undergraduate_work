from cvxopt import matrix
from cvxopt import solvers
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy import sparse as spar
%matplotlib inline

# Problem 2
def problem2(n):
    def laplacian(n):
        """
        Construct the discrete Dirichlet energy matrix H for an n x n grid.
        Inputs:
            n -- side length of grid
        Returns:
            dense array of shape n^2 x n^2
        """
        data = -1*np.ones((5, n**2))
        data[2,:] = 4
        data[1, n-1::n] = 0
        data[3, ::n] = 0
        diags = np.array([-n, -1., 0, 1., n])
        return spar.spdiags(data, diags, n**2, n**2).todense()

    #create the tent pole configuration
    L = np.zeros((n,n))
    L[n/2-1:n/2+1,n/2-1:n/2+1] = .5
    m = [n/6-1, n/6, int(5*(n/6.))-1, int(5*(n/6.))]
    mask1, mask2 = np.meshgrid(m, m)
    L[mask1, mask2] = .3
    L = -1*L.ravel()

    #initial guess
    x = np.ones((n,n))
    x = x.ravel()
    y = np.ones(n**2)
    l = np.ones(n**2)

    A = -1*np.eye(n**2)
    H = laplacian(n)
    c = -1.0/(n-1)**2 * np.ones(n**2)
    
    H = matrix(H)
    A = matrix(A)
    c = matrix(c)
    L = matrix(L)
    
    sol = solvers.qp(H,c,A,L)
    z = np.array(sol['x']).reshape(n,n)
    
    #plot the solution
    dom = np.arange(n)
    X, Y = np.meshgrid(dom, dom)
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.plot_surface(X, Y, z, rstride=1, cstride=1, color='r')
    plt.show()
    
def problem3():
    data = np.loadtxt('portfolio.txt')
    mu = np.mean(data, axis=0)
    
    n = 8
    noyrs = data[:,1:].T
    Q = matrix(np.cov(noyrs))
    c = matrix(np.zeros(n))
    G = matrix(-np.eye(n))
    h = matrix(np.zeros(n))
    A = matrix(np.vstack((np.ones(n),mu[1:])))
    b = matrix(np.array([[1],[1.13]]))
    
    sol = solvers.qp(Q,c,G,h,A,b)
    return np.array(sol['x'])
    
problem3()