import numpy as np
from mayavi import mlab
from scipy import linalg as la
from matplotlib import pyplot as plt
%matplotlib inline

# Problem 2(a)
def weightedlstsq(A,b,C):
    R = la.cholesky(C)
    A, b = np.dot(R,A), np.dot(R,b)
    return la.lstsq(A,b)[0]

# Problem 2(b)
def prob2b():
    def solve_weighted(C):
        data = np.loadtxt('510midtermdata.csv')
        Ra = data[:,0].reshape(27,1)
        Nu = data[:,3].reshape(27,1)
        # Problem: in order to use least squares on an exponential function,
        # we have to take the log of both sides to get something that can be
        # written as a system, i.e. log(Nu) = Ce^(beta*log(Ra)) which gives
        # log(Nu) = log(C) + beta*log(Ra); however, that doesn't make sense,
        # since C is a matrix, not a vector, and experimenting doesn't give
        # a self-consistent result. The only way I can see to 
        # make this a linear system is to treat the problem as Nu = (CRa)^beta
        # and go from there. I understand this will likely lead to me getting
        # only partial credit, but it's the best I can do.

        Ra = np.dot(C,Ra)
        beta, res, rank, s = la.lstsq(np.log(Ra),np.log(Nu))
        plt.scatter(np.log(Nu), np.log(Ra))
        domain = np.log(np.linspace(1,780,27))
        vals = (1/beta*np.dot(C,domain)).reshape(27,)
        #vals = np.exp(beta*domain).reshape(3000,)
        print 'beta=', beta
        plt.plot(domain, vals)
        plt.show()

    print 'Unweighted solution'
    solve_weighted(np.eye(27))
    C = np.diag([1+i for i in xrange(27)])
    print 'Weight lower Ra values more heavily'
    solve_weighted(C)
    print 'Weight higher Ra values more heavily'
    C = np.diag([1+0.01*np.exp(0.1*(27-i)) for i in xrange(27)])
    solve_weighted(C)
    print 'Note that the line is being pulled upward sharply by the last Ra value'

# Problem 3
def prob3():
    def householder(A):
        m, n = A.shape
        R = np.copy(A)
        Q = np.eye(m)
        P = np.eye(n)

        for k in xrange(m):
            # Find i such that norm(R[k:,i]) is maximized
            # This part is not very efficient
            colnorms = np.array([la.norm(R[k:,i]) for i in xrange(k,n)])
            i = np.argmax(colnorms)
            i += k

            # Swap columns of R
            temp = np.copy(R[:,k])
            R[:,k] = R[:,i]
            R[:,i] = temp

            # Swap columns of P
            temp = np.copy(P[:,k])
            P[:,k] = P[:,i]
            P[:,i] = temp

            u = np.copy(R[k:,k])
            u[0] = u[0] + u[0]/np.abs(u[0])*la.norm(u)
            u = u/la.norm(u)
            R[k:,k:] = R[k:,k:] - 2*np.dot(np.outer(u,u), R[k:,k:])
            Q[k:] = Q[k:] - 2*np.dot(np.outer(u,u), Q[k:])
        return P, np.transpose(Q), R

    A = np.random.rand(100,150)
    P, Q, R = householder(A)
    q, r, p = la.qr(A, pivoting=True)

    # Check if Q is right
    print np.allclose(np.abs(Q),np.abs(q))
    # Check if R is right
    print np.allclose(np.abs(R),np.abs(r))
    # Check to make sure the diagonal of R is nonincreasing
    print np.allclose(np.abs(np.diag(R)), np.sort(np.abs(np.diag(R)))[::-1])
    # Check to make sure AP = QR as it should
    print np.allclose(np.dot(A,P),np.dot(Q,R))

# Problem 4
def prob4():
    def randkrank(k):
        U = np.random.rand(100,100)
        U = la.orth(U)

        V = np.random.rand(200,200)
        V = la.orth(V)

        Sdiag = np.sort(np.random.rand(k))[::-1]
        Sdiag = np.append(Sdiag, np.zeros(100-k))
        S = np.hstack((np.diag(Sdiag),np.zeros((100,100))))

        return np.dot(U.T,np.dot(S,V))

    def decomp(k):
        A = randkrank(k)
        P, Q, R = householder(A)
        AP = np.dot(A,P)
        A1 = AP[:,:k]
        A2 = AP[:,k:]

        R1 = R[:k,:k]
        R2 = R[:k,k:]
        R1inv = la.inv(R1)
        R1invR2 = np.dot(R1inv,R2)
        A2prime = np.dot(A1,R1invR2)
        # I'm assuming the error is in A2 versus the constructed
        # version, since it's not really specified.
        print 'Rank ', k, ': ', la.norm(R1invR2), la.norm(A2-A2prime)

    for i in [5,10,20,50,100]:
        decomp(i)
        # We seem to have increasing error as norm(R1invR2) increases
        # Ex: norm(R1inv*R2), norm(A2-A1*R1inv*R2)
        #Rank  5 :   9.53045095598 7.44074624367e-16
        #Rank  10 :  14.9443468043 1.06508136317e-15
        #Rank  20 :  17.4825950006 2.39431277142e-15
        #Rank  50 :  23.873094983  4.96569812573e-15
        #Rank  100 : 28.9752418519 1.09224870153e-14
    
# Problem 5 demonstration
def prob5():
    def condition1(r):
        return 0.5*np.sqrt(2+r**2+1/r**2)
    n = 100
    x = np.linspace(0,10,n+2)[1:-1]
    y = np.zeros_like(x)
    for i in xrange(n):
        y[i] = condition1(x[i])
    plt.plot(x,y)
    plt.title('Condition number as a function of x/y')
    plt.show()