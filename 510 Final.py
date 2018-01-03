
# coding: utf-8

# In[6]:

import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg as la
import scipy.sparse as sps
get_ipython().magic(u'matplotlib inline')

# Problem 1

def makeA(n):
    return np.diag(2*np.ones(n))+np.diag(-np.ones(n-1),-1)+np.diag(-np.ones(n-1),1)
    
def GMRES(A,b,maxiter=500,precond=False):
    # Since diag(A) = 2*I, M_inverse = 0.5*I
    if precond:
        A = 0.5*A
        b = 0.5*b
        
    m = A.shape[0]
    if maxiter > m:
        maxiter = m-1
        
    Q = np.zeros((m,m),dtype=np.complex128)
    H = np.zeros((m,m),dtype=np.complex128)
    
    bnorm = la.norm(b)
    Q[:,0] = b/bnorm
    
    for n in xrange(1,maxiter):
        print n,'\r',
        # Step n of arnoldi
        Q[:,n] = np.dot(A,Q[:,n-1])
        for i in xrange(n):
            H[i,n-1] = np.dot(Q[:,i].conj(),Q[:,n])
            Q[:,n] -= H[i,n-1]*Q[:,i]
            
        #print Q[:,n]
        H[n,n-1] = la.norm(Q[:,n])
        Q[:,n] /= H[n,n-1]
        
        # Find y to minimize ||H_n*y - la.norm(b)e1||, i.e. r_n
            # H_n is the n+1 by n upper left section of H
        e1 = np.zeros(n+1)
        e1[0] = bnorm
        
        y = la.lstsq(H[:n+1,:n],e1)[0]
        
        r = la.norm(np.dot(H[:n+1,:n],y)-e1)/bnorm
        if r < 1e-8:
            return np.dot(Q[:,:n],y), n
        
    return np.dot(Q[:,:n],y), n

def conjugate_gradient(A,b,maxiter=100,precond=False):
    if precond:
        A = 0.5*A
        b = 0.5*b
        
    m = len(b)
    x = np.ones_like(b)
    r = np.dot(A,x) - b
    p = -r
    n = 0
    done = False
    #maxiter = A.shape[0]
    #betas = np.zeros(maxiter)
    
    while not done and n < maxiter:
        print n,'\r',
        
        Ap = np.dot(A,p)
        rTr = np.dot(r,r)
        
        alpha = rTr/np.dot(p,Ap)
        x = x + alpha*p
        R = r + alpha*Ap
        #r = rprev - alpha*Ap
        beta = np.dot(R,R)/rTr
        #betas[n] = beta
        p = -R + beta*p
        r = R
        
        n += 1
        if np.allclose(r,np.zeros_like(r)):
            done = True
    #plt.plot(np.arange(maxiter),betas)
    #plt.show()
    return x, n

def prob1():
    print "No preconditioning:"
    for n in [50,200,1000]:
        A = makeA(n)
        b = np.ones(n)
        print "\nn = ", n
        
        
        x = la.solve(A,b)
        
        x1, n1 = GMRES(A,b,maxiter=n-1)
        x2, n2 = conjugate_gradient(A,b,maxiter=n)
        
        print "GMRES Error: ", la.norm(np.dot(A,x1)-b)
        print "Conjugate Gradient error: ", la.norm(np.dot(A,x2)-b), '\n'
        
        print "GMRES Iterations: ", n1
        print "Conjugate Gradient Iterations: ", n2
        
    print "\nWith preconditioning:"
    for n in [50,200,1000]:
        A = makeA(n)
        b = np.ones(n)
        print "\nn = ", n
        
        x = la.solve(A,b)
        
        x1, n1 = GMRES(A,b,maxiter=n-1,precond=True)
        x2, n2 = conjugate_gradient(A,b,maxiter=n,precond=True)
        
        print "GMRES Error: ", la.norm(np.dot(A,x1)-b)
        print "Conjugate Gradient error: ", la.norm(np.dot(A,x2)-b), '\n'
        
        print "GMRES Iterations: ", n1
        print "Conjugate Gradient Iterations: ", n2
        
    print "There does not seem to be any improvement with the preconditioner."
        
"""Results:
No preconditioning:

n =  50
GMRES Error:  1.2889840596e-12
Conjugate Gradient error:  3.51126035997e-13 

GMRES Iterations:  25
Conjugate Gradient Iterations:  24

n =  200
GMRES Error:  1.49324462248e-10
Conjugate Gradient error:  1.86765451876e-11 

GMRES Iterations:  100
Conjugate Gradient Iterations:  99

n =  1000
GMRES Error:  3.62179694833e-08
Conjugate Gradient error:  2.14219688697e-09 

GMRES Iterations:  500
Conjugate Gradient Iterations:  500
With preconditioning:

n =  50
GMRES Error:  1.2889840596e-12
Conjugate Gradient error:  3.51126035997e-13 

GMRES Iterations:  25
Conjugate Gradient Iterations:  24

n =  200
GMRES Error:  1.49324462248e-10
Conjugate Gradient error:  1.86765451876e-11 

GMRES Iterations:  100
Conjugate Gradient Iterations:  99

n =  1000
GMRES Error:  3.62179694833e-08
Conjugate Gradient error:  2.14219688697e-09 

GMRES Iterations:  500
Conjugate Gradient Iterations:  500
"""
        
        
# Problem 3

def shiftedQR(B):
    # Since B is symmetric, hessenberg form will be tridiagonal
    w, v = la.eig(B)
    print np.sort(w)[0]
    A, Q = la.hessenberg(B,calc_q=True)
    done = False
    i = 0
    while not done and i < 500:
        mu = A[-1,-1]
        temp = np.diag(A)
        Q, R = la.qr(A-mu*np.eye(A.shape[0]))
        A = np.dot(R,Q)+mu*np.eye(A.shape[0])
        if la.norm(temp-np.diag(A)) < 1e-7:
            done = True
        i += 1
        print i,'\r',
    print i
    return np.sort(np.diag(A))[0]

def prob3():
    n = 1000
    diag = np.array([np.sqrt(i) for i in xrange(1,n+1)])
    sub1 = np.ones(n-1)
    sub100 = np.ones(n-100)
    A = np.diag(diag)+np.diag(sub1,1)+np.diag(sub1,-1)+np.diag(sub100,100)+np.diag(sub100,-100)
    
    #H, Q = la.hessenberg(A,calc_q=True)
    print shiftedQR(A)
    # Note that since only the smallest eigenvalue is needed, 
    
    # Returns -0.91398428 when indexing by zero (i.e. A[0,0] = 0)
    # Returns -0.30096265 when indexing by one (i.e. A[0,0] = 1)
    # Want to be able to check the answer
    #print la.eigh(A,eigvals_only=True)
    
    # Zero indexing result for shiftedQR():
    #(-0.913984284215+0j)
    # (wait for convergence)
    #-0.913984284215
    
    # One indexing result for shiftedQR():
    # (-0.300962645776+0j)
    # wait for convergence
    # -0.300962645776
    
# Eigenvalues of M^-1*A are clustered around 0 and 2. 
# Not exactly close to 1
w = la.eig(0.5*makeA(50))[0]
print la.norm(0.5*makeA(50)-np.eye(50))
plt.scatter(w.real, w.imag)
plt.show()


# In[ ]:



