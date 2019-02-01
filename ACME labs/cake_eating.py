import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats as st
from scipy import linalg as la

def eatCake(beta,N,W_max=1.,finite=True,T=None,plot=True):
    """
    
    Parameters
    ----------
    beta : float
        The discount factor, 0 < beta < 1
    N : integer
        The number of discrete cake values to consider
    W_max : float
        The original size of the cake, default 1.0
    finite : bool
        Whether or not to consider a finite horizon
    T : integer
        The number of time periods, None if finite is False
    plot : bool
        Whether or not to plot the results.
        
    Returns
    -------
    psi : ndarray of shape (N,T+1)
        The value function
    V : ndarray of shape (N,T+2)
        The policy function
    
    """
    W = np.linspace(0,W_max,N)
    f = lambda i,j:np.sqrt(W[i]-W[j]) if i >= j else -np.inf
    u = np.array([[f(i,j) for i in xrange(N)] for j in xrange(N)])

    if finite:
        V = np.zeros((N,T+2))
        psi = np.empty((N,T+1))
        # Iteratively set the proper values of the columns of V and psi in reverse order.
        # np.max and np.argmax will be useful here, using the proper axis as a keyword.
        # For psi you will also have to use the array W.
        nrange = np.arange(N)
        for i in xrange(T+1,0,-1):
            # We compute u(W_T-W_T+1)+beta*V_T+1(W_T+1) for all W_T and W_T+1
            Wtemp = u + beta*V[:,i].reshape((1,N))
            argmaxs = np.argmax(Wtemp,axis=1)
            V[:,i-1] = Wtemp[nrange,argmaxs].reshape((1,N))
            psi[:,i-1] = W[argmaxs]
        
        if plot:
            x = np.arange(0,N)
            y = np.arange(0,T+2)
            X, Y = np.meshgrid(x,y)
            fig1 = plt.figure()
            ax1 = Axes3D(fig1)
            ax1.plot_surface(W[X],Y,np.transpose(V),cmap=cm.coolwarm)
            plt.show()
            
            fig2 = plt.figure()
            ax2 = Axes3D(fig2)
            y = np.arange(0,T+1)
            X,Y = np.meshgrid(x,y)
            ax2.plot_surface(W[X],Y,np.transpose(psi),cmap=cm.coolwarm)
            plt.show()
    
    else:
    
        Vprev = np.zeros(N)
        psi = np.zeros(N)
        distance = 1
        nrange = np.arange(N)
        while distance >= 10e-9:
            Wtemp = u + beta*Vprev.reshape((1,N))
            argmaxs = np.argmax(Wtemp,axis=1)
            V = Wtemp[nrange,argmaxs].reshape((1,N))
            psi = W[argmaxs]
            distance = la.norm(V-Vprev)
            Vprev = V
            i += 1
        V = Vprev.reshape(N,)
        
        if plot:
            plt.plot(nrange,V)
            plt.show()
            plt.plot(nrange,psi)
            plt.show()
            
    return V, psi

def prob1():
    eatCake(.9,100,T=10)
    
def prob2():
    eatCake(.9,100,T=1000)

def prob3():
    eatCake(.9,100,finite=False)