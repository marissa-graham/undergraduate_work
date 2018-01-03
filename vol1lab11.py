'''
Solutions file for Volume 1 Lab 11
Call your solutions file 'lab11.py'
'''
import numpy as np
from scipy import linalg as la
from matplotlib import pyplot as plt

def svd_approx(A, k):
    '''
    Calculate the best rank k approximation to A with respect to the induced
    2-norm. Use the SVD.
    Inputs:
        A -- array of shape (m,n)
        k -- positive integer
    Returns:
        Ahat -- best rank k approximation to A obtained via the SVD
    '''
    U, s, Vt = la.svd(A, full_matrices=False)
    S = np.diag(s[:k])
    Ahat = U[:,:k].dot(S).dot(Vt[:k,:])
    return Ahat
    
    
def lowest_rank_approx(A,e):
    '''
    Calculate the lowest rank approximation to A that has error strictly less 
    than e (in induced 2-norm).
    Inputs:
        A -- array of shape (m,n)
        e -- positive floating point number
    Returns:
        Ahat -- the best rank s approximation of A constrained to have error 
                less than e, where s is as small as possible.
    '''
    # 2norm(A-Ahat_s) = sigma[s+1]
    U, s, Vt = la.svd(A, full_matrices=False)
    k = np.where(s < e)[0][0]
    S = np.diag(s[:k])
    return U[:,:k].dot(S).dot(Vt[:k,:])


# In[ ]:



