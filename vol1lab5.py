import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg as la

def plot_old_new(old, new):
    '''
    This function is for you to look at your transformations. 
    Do NOT call it in your final submission of the lab.
    
    Inputs:
    new -- a (2,n) numpy array containing x-coordinates on the first row and 
           y-coordinates on the second row.
    old -- a (2,n) numpy array containing x-coordinates on the first row and 
           y-coordinates on the second row.
    '''
    plt.subplot(2, 1, 1)
    plt.scatter(old[0], old[1])
    plt.axis('equal')
    plt.subplot(2, 1, 2)
    plt.scatter(new[0], new[1])
    plt.show()

######################################################### Problem 1
def dilation(pts, stretch):
    '''
    Do NOT plot
    Hint: Use helper function plot_old_new to view your transformations
    Inputs: 
    pts     -- (2,n) array where the first row represents the x-coordinates and the 
               second row the y-coordinates
    stretch -- one-dimensional numpy array with two entries for dilation factor in 
               x- and y-directions respectively
    Return: (2,n) array of transformed coordinates
    '''
    A = np.diag(stretch)
    return np.dot(A, pts)

######################################################### Problem 2
def rotation(pts, angle):
    '''
    Do NOT plot
    Hint: Use helper function plot_old_new to view your transformations
    Inputs: 
    pts   -- (2,n) array where the first row represents the x-coordinates and the 
             second row the y-coordinates
    angle -- floating point number representing the angle of rotation in radians
    Return: (2,n) array of transformed coordinates
    '''
    a = np.cos(angle)
    b = np.sin(angle)
    R = np.array([[a, -b, b, a]]).reshape((2,2))
    return np.dot(R, pts)

######################################################### Problem 5
def translation(pts, shift):
    ''' 
    Do NOT plot
    Hint: Use helper function plot_old_new to view your transformations
    Inputs: 
    pts   -- (2,n) array where the first row represents the x-coordinates and the 
             second row the y-coordinates

    shift -- array with two entries indicating shift in x- and y-directions respectively
    Return: (2,n) array of transformed coordinates
    '''
    x_shift = [shift[0]]*len(pts[0])
    y_shift = [shift[1]]*len(pts[1])
    new = pts
    new[0] = pts[0] + x_shift
    new[1] = pts[1] + y_shift
    return new

######################################################### Problem 6
def trajectory(t, omega, v, s):
    ''' 
    Do NOT plot in this function!
    Use helper function plot_trajectory to view your particle
    Hint: you will use a composition of rotation and translation in this problem.
    Inputs: 
    t     -- time (seconds)
    omega -- angular velocity (radians/second)
    v     -- array indicating x- and y-directions of motion
    s     -- speed (meters/second)
    Return: array indicating the new position of p1 at time t
    '''
    p_2 = (s*t/np.sqrt(v[0]**2+v[1]**2))*v
    p_1 = np.array([1,0]).reshape((2,1))
    p_1 = rotation(p_1, t*omega)
    p_1 = translation(p_1, p_2)
    return p_1
    
def plot_trajectory(times = np.arange(0,10,.1), omega = np.pi, v = np.array([1,1]), s = 3):
    ''' 
    Plot the trajectory of p_1 as described in problem
    '''
    pos = np.zeros((2,len(times)))
    for i in xrange(len(times)):
        pos[:,i] = trajectory(times[i], omega, v, s).flatten()
    plt.plot(pos[0], pos[1])
    plt.show()

######################################################### Problem 7 
def REF(A):
    '''
    Hint: Use helper functions to perform elementary row operations 
    Inputs: 
    A -- Square matrix (numpy array) such that:
         - A is invertible
         - A zero will never appear on main diagonal during row reduction
         - Round-off errors may be ignored
           
    Change the matrix in place to compute the REF. 
    Do NOT return a value.
    '''
    n = A.shape[0]
    for c in xrange(n-1):
        for r in xrange(c+1, n):
            A[r,c:] = A[r,c:] - A[r,c]/float(A[c,c])*A[c,c:]
    return A
    
######################################################### Problem 8
def LU_decomp1(A):
    '''
    The matrix A should not be changed
    Inputs: 
    A -- Invertible matrix      
    Return: L,U as described in problem
    '''
    n = A.shape[0]
    U = np.copy(A)
    L = np.eye(n)
    for c in xrange(n-1):
        for r in xrange(c+1, n):
            L[r,c] = U[r,c]/float(U[c,c])
            U[r,c:] -= U[r,c]/float(U[c,c])*U[c,c:]
    return L, U

######################################################### Problem 10
def print_times():
    '''
    Fill in the times for questions 3 and 4. Discuss question 5 with your neighbor.
    '''

    #ANSWER 3
    lu_algorithm_time =  0.000969      #hard-code in the time you get (in seconds)
    
    #ANSWER 4
    inv_algorithm_time = 0.00732     #hard-code in the time you get (in seconds)
    
    print "LU solve:  %f" %lu_algorithm_time
    print "Invert:    %f" %inv_algorithm_time