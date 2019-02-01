"""
Programming Homework 2
Math 411
Winter Semester 2017

Note that a small portion of your grade will be based on the readability of 
your code. The spec files given to you are intended to be an example for you
to follow, but for additional information, you can refer to the 
PEP8 Python Style Guide here: https://www.python.org/dev/peps/pep-0008/

The line "pass" you see in the functions indicates that the function is empty,
and you should delete it once you start working.
"""

import numpy as np
from scipy import linalg as la


############### get_cheby_grid #########################################

def get_cheby_grid(a,b,N):
    """
    Returns the Chebyshev nodes of the second kind.
    
    INPUTS
    ------
    a - A float indicating the left end point of the real interval on which a 
        function is to be interpolated.
    b - A float indicating the right end point of the real interval on which a 
        function is to be interpolated.
    N - An int giving the bound on the degree of the Chebyshev polynomial interpolant, 
        indicating that there are N+1 spatial nodes for the interpolant.

    OUTPUTS
    -------
    x_nodes - An array of floats containing the Chebyshev points of the second
            kind: x_j = (a+b)/2+(a-b)*cos(j*pi/N)/2, where j = 0:1:N.
    """
    
    pass
   
############### get_cheby_coef ################################################

def get_cheby_coef(fx):
    """
    Find the coefficients of the Chebyshev interpolant 

    INPUTS
    ------
    fx - An array of floats containing f(x_j), j = 0,...,N where x_j are the
        Chebyshev nodes of the second kind and f is a function from [a,b] to 
        the complex numbers. 

    OUTPUTS
    -------
    coef    - An array of floats containg the coefficients c_j, j = 0,...,N to 
            the Chebyshev interpolant of f, p(x) = sum_{j=0}^N c_j T_j(x), 
            where T_j is the jth Chebyshev polynomial.
    """
    pass

############### evaluate_chebyshev_interpolant ################################

def evaluate_chebyshev_interpolant(coef,xpts,a,b):
    """
    Evaluates the chebyshev interpolant with coefficients coef at the points x

    INPUTS
    ------
    coef    - Array of doubles containing the coefficients of the Chebyshev 
            interpolant.
    xpts    - Array of doubles indicating the points at which the Chebyshev
            interpolant is to be evaluated.
    a       - A double indicating the left end point of the real interval on 
            which the function is interpolated by the Chebyshev interpolant.
    b       - A double indicating the right end point of the real interval on 
            which the function is interpolated by the Chebyshev interpolant.
    

    OUTPUTS
    -------
    px - A numpy array of floats which contains the values of the Chebyshev 
        interpolant evaluated at the points given in xpts.
    """
    pass
    
############### evaluate_chebyshev_interpolant ################################

def polynomial_interpolation(x,fx,xpts):
    """
    Evaluates the polyomial interpolant of a function at specified points using
    either the Newton or the Lagrange method.

    INPUTS
    ------
    x       - An array of floats indicating the interpolation nodes. Entries 
            should be distinct.
    fx      - An array of floats containing the values of a function f evaluated
            at the points given in x.
    xpts    - An array of floats containing the points at which the polynomial 
            interpolant is to be evaluated.
    
    OUTPUTS
    -------
    px - A numpy array of floats which contains the values of the polynomial
        interpolant evaluated at the points given in xpts
    """
    pass  
    
############### test ##########################################################

def test(N,num_points):
    """
    Evaluates the polyomial interpolant of a function at specified points using
    either the Newton or the Lagrange method.

    INPUTS
    ------
    N           - An int indicating the degree of polynomial interpolants. 
    num_points  - An int indicating the number of evenly spaced poitns on which
                the interpolants and functions will be graphed.
                
    INFORMATION
    -----------
    
    Runge's function is given by f(x) = 1/(1+25*X**2).
    
    OUTPUTS
    -------
    Plots all on the same graph, Runge's function (black) on the interval
    [-1,1], a Chebyshev interpolant (red) of Runge's function of degree N or 
    less, and the polynomial interpolant (blue) of Runge's function of degree N 
    or less that uses evenly spaced nodes. 
    """
    pass  
    
    
    
    
    
    
