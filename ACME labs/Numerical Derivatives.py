import numpy as np
from matplotlib import pyplot as plt

def cent_diff_quotients(f, pts, h=1e-5):
    """ 
    Calculate the centered difference quotient of a function at some given points.
    Inputs: 
        f   -- a callable function
        pts -- an array of points
        h   -- floating point
    Returns:
        array of centered difference quotients
    """
    return (1.0/2.0f(pts+h)-1.0/2.0f(pts-h))/h
    
def jacobian(f, m, n, pt, h = 1e-5):
    """
    Compute the approximate Jacobian matrix of a function at a point using the 
    centered coefficients difference quotient.
    Discuss part 2 with a partner.
    Inputs:
        f  -- a callable function
        m  -- dimension of the range of f
        n  -- dimension of the domain of f
        pt -- n-dimensional array
        h  -- floating point
    Returns:
        Jacobian matrix (numpy array)
    """
    pass

################################################################## Filters ####
gaussian_blur = 1./159 * np.array([[2.,  4.,  5.,  4.,  2.],
                                   [4.,  9.,  12., 9.,  4.],
                                   [5.,  12., 15., 12., 5.],
                                   [4.,  9.,  12., 9.,  4.],
                                   [2.,  4.,  5.,  4.,  2.]])

S = 1./8 * np.array([[-1.,0.,1.],[-2.,0.,2.],[-1.,0.,1.]])
###############################################################################


def Filter(image, filter):
    """
    Apply a filter to an image.
    Try question 2 and discuss with a partner.
    Inputs:
        image  -- an array of intensities representing an image
        filter -- an array representing the filter to apply
    Returns:
        array of the filtered image intensities
    """
    m, n = image.shape
    l, k = filter.shape
    
    # Create an array of zeros of the appropriate size
    
    
    # Make the interior of the image_pad equal to image
    
    
    C = np.empty_like(image)
    for i in xrange(n):
        for j in xrange(m):
            # Compute C[i,j]
            pass
    
    return C

def sobel(image):
    """
    Apply the Sobel filter to an image.
    Calculate the cutoff gradient (M) by using four times the average value of 
    the Euclidean norm of the partial derivatives.
    Inputs:
        image -- an array of intensities representing an image
    Returns:
        array after applying the sobel filter
    """
    pass


# In[ ]:



