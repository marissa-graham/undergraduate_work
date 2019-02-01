import numpy as np
import h5py
import copy
import sympy
import itertools
import operator
import string
import gauss_integrate
from fractions import Fraction

from scipy import linalg as la
from scipy.misc import comwb
from scipy.sparse import coo_matrix
from scipy.interpolate import griddata
import scipy.integrate
import scipy.special

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.transforms import offset_copy
from mpl_toolkits.mplot3d import axes3d

from NURBStools import computeSympyExtractionOperator

# no idea
def color_blind_list():
    return [
        (0, 0.45, 0.7),  # blue
        (0.0, 0.6, 0.5),  # bluish green
        (0.8, 0.4, 0.0),  # vermillion
        (0.35, 0.7, 0.9),  # sky blue
        (0.9, 0.6, 0.0),  # orange
        (0.8, 0.6, 0.7),  # reddish purple
        (0.95, 0.9, 0.25)  # yellow
        ]

# no idea
def ratio_markers(ratios, rot_angle = 0):
    markers = []
    accum = 0
    for i in range(len(ratios)):
        prev = accum
        accum += ratios[i]
        x = [0] + np.cos(np.linspace( 2 * np.pi * prev, 2 * np.pi * accum )).tolist()
        y = [0] + np.sin(np.linspace( 2 * np.pi * prev, 2 * np.pi * accum )).tolist()
        markers.append((list(zip(x,y)), rot_angle))
    return markers

#return a list with n copies of x
def dup(x, n):  
    return [ x for i in range(n) ]

# So it's taking a local knot vector, and just... I dunno
def local_knots_to_segment_rows( knots, fraction = False, rational = False, filter_zlen = True ):
    
    # calculuate the polynomial degree of the basis function
    n = len(knots) - 2
    # construct the extended knot vector
    k = dup(knots[0], n-2) + knots + dup(knots[-1], n-2)
    # initialize a list of zeros to hold the coefficients
    # the coefficients for the basis functions?
    p = dup(0, n * (n-1) + 1)
    for m in range(1, n):
        p_i = (m-1) * n
        p[p_i] = 0
        p[p_i + n - m] = 1
        # get the locations for knot insertion
        c1 = k[(m-1) + n -1]
        c2 = k[(m-1) + n ]
        # perform knot insertion n - 1 times at c1
        for i in range(n-1):     #i =[0...n-1)
            for j in range(n-i):   #j =[0...n-i)
                kij = k[(m-1) + i + j]
                knj = k[(m-1) + n + j]
                if kij != c1:
                    if fraction:
                        p[p_i + j] = frac((knj - c1) * p[p_i + j] + (c1 - kij) * p[p_i + j + 1],knj - kij)
                    elif rational:
                        p[p_i + j] = sympy.Rational((knj - c1) * p[p_i + j] + (c1 - kij) * p[p_i + j + 1],knj - kij)
                    else:
                        p[p_i + j] = ((knj - c1) * p[p_i + j] + (c1 - kij) * p[p_i + j + 1]) / float(knj - kij)
        # perform knot insertion n - 1 times at c2
        for i in range(n-1):          #i =[0...n-1)
            for j in range(n, i+1, -1): #j =[n..i+1) descending
                knij = k[(m-1) + n - i + j -1]
                if knij != c1:
                    if fraction:
                        p[p_i + j] = frac((knij - c2) * p[p_i + j -1] + (c2 - c1) * p[p_i + j], knij - c1)
                    elif rational:
                        p[p_i + j] = sympy.Rational((knij - c2) * p[p_i + j -1] + (c2 - c1) * p[p_i + j], knij - c1)
                    else:
                        p[p_i + j] = ((knij - c2) * p[p_i + j -1] + (c2 - c1) * p[p_i + j]) / float(knij - c1)
    # handle the linear case
    if n == 1: p = [1]

    # convert abbreviated extraction row to segment extraction rows
    extract_row = dup(0, n) + p + dup(0, n)
    segs = [ extract_row[i:i+n+1] for i in range(0, n * (n+1), n) ]

    # filter out zero length segments
    if filter_zlen:
        return [ seg for k1, k2, seg in zip(knots[0:], knots[1:], segs) if k1 != k2 ]
    else:
        return segs

# duh
def equal( a, b, tol = 1e-12 ):
    try:
        return all(abs(a-b) < tol)
    except TypeError:
        return abs(a-b) < tol

# B_i^n(x), easy enough
def bernstein_basis(n,i,x):  
    return comb(n,i) * ((x+1)/2.0) ** i * ((1-x)/2.0) ** (n-i)

# duh
def get_unique_knots( knot_vector ):
    unique_knots = [knot_vector[0]]
    for knot in knot_vector:
        if knot != unique_knots[-1]:
            unique_knots.append(knot)
    unique_knots = unique_knots
    return unique_knots

# duh
def get_knot_multiplicities( knot_vector ):
    unique_knots = get_unique_knots(knot_vector)
    multiplicities = [0 for i in range(len(unique_knots))]
    for knot in knot_vector:
        multiplicities[ unique_knots.index(knot) ] += 1
    # it would be more efficient for this one to return the unique knots as well
    return multiplicities

# Object is a set of rows, includes a method that gives the index where a certain row
# occurs in that set of rows 
class ExtractCache():
    
    def __init__(self, init_extract_rows = [], tol = 1e-12):
        self.extract_rows = init_extract_rows
        self.tol = tol
    def row_index( self, input_row ):
        
        for row, irow in zip(self.extract_rows, itertools.count()):
            # If all the elements in input_row and row are equal, return 
            if all(map(lambda x, y : equal(x,y,tol=self.tol), input_row, row)):
                return irow
        self.extract_rows.append(input_row)
        return len(self.extract_rows)-1

# Seems like a class for working with B-spline basis functions, but it doesn't seem
# to be being used elsewhere
class BFunc():
    def __init__(self, *knot_vectors):
        self.knot_vectors = knot_vectors
        self.dim = len(knot_vectors)
        self.extraction_rows = []
        self.degree = []
        self.unique_knots = []
        self.internal_funcs = []
        for idim in range(self.dim):
            self.degree.append(len(self.knot_vectors[idim]) - 2)
            self.unique_knots.append(get_unique_knots(self.knot_vectors[idim]))
            self.extraction_rows.append(local_knots_to_segment_rows( list(self.knot_vectors[idim]) ))
        
        def temp_func(x, idim):
            for i in range(len(self.unique_knots[idim])-1):
                x_min = self.unique_knots[idim][i]
                x_max = self.unique_knots[idim][i+1]
                if x >= x_min and x <= x_max:
                    dx = float(x_max - x_min)
                    result = 0
                    for j in range(0, self.degree[idim] + 1):
                        s = 2 * (x - x_min) / dx - 1
                        result += self.extraction_rows[idim][i][j] * bernstein_basis( self.degree[idim], j, s )
                    return result
                elif x < self.unique_knots[idim][0] or x > self.unique_knots[idim][-1]:
                    return 0.0
        self.internal_func = np.vectorize(temp_func)
    def __call__(self, *args):
        
        
        if len(args) != len(self.degree):
            raise ValueError("Number of input arguments must match the number of dimensions.")
        val_gen = ( self.internal_func(*vars) for vars in zip(args,itertools.count()))
        return reduce(operator.mul, val_gen, 1)
    def greville(self):
        pt = []
        for kv in self.knot_vectors:
            pt.append( np.mean( kv[1:-1] ) )
        return pt
    def anchor(self, type='center'):
        pt = []
        if type == 'center':
            for i, kv in enumerate(self.knot_vectors):
                center_start = (self.degree[i] + 1)/2
                if self.degree[i] & 1:
                    pt.append(kv[center_start])
                else:
                    pt.append( np.mean( kv[center_start:center_start+2] ) )
            return pt
        else:
            for kv in self.knot_vectors:
                pt.append( kv[0] )
            return pt
    def shifted_anchor( self, type = "center", split_space = 0.1 ):
        anchor_vec = self.anchor(type = type)
        ret_val = []
        for idim in range(self.dim):
            anchor = anchor_vec[ idim ]
            local_kv = self.knot_vectors[ idim ]
            knot_mults = get_knot_multiplicities( local_kv )
            unique_knots = get_unique_knots( local_kv )
            uind = None
            deg = len(local_kv) - 2 
            for i, uk in enumerate(unique_knots):
                if abs( uk - anchor ) < 1e-8:
                    uind = i
                    break
            if uind == None or knot_mults[ uind ] == 1:
                ret_val.append(anchor)
                continue
            elif uind == 0:
                mult = knot_mults[ uind ]
                mult_start = unique_knots[ uind ] + split_space / 2.0
                ret_val.append(mult_start + ( mult - deg + deg & 1 - (deg+1)&1) * split_space)
                continue
            elif uind == len(unique_knots)-1:
                mult = knot_mults[ uind ]
                mult_start = unique_knots[ uind ] - split_space / 2.0
                ret_val.append( mult_start - ( mult - deg + deg & 1 - (deg+1)&1) * split_space )
                continue
            elif uind == len(unique_knots) / 2 and not deg & 1:
                ret_val.append( anchor )
                continue
            #     if 
            else:
                mult = knot_mults[ uind ]
                offset = -1
                if knot_mults[uind+1] > knot_mults[uind-1]:
                    offset = 1
                mult_start = unique_knots[ uind ]
                ret_val.append( mult_start  + split_space * offset * mult / 4.0 )
                continue
        return ret_val

# Using sympy, make the matrix for transitioning from degree q to p, I think.
def degree_projection_op(p,q):
    x = sympy.Symbol('x')
    A = sympy.Matrix(q+1,q+1, lambda i,j:0)
    for i in range(q+1):
        for j in range(q+1):
            A[i,j] = sympy.integrate( bernstein_basis(q, i, x) * bernstein_basis(q, j, x), (x, -1, 1))
    M = sympy.Matrix(p+1,q+1, lambda i,j:0)
    for i in range(p+1):
        v = sympy.Matrix(q+1, 1, lambda i,j:0)
        for j in range(0, q+1):
            v[j,0] = sympy.integrate(bernstein_basis(p, i, x) * bernstein_basis(q, j, x), (x, -1, 1))
        for j in range(0, q+1):
            M[i,j] = (A.inv()*v)[j,0]
    return np.array(np.array(M), np.float)

# get local extraction operators. don't be fooled by the function name being singular
def computeExtractionOperator(U,p):
    
    """Compute the extraction operator relating a NURBS basis to the Bernstein polynomials over a single patch.

    Arguments:
    - `U`: Knot vector
    - `p`: Curve degree
    """
    m = len(U)
    a = p + 1
    b = a + 1
    nb = 0
    C = []
    C.append(np.eye(p+1))
    Unew = []
    for knot in U:
        Unew.append(knot)
    U = np.array(Unew)
    while b<m:
        C.append(np.eye(p+1))
        i = b
        mult = 1
        while b<m and equal(U[b],U[b-1]):
            mult = mult+1
            b = b+1
        # mult = b - i + 1
        if mult < p:
            numer = U[b-1]-U[a-1]
            alphas = np.zeros((p+1))
            for j in reversed(range(mult+1,p+2)):
                alphas[j-mult-1] = numer / float( U[a + j -1] - U[a-1] )
            r = p-mult
            for j in range(1,r+1):
                save = r - j
                s = mult + j
                for k in reversed(range(s,p+1)):
                    alpha = alphas[k-s]
                    C[nb][:,k] = alpha * C[nb][:,k] + (1 - alpha) * C[nb][:, k-1]
                if b<m:
                    C[nb+1][save:j+save+1,save] = C[nb][p-j:p+1,p]
            nb = nb + 1
        if b < m:
            a = b
            b = b + 1
    C.pop()
    return C

# Symbolic representation of a bernstein basis
class bernstein_basis_poly(sympy.Function):
    nargs = 3
    @classmethod
    def eval(cls, n, i, x):
        return sympy.binomial(n,i) * ((x+1)/2) ** i * ((1-x)/2) ** (n-i)

# Vector of bernstein basis functions (sympy), length p+1
class bernstein_basis_vector(sympy.Function):
    nargs = 2
    @classmethod
    def eval(cls, p, x):
        return sympy.Matrix(p+1,1,lambda i,j: bernstein_basis_poly(p,i, x))

# get the basis function defined by a local knot vector? but B-spline, not bezier, I think
class spline_basis():
    def __init__(self,local_knot_vector):
        
        
        self.extraction_rows = local_knots_to_segment_rows( list(local_knot_vector) )
        self.degree = len(local_knot_vector) - 2
        self.knot_vector = local_knot_vector
        self.unique_knots = [local_knot_vector[0]]
        for knot in local_knot_vector:
            if knot != self.unique_knots[-1]:
                self.unique_knots.append(knot)
        self.unique_knots = np.array(self.unique_knots)
        def temp_func(x):
            for i in range(len(self.unique_knots)-1):
                x_min = self.unique_knots[i]
                x_max = self.unique_knots[i+1]
                if x >= x_min and x <= x_max:
                    dx = float(x_max - x_min)
                    result = 0
                    for j in range(0, self.degree + 1):
                        s = 2 * (x - x_min) / dx - 1
                        result += self.extraction_rows[i][j] * self.bernstein( s, j )
                    return result
                elif x < self.unique_knots[0] or x > self.unique_knots[-1]:
                    return 0.0
        self.internal_func = np.vectorize(copy.copy(temp_func))
    def bernstein(self, t, i):
        
        if i>=0:
            return 1.0 / 2**self.degree * comb(self.degree, i)*(1-t)**(self.degree-(i)) * (1+t)**(i)
        else:
            return 0.0

    def __call__(self, x):
        return self.internal_func(x)

# get all the local knot vectors? like, U[:4],U[1:5], etc., I think
def window(iterable, size):
    iters = itertools.tee(iterable, size)
    for i in xrange(1, size):
        for each in iters[i:]:
            next(each, None)
    return itertools.izip(*iters)

# Get the basis functions from all the local knot vectors
# Oh, we're doing it the local knot vector way instead of the 
def spline_basis_funcs( degree, knot_vector ):
    funcs = []
    for local_kv in window(knot_vector, degree + 2):
        funcs.append(spline_basis(local_kv))
    return funcs

# This is where we get the indices of the basis functions that contribute to each
# segment! Like in getting the extraction operators! I GOT IT
def build_1d_connects(knot_vector, degree):
    connects = []
    mults = get_knot_multiplicities(knot_vector)
    ukv = get_unique_knots(knot_vector)
    count = 0
    # for each unique knot
    for i in xrange(len(ukv)-1):

        connects.append([j for j in range(count, count + degree + 1)])
        count += mults[i+1]
    return connects

# Tells you which of the segments each basis function contributes to, I think
# connects should be the output from build_1d_connects
def build_1d_rev_connects(connects):
    # reverse connect?
    rev_connects = [[] for i in range(connects[-1][-1]+1)]

    # for i in xrange(len(connects)): c = connects[i] equivalent
    for i, c in enumerate(connects):
        # so each item in connects is a list
        for j in c:
            rev_connects[j].append(i)
    return rev_connects

# Get the dual extraction operator
# I didn't see this in the paper, how relevant is it? Nothing else calls it
def spline_dual_ext_ops( knot_vector, degree ):
    ukv = get_unique_knots(knot_vector)
    ext_ops = computeExtractionOperator(knot_vector, degree)
    G = np.matrix(sympy_gramian(degree))
    Ginv = np.linalg.inv(G)
    dual_ext_ops = []
    for C in ext_ops:
        Cinv = np.linalg.inv(np.matrix(C))
        dual_ext_ops.append(Cinv.T * Ginv.T)
    connects = build_1d_connects(knot_vector, degree)
    rev_connects = build_1d_rev_connects(connects)
    funcs = []
    weights = [[0 for j in range(degree+1)] for i in range(len(connects))]
    sum_vals = [0 for i in range(len(rev_connects))]
    for i, (c, C) in enumerate(zip(connects, ext_ops)):
        elem_size = ukv[i+1] - ukv[i]
        for j, jfunc in enumerate(c):
            row_sum = sum(C[j,:])
            weights[i][j] += row_sum * elem_size
            sum_vals[jfunc] += row_sum * elem_size
    for i, c in enumerate(connects):
        for j, jfunc in enumerate(c):
            weights[i][j] /= sum_vals[jfunc]
    for i, rc in enumerate(rev_connects):
        for e in rc:
            ifunc_l = connects[e].index(i)
            dual_ext_ops[e][ifunc_l,:] = dual_ext_ops[e][ifunc_l,:] * weights[e][ifunc_l] * 1 / float(ukv[e+1]-ukv[e])
    return dual_ext_ops
    
# Symbolic version of spline_dual_ext_ops, also not 
def sympy_spline_dual_ext_ops( knot_vector, degree ):
    
    
    ukv = get_unique_knots(knot_vector)
    ext_ops = computeSympyExtractionOperator(knot_vector, degree)
    G = sympy_gramian(degree)
    Ginv = G.inv()
    dual_ext_ops = []
    for C in ext_ops:
        Cinv = C.inv()
        dual_ext_ops.append(Cinv.T * Ginv.T)
    connects = build_1d_connects(knot_vector, degree)
    rev_connects = build_1d_rev_connects(connects)
    funcs = []
    weights = [[0 for j in range(degree+1)] for i in range(len(connects))]
    sum_vals = [0 for i in range(len(rev_connects))]
    for i, (c, C) in enumerate(zip(connects, ext_ops)):
        elem_size = ukv[i+1] - ukv[i]
        for j, jfunc in enumerate(c):
            row_sum = sum(C[j,:])
            weights[i][j] += row_sum * elem_size
            sum_vals[jfunc] += row_sum * elem_size
    for i, c in enumerate(connects):
        for j, jfunc in enumerate(c):
            weights[i][j] /= sum_vals[jfunc]
    for i, rc in enumerate(rev_connects):
        for e in rc:
            ifunc_l = connects[e].index(i)
            dual_ext_ops[e][ifunc_l,:] = dual_ext_ops[e][ifunc_l,:] * weights[e][ifunc_l] * sympy.Rational(1, (ukv[e+1]-ukv[e]))
    return dual_ext_ops

# getting the dual basis functions, I don't know why that's useful
# something to do with the inverse gramian, maybe?
def spline_dual_funcs( degree, knot_vector, default_val = 0.0 ):
    
    ukv = get_unique_knots(knot_vector)
    ext_ops = computeExtractionOperator(knot_vector, degree)
    G = gramian(degree)
    Ginv = np.linalg.inv(G)
    dual_ext_ops = []
    for C in ext_ops:
        Cinv = np.linalg.inv(C)
        dual_ext_ops.append(np.dot(Cinv.T, Ginv.T))
    connects = build_1d_connects(knot_vector, degree)
    rev_connects = build_1d_rev_connects(connects)
    funcs = []
    weights = [[0 for j in range(degree+1)] for i in range(len(connects))]
    sum_vals = [0 for i in range(len(rev_connects))]
    for i, (c, C) in enumerate(zip(connects, ext_ops)):
        elem_size = ukv[i+1] - ukv[i]
        for j, jfunc in enumerate(c):
            row_sum = sum(C[j,:])
            weights[i][j] += row_sum * elem_size
            sum_vals[jfunc] += row_sum * elem_size
    for i, c in enumerate(connects):
        for j, jfunc in enumerate(c):
            weights[i][j] /= sum_vals[jfunc]
    for i, rc in enumerate(rev_connects):
        func_ext_rows = []
        elem_lines = []
        for e in rc:
            ifunc_l = connects[e].index(i)
            func_ext_rows.append(dual_ext_ops[e][ifunc_l,:] * weights[e][ifunc_l] * 2.0 / (ukv[e+1]-ukv[e]))
            elem_lines.append(ukv[e])
        elem_lines.append(ukv[rc[-1]+1])
        funcs.append(ExtractionFunc(degree, func_ext_rows, elem_lines, default_val))
    return funcs

# Object is a bernstein basis of a certain degree, you can call it on B_i^degree(t) 
# and there's an integrate function that may or may not be useful
class BernsteinBasis():
    def __init__(self, degree):
        self.degree = degree
    def __call__(self, t, i):
        
        if i>0:
            return 1.0 / 2**self.degree * comb(self.degree, i-1)*(1-t)**(self.degree-(i-1)) * (1+t)**(i-1)
        else:
            return 0.0

    def integrate(self, i):
        return 2 / float(self.degree + 1)

# Evaluate a Bezier curve at a point x with control points coeffs
def bern_eval(coeffs, x):
    deg = len(coeffs) - 1
    B = BernsteinBasis(deg)
    ret = x*0
    for i in range(deg+1):
        ret += B(x,i+1) * coeffs[i]
    return ret

# Evaluate all the bernstein basis functions of a certain degree at all of the pts
def eval_bernstein_pts( degree, pts ):
    
    vals = np.zeros((len(pts), degree + 1))
    B = BernsteinBasis(degree)
    for j in range(degree + 1):
        for i in range(len(pts)):
            vals[i,j] = B(pts[i],j+1)
    return vals

# Bezier elements, used in l2_project etc.
# not sure what elem_connect and elem_ext are
class Elem():
    def __init__(self, x_min, x_max, elem_connect, elem_ext):
        self.x_min = float(x_min)
        self.x_max = float(x_max)
        self.elem_connect = elem_connect
        self.elem_ext = elem_ext
        self.degree = elem_ext.shape[1] - 1
        self.globalFuncN = elem_ext.shape[0]
        self.dx = self.x_max - self.x_min

    def eval_basis(self, t):
        B = BernsteinBasis(self.degree)
        Bval = np.zeros(self.degree+1)
        for i in range(self.degree + 1):
            Bval[i] = B(t,i+1)
        return np.dot(self.elem_ext,Bval)

    def eval_func(self, t, ifunc):
        if self.elem_connect.count(ifunc) == 0:
            return t * 0
        else:
            ifunc_l = self.elem_connect.index( ifunc )
            B = BernsteinBasis(self.degree)
            Bval = np.zeros(self.degree+1)
            for i in range(self.degree + 1):
                Bval[i] = B(t,i+1)
            return np.dot(self.elem_ext[ifunc_l],Bval)

    def integrate_func(self,ifunc_l):
        return sum(self.elem_ext[ifunc_l]) * self.dx * 2 / float(self.degree+1)

    def integrate_func_sq(self, ifunc_l):
        
        ext_row = self.elem_ext[ifunc_l]
        bern_prod_coeffs = []
        for i in range(self.degree + 1):
            for j in range(i, self.degree + 1):
                bern_prod_coeffs.append( comb(self.degree, i) * comb(self.degree, j) / float(comb(2 * self.degree, i + j) ) * ext_row[i] * ext_row[j])
        return sum(bern_prod_coeffs) * self.dx * 2 / float( 2 * self.degree+1)

    def global_coord( self, lpt ):
        return (self.x_max - self.x_min) * ( lpt + 1 ) / 2.0 + self.x_min

    def local_coord( self, global_pt):
        return 2 * (global_pt - self.x_min) / ( self.x_max - self.x_min ) - 1

    def localize( self, global_coeffs ):
        local_coeffs = np.zeros(len(self.elem_connect))
        for i in range(len(self.elem_connect)):
            local_coeffs[i] = global_coeffs[self.elem_connect[i]]
        return local_coeffs

    def eval_sol( self, local_coeffs, num_pts = 3, global_pts = None, local_pts = None ):
        if global_pts != None:
            sol = np.zeros(global_pts.shape)
            for i in range(len(global_pts)):
                sol[i] = np.dot(self.eval_basis(global_pts[i]), local_coeffs)
            return global_pts, sol

        if local_pts != None:
            sol = np.zeros(len(local_pts))
            for i in range(len(local_pts)):
                sol[i] = np.dot(self.eval_basis(local_pts[i]), local_coeffs)
            return self.global_coord(np.array(local_pts)), sol

        else:
            t = np.linspace(-1,1,num_pts)
            sol = np.zeros(t.shape)
            for i in range(len(t)):
                sol[i] = np.dot(self.eval_basis(t[i]), local_coeffs)
            return self.global_coord(t), sol

# OHHHHHH it's a tensor product, I see. Tensor product all the vecs
def ten_prod(*vecs):
    b_vals = []
    
    letters = (l for l in string.lowercase)
    einsum_string = ""
    for i in range(len(vecs)):
        if len(einsum_string) == 0:
            einsum_string =  "" + letters.next()
        else:
            einsum_string = einsum_string + "...," + letters.next()
    einsum_string += "..."
    return np.einsum(einsum_string, *vecs)

# Spline function object. Has basis functions, knot vectors, control pts, weights, etc.
# and then you can call it to get function values
class SplineFunc():
    
    def __init__(self, degrees, kvs, pts = None, homogeneous_pts = None, pts_wts = None):
        self.degrees = degrees
        self.kvs = kvs
        self.basis = []
        self.num_bfuncs = 1
        self.num_bfuncs_dim = []
        for deg, kv in zip(degrees,kvs):
            self.basis.append(spline_basis_funcs(deg, kv))
            self.num_bfuncs_dim.append(len(kv) - (deg + 1))
            self.num_bfuncs *= self.num_bfuncs_dim[-1]
        if pts != None:
            self.pts = pts
            self.wts = None
        if homogeneous_pts != None:
            self.pts = [np.array(pt[0:-1]) / pt[-1] for pt in homogeneous_pts]
            self.wts = [pt[-1] for pt in homogeneous_pts]
        if pts_wts != None:
            self.pts = [np.array(pt[0:-1]) for pt in pts_wts]
            self.wts = [pt[-1] for pt in pts_wts]
        
    def set_pts(self, pts):
        self.pts = pts
    def set_wts(self, wts):
        self.wts = wts
    def __call__(self, *param_pos):
        
        ret_val = np.zeros(param_pos[0].shape + np.array(self.pts[0]).shape)
        wt_val = 0
        for ifunc in range(self.num_bfuncs):
            ifunc_dims = np.unravel_index(ifunc, self.num_bfuncs_dim, order = 'F')
            func_val = self.basis[0][ifunc_dims[0]](param_pos[0])
            for idim in range(1,len(self.degrees)):
                func_val *= self.basis[idim][ifunc_dims[idim]](param_pos[idim])
            if self.wts == None:
                ret_val += func_val[..., np.newaxis] * np.array(self.pts[ifunc])
            else:
                ret_val += self.wts[ifunc ] * func_val[..., np.newaxis] *  np.array(self.pts[ifunc])
                wt_val += self.wts[ifunc] * func_val

        if self.wts != None:
            ret_val = np.einsum('a...,a...->a...', ret_val, 1 / wt_val)
        if ret_val.shape[0] == 1:
            return ret_val[0]
        else:
            return ret_val

# B-spline element objects with lots of features
class BElem():
    def __init__(self, spatial_dim, param_dim, elem_info, func_cache = None):
        # elem_info has a lot of stuff in it that's not being mentioned here
        # but it's a dictionary that includes degrees, ext_op, global_func_n, lnodes
        # elem_connect, bounds

        # ext_info includes
        self.elem_info = copy.deepcopy(elem_info)
        self.local_func_n = 1
        for deg in self.elem_info['degrees']:
            self.local_func_n *= (deg + 1)
        self.param_dim = param_dim
        self.spatial_dim = spatial_dim
        if func_cache == None:
            self.func_cache = {}
        else:
            self.func_cache = func_cache
        gauss_1d = []
        wts_1d = []
        for deg in self.elem_info['degrees']:
            gpts, wts = gauss_integrate.get1dGaussRule(deg + 1)
            gauss_1d.append(gpts)
            wts_1d.append(wts)
        self.wts_1d = wts_1d
        self.gauss_1d = gauss_1d
        self._bounds = []
        for i in range(self.param_dim):
            param_pos = [[-1] for j in range(self.param_dim)]
            ll = self.eval(*param_pos)
            s_min = np.squeeze(ll)[i]
            param_pos = [[1] for j in range(self.param_dim)]
            ur = self.eval(*param_pos)
            s_max = np.squeeze(ur)[i]
            self._bounds.append([s_min, s_max])

    def global_coord(self, lpt):
        ret=[]
        for idim in range(self.param_dim):
            ret.append((self.bounds(idim)[1]-self.bounds(idim)[0])*(lpt[idim]+1)/2.0+self.bounds(idim)[0])
        return ret

    def local_coord(self, global_pt):
        ret=[]
        for idim in range(self.param_dim):
            ret.append(2.0*(global_pt[idim]-self.bounds(idim)[0])/float(self.bounds(idim)[1]-self.bounds(idim)[0])-1)
        return ret

    def ext_op(self):
        return self.elem_info['ext_op']

    def bernstein_no_vec(self,degree, i, s):
        # compute the value (at a single point)
        u1 = 0.5 * ( 1 - s )
        u2 = 0.5 * ( 1 + s )
        temp = np.zeros(degree + 1)
        temp[ degree - i ] = 1.0
        for ii in range(1, degree+1):
            for jj in range(degree, ii-1, -1):
                temp[ jj ] = u1 * temp[ jj ] + u2 * temp[ jj - 1 ]

        val = temp[ degree ]
        return val

    def bernstein(self,degree,i,t):
        return np.vectorize(self.bernstein_no_vec, excluded = ('self', 'degree', 'i'))(degree, i, t)
    
    # I think this is equipped for 2d, so it does the first dimension, and then calls
    # that on the second dimension
    def bernstein_d1(self,degree,i,s):
        val = np.zeros(np.shape(s))
        if i > 0:   val += self.bernstein( degree - 1,i - 1,  s )
        if i < degree: val -= self.bernstein( degree - 1, i, s )
        return val * degree / 2.0;

    def bernstein_d2(self,degree,i,s):
        val = np.zeros(np.shape(s))
        if i > 0:   val += self.bernstein_d1( degree - 1, i - 1, s )
        if i < degree: val -= self.bernstein_d1( degree - 1, i, s )
        return val * degree / 2.0;

    def eval_local_basis_no_weight(self, ifunc, *param_pos, **kwargs):
        if kwargs.has_key('derivs'):
            derivs = kwargs.pop('derivs')
        else:
            derivs = [0]
        if len(derivs) == 0:
            derivs = [0 for i in range(self.param_dim)]
        elif len(derivs) < self.param_dim:
            for i in range(self.param_dim - len(derivs)):
                derivs.append(0)
        cache = False
        if len(param_pos) == 0:
            if self.func_cache.has_key((ifunc,tuple(derivs))):
                return self.func_cache[(ifunc,tuple(derivs))]
            else:
                param_pos = self.gauss_1d
                cache = True
        ifunc_dims = np.unravel_index(ifunc, [p + 1 for p in self.elem_info['degrees']], order = 'F')
        b_vals = []
        
        letters = (l for l in string.lowercase)
        # einsum_string = ""
        for i in range(self.param_dim):
            # if len(einsum_string) == 0:
            #     einsum_string =  "" + letters.next()
            # else:
            #     einsum_string = einsum_string + "," + letters.next()
            xi = param_pos[i]
            deg = self.elem_info['degrees'][i]
            if derivs[i] == 0:
                b_vals.append(self.bernstein(deg, ifunc_dims[i], xi))
            elif derivs[i] == 1:
                b_vals.append(self.bernstein_d1(deg, ifunc_dims[i], xi))
            elif derivs[i] == 2:
                b_vals.append(self.bernstein_d2(deg, ifunc_dims[i], xi))
        b_vals.reverse()
        ret_val = ten_prod(*b_vals)
        if cache:
            self.func_cache[ifunc,tuple(derivs)] = ret_val
        return ret_val
    def eval_weight(self, *param_pos, **kwargs):
        if kwargs.has_key('derivs'):
            derivs = kwargs.pop('derivs')
        else:
            derivs = [0]

        # ret_val = np.zeros(np.shape(param_pos[0]))
        ret_val = None
        for lnode, i in zip(self.elem_info['lnodes'], itertools.count()):
            if ret_val == None:
                ret_val = lnode[-1] * self.eval_local_basis_no_weight( i, derivs = derivs , *param_pos)
            else:
                ret_val += lnode[-1] * self.eval_local_basis_no_weight( i, derivs = derivs , *param_pos )
        return ret_val
    def eval_local_basis(self, ifunc, *param_pos, **kwargs):
        if kwargs.has_key('derivs'):
            derivs = kwargs.pop('derivs')
        else:
            derivs = [0]
        ifunc_dims = np.unravel_index(ifunc, [p + 1 for p in self.elem_info['degrees']], order = 'F')
        N = self.eval_local_basis_no_weight( ifunc, *param_pos )
        W = self.eval_weight(*param_pos)
        if sum(derivs) == 0:
            return N / W
        elif sum(derivs) == 1:
            dN = self.eval_local_basis_no_weight( ifunc, derivs = derivs, *param_pos )
            dW = self.eval_weight( *param_pos, derivs = derivs )
            return dN / W - N * dW / W ** 2
    def eval_global_basis(self, ifunc, *param_pos, **kwargs):
        if kwargs.has_key('derivs'):
            derivs = kwargs.pop('derivs')
        else:
            derivs = [0]
        
        ret_val = None
        for i in range(self.local_func_n):
            if ret_val == None:
                ret_val = self.elem_info['ext_op'][ifunc, i] * self.eval_local_basis( i, derivs = derivs, *param_pos )
            else:
                ret_val += self.elem_info['ext_op'][ifunc, i] * self.eval_local_basis( i, derivs = derivs, *param_pos )
        return ret_val
    def eval(self, *param_pos, **kwargs ):
        if kwargs.has_key('derivs'):
            derivs = kwargs.pop('derivs')
        else:
            derivs = [0]
        if len(param_pos) == 0:
            ret_val = [np.zeros([len(pp) for pp in self.gauss_1d]) for i in range(len(self.elem_info['lnodes'][0])-1)]
        else:
            ret_val = [np.zeros([len(pp) for pp in param_pos]) for i in range(len(self.elem_info['lnodes'][0])-1)]
        for ifunc in range(self.local_func_n):
            func_val = self.eval_local_basis(ifunc, derivs = derivs, *param_pos)
            lnode = self.elem_info["lnodes"][ifunc]
            for i in range(len(lnode)-1):
                ret_val[i] += lnode[i] * func_val
        return ret_val
    def eval_jac(self, *param_pos):
        # dx/ds
        spatial_dim = self.spatial_dim
        param_dim = self.param_dim
        deriv_vals = []
        for j in range(param_dim):
            derivs = [0 for i in range(param_dim)]
            derivs[j] = 1
            deriv_vals.append(self.eval(derivs = derivs, *param_pos))
        jac = []
        for i in range(spatial_dim):
            new_row = []
            for j in range(param_dim):
                new_row.append(deriv_vals[j][i])
            jac.append(new_row)
        return np.array(jac)
    def eval_jac_det(self, *param_pos):
        jac = self.eval_jac(*param_pos)
        return np.linalg.det(jac.T).T
    def localize(self, coeffs):
        return coeffs[tuple(self.elem_info['elem_connect'])]
    def eval_sol(self, coeffs, *param_pos):
        local_coeffs = self.localize(coeffs)
        ret_val = None
        gfunc_n = self.elem_info['global_func_n']
        for i in range(gfunc_n):
            if ret_val == None:
                ret_val = self.eval_global_basis(i, *param_pos)
            else:
                ret_val += self.eval_global_basis(i, *param_pos)
        return ret_val
    def bounds(self, idim = None):
        if idim == None:
            return self._bounds
        else:
            return self._bounds[idim]

    def param_vol(self):
        pvol = 1
        for idim in range(self.param_dim):
            p_bounds = self.bounds(idim)
            pvol *= (p_bounds[1] - p_bounds[0])
        return pvol

# This is like (60) from page 72, in section 3.1.1, I think
# I dunno what the jacobian has to do with anything, it doesn't show up in the paper
def LocalAssemblyL2New( elem, func ):
    # Assuming we're taking in an element built by BElem()? Or Elem? Or either?
    gfunc_n = elem.elem_info['global_func_n']
    k = np.zeros((gfunc_n,gfunc_n))
    f = np.zeros(gfunc_n)
    gauss_1d = []
    wts_1d = []
    for deg in elem.elem_info['degrees']:
        gpts, wts = gauss_integrate.get1dGaussRule(deg + 1)
        gauss_1d.append(gpts)
        wts_1d.append(wts)
    wt = ten_prod(*wts_1d)
    spatial_coords = elem.eval(*gauss_1d)
    f_val = func(*spatial_coords)
    jac = elem.eval_jac_det()
    for a in range(gfunc_n):
        # Not quite sure why the jacobian is involed
        N_a = elem.eval_global_basis(a)
        f_prod = N_a * f_val * jac * wt
        f[a] = np.sum(f_prod)
        for b in range(gfunc_n):
            N_b = elem.eval_global_basis(b)
            k_prod = N_a * N_b * jac * wt
            k[a,b] = np.sum(k_prod)
    return k, f

# Not entirely sure what this does or why it's helpful
def SparseAssembleL2( ext_info, local_assembly, func ):
    k_row = []
    k_col = []
    k_data = []
    f_row = []
    f_col = []
    f_data = []
    node_n = ext_info['nodeN']
    func_cache = {}
    K = np.zeros((node_n,node_n))
    F = np.zeros(node_n)
    for elem_info in ext_info['elems']:
        elem = BElem(ext_info['spatial_dim'],ext_info['param_dim'],elem_info, func_cache = func_cache)
        # I'm assuming this would be, like, LocalAssemblyL2New or something
        k, f = local_assembly( elem, func )
        for a in range(elem.elem_info['global_func_n']):
            A = elem.elem_info['elem_connect'][a]
            f_row.append(A)
            f_col.append(0)
            f_data.append(f[a])
            F[A] += f[a]
            for b in range(elem.elem_info['global_func_n']):
                
                B = elem.elem_info['elem_connect'][b]
                k_row.append(A)
                k_col.append(B)
                k_data.append(k[a,b])
                K[A,B] += k[a,b]
    
    # coo_matrix: a sparse matrix in coordinate format
    return coo_matrix((k_data, (k_row, k_col))).tocsr(), coo_matrix((f_data, (f_row, f_col))).todense()
    return K,F

# Seems like it's trying to put all of the elements into a global function
def func_from_ext( ext_info, ifunc, num_pts = 20, **kwargs ):
    ax = plt.gca()
    param_vals_1d = []
    # try:
    #     for np in num_pts:
    #         param_vals_1d.append(np.linspace(-1,1,np))
    # except TypeError:
    for i in ext_info['elems'][0]['degrees']:
        param_vals_1d.append(np.linspace(-1,1,num_pts))
    param_vals = np.meshgrid(*param_vals_1d)
    ielems = ext_info['func_map'][ifunc]['elems']
    ifuncs_l = ext_info['func_map'][ifunc]['local_ifuncs']
    # x = None
    x = []
    y = []
    vals = []
    for ielem, ifunc_l in zip(ielems, ifuncs_l):
        elem = BElem(ext_info['spatial_dim'], ext_info['param_dim'], ext_info['elems'][ielem])
        b_val = elem.eval_global_basis(ifunc_l, *param_vals_1d)
        spatial_coords = elem.eval(*param_vals_1d)
        x.extend([xi for xi in spatial_coords[0].flat])
        y.extend([xi for xi in spatial_coords[1].flat])
        vals.extend([v for v in b_val.flat])
    x = np.array(x)
    y = np.array(y)
    vals = np.array(vals)
    
    xn = np.linspace(x.min(), x.max(), num_pts)
    yn = np.linspace(y.min(), y.max(), num_pts)
    xa,ya = np.meshgrid(xn,yn)
    spl_func = lambda xa, ya : griddata( (x, y), vals, (xa, ya), fill_value = 0)
    return spl_func

# The main difference between this and func_from_ext is it uses coeffs instead of
# ifunc, so whatever that's supposed to mean
def func_from_sol( ext_info, coeffs, num_pts = 10, **kwargs ):
    
    ax = plt.gca()
    param_vals_1d = []
    # try:
    #     for np in num_pts:
    #         param_vals_1d.append(np.linspace(-1,1,np))
    # except TypeError:
    for i in ext_info['elems'][0]['degrees']:
        param_vals_1d.append(np.linspace(-1,1,num_pts))
    param_vals = np.meshgrid(*param_vals_1d)
    # x = None
    x = []
    y = []
    vals = []
    ax = plt.gca()
    for elem_info in ext_info['elems']:
        elem = BElem(ext_info['spatial_dim'], ext_info['param_dim'], elem_info)
        sol_val = None
        # elem_fig = plt.figure()
        for a in range(elem.elem_info['global_func_n']):
            # ax_func = elem_fig.add_subplot(ext_info['degrees'][0] + 1, ext_info['degrees'][0] + 1, a + 1)
            A = elem.elem_info['elem_connect'][a]
            N_a = elem.eval_global_basis(a, *param_vals_1d)
            if sol_val == None:
                sol_val = coeffs[A] * N_a
            else:
                sol_val += coeffs[A] * N_a
        spatial_coords = elem.eval(*param_vals_1d)

        # you're just extending the list it's not some special function
        x.extend([xi for xi in spatial_coords[0].flat])
        y.extend([xi for xi in spatial_coords[1].flat])
        vals.extend([sv for sv in sol_val.flat])
    x = np.array(x)
    y = np.array(y)
    vals = np.array(vals)
    
    xn = np.linspace(x.min(), x.max(), num_pts)
    yn = np.linspace(y.min(), y.max(), num_pts)
    xa,ya = np.meshgrid(xn,yn)
    spl_func = lambda xa, ya : griddata( (x, y), vals, (xa, ya), fill_value = 0, method = "cubic")
    return spl_func

# Same as func_from_ext, just plot it as well as returning it. Pretty sure
def plot_func_from_ext( ext_info, ifunc, num_pts = 20, **kwargs ):
    ax = plt.gca()
    param_vals_1d = []
    # try:
    #     for np in num_pts:
    #         param_vals_1d.append(np.linspace(-1,1,np))
    # except TypeError:
    for i in ext_info['elems'][0]['degrees']:
        param_vals_1d.append(np.linspace(-1,1,2))
    param_vals = np.meshgrid(*param_vals_1d)
    ielems = ext_info['func_map'][ifunc]['elems']
    ifuncs_l = ext_info['func_map'][ifunc]['local_ifuncs']
    # x = None
    x = []
    y = []
    vals = []
    for ielem, ifunc_l in zip(ielems, ifuncs_l):
        elem = BElem(ext_info['spatial_dim'], ext_info['param_dim'], ext_info['elems'][ielem])
        b_val = elem.eval_global_basis(ifunc_l, *param_vals_1d)
        spatial_coords = elem.eval(*param_vals_1d)
        x.extend([xi for xi in spatial_coords[0].flat])
        y.extend([xi for xi in spatial_coords[1].flat])
        vals.extend([v for v in b_val.flat])
    x = np.array(x)
    y = np.array(y)
    vals = np.array(vals)
    xn = np.linspace(x.min(), x.max(), num_pts)
    yn = np.linspace(y.min(), y.max(), num_pts)
    xa,ya = np.meshgrid(xn,yn)
    spl_func = func_from_ext( ext_info, ifunc, num_pts = num_pts)
    za = spl_func(xa,ya)

    if ax.properties().has_key("ylim3d"):
        ax.plot_surface(xa,ya, za, **kwargs)
    else:
        ax.pcolormesh(xa,ya,za, vmin = 0, vmax = 1, **kwargs)
    return spl_func

# Looks very similar LocalAssemblyL2 but with slightly different datatypes?
# Both of them are used in other functions, but are they interchangeable?
def LocalAssemblyL2( elem, func, gauss_pts, wts ):
    gfunc_n = elem.globalFuncN
    jac = float(elem.x_max - elem.x_min)

    k = np.zeros((gfunc_n,gfunc_n))
    f = np.zeros(gfunc_n)

    func_val = np.zeros(len(gauss_pts))
    N_vals = np.zeros((gfunc_n, len(gauss_pts)))
    for i in range(len(gauss_pts)):
        gpt = gauss_pts[i]
        global_pt = elem.global_coord(gpt)
        func_val[i] = func(global_pt)
        N_vals[:,i] = elem.eval_basis(gpt)
    for b in range(gfunc_n):
        for i in range(len(gauss_pts)):
            f[b] += func_val[i] * N_vals[b,i] * jac * wts[i]
        for a in range(gfunc_n):
            for i in range(len(gauss_pts)):
                k[a,b] += N_vals[a,i] * N_vals[b,i] * jac * wts[i]
    return k, f

# Think we're just getting error approximations here
def LocalSqDiffL2( elem, local_coeffs, func, gauss_pts, wts ):
    gfunc_n = elem.globalFuncN
    jac =  float( elem.x_max - elem.x_min )

    sq_diff = 0

    for i in range(len(gauss_pts)):
        gpt = gauss_pts[i]
        global_pt = elem.global_coord(gpt)
        func_val = func(global_pt)
        pts, sol = elem.eval_sol( local_coeffs, local_pts = [gpt] )
        sol_val = sol[0]
        sq_diff += (func_val - sol_val) ** 2 * jac * wts[i]
    return sq_diff

# I guess it's just putting all the elements together back into a global thing
def AssembleL2( ndofs, elems, func, gauss_rule ):
    K = np.zeros((ndofs,ndofs))
    F = np.zeros(ndofs)
    gpts, wts = gauss_integrate.get1dGaussRule(gauss_rule)
    for elem in elems:
        ke, fe = LocalAssemblyL2( elem, func, gpts, wts )
        for b in range(elem.globalFuncN):
            B = elem.elem_connect[b]
            F[ B ] += fe[b]
            for a in range(elem.globalFuncN):
                A = elem.elem_connect[a]
                K[A,B] += ke[a,b]
    return K, F

# pretty self-explanatory
def L2Error( elems, coeffs, func, gauss_rule ):
    
    sq_diff = 0
    gpts, wts = gauss_integrate.get1dGaussRule(gauss_rule)
    x_min = elems[0].x_min
    x_max = elems[0].x_max
    for elem in elems:
        local_coeffs = elem.localize(coeffs)
        sq_diff += LocalSqDiffL2( elem, local_coeffs, func, gpts, wts )
        x_min = min( x_min, elem.x_min )
        x_max = max( x_max, elem.x_max )
    return np.sqrt(sq_diff / (x_max - x_min))

# Projection to get the coefficients for an element (probably not a spline)
# (but maybe a spline, if you're going to a smaller space)
def l2_project( func, domain, degree, num_elems = 10, knot_vector = None ):
    if knot_vector == None:
        knots = np.linspace(domain[0], domain[1], num_elems + 1)
        knot_vector = np.array(dup(knots[0], degree) + list(knots) + dup(knots[-1], degree))
    else:
        knots = get_unique_knots(knot_vector)
    ext_ops = computeExtractionOperator( knot_vector, degree )
    ndofs = len(knot_vector) - degree - 1
    elems = []
    for i in range(len(ext_ops)):
        x_min = knots[i]
        x_max = knots[i+1]
        elems.append( Elem( x_min, x_max, range(i, i + degree+1), ext_ops[i] ) )
    K, F = AssembleL2( ndofs, elems, func, degree + 2 )
    d = np.linalg.solve(K,F)
    return elems, d

# I don't really understand what's going on here
def get_coeff_weights(degree, num_coeffs):
    # Since it's a bernstein basis the integrals work out in a straightforward manner
    B = BernsteinBasis(degree)
    b_vals = np.zeros(degree + 1)
    for i in range(degree + 1):
        b_vals[i] = B(0, i+1)
    weights = np.zeros(num_coeffs)
    stop_ind = len(b_vals)-num_coeffs + 1
    weights[0] = np.sum(b_vals[:stop_ind])
    for i in range(len(b_vals)-stop_ind):
        weights[1+i] = b_vals[stop_ind + i]
    return weights

# Get elements using the extraction operator
def get_elems( knot_vector, degree ):
    elems = []
    knots = get_unique_knots( knot_vector )
    knot_mults = get_knot_multiplicities( knot_vector )
    ext_ops = computeExtractionOperator( knot_vector, degree )
    elem_connect = np.arange(0, degree + 1)
    for i in range(len(ext_ops)):
        x_min = knots[i]
        x_max = knots[i+1]
        elems.append( Elem( x_min, x_max, copy.copy(elem_connect), ext_ops[i] ) )
        elem_connect += knot_mults[i+1]
    return elems

# THIS is the part from the projection paper I was thinking of earlier.
# (60) from page 72
def get_elem_func_weights( degree, knot_vector ):

    unique_knots = get_unique_knots( knot_vector )
    elems = get_elems( knot_vector, degree )
    num_elems = len(elems)
    funcs = spline_basis_funcs( degree, knot_vector )
    weights = [[] for i in range(num_elems)]
    func_integrals = []
    for func in funcs:
        func_sq = lambda x : func(x) ** 2
        func_integral, err = scipy.integrate.quad( func_sq, func.unique_knots[0], func.unique_knots[-1] )
        func_integrals.append(func_integral)
    for elem, ielem in zip(elems,range(len(elems))):
        for a in range(elem.globalFuncN):
            A = elem.elem_connect[a]
            func = funcs[A]
            func_sq = lambda x : func(x) ** 2
            elem_func_integral, err = scipy.integrate.quad( func_sq, elem.x_min, elem.x_max )
            func_integral = func_integrals[A]
            weights[ielem].append(elem_func_integral / func_integral )
    return weights

# helper function that integrates over elements
def integrate_funcs( elems, ndofs ):
    # wtf is ndofs supposed to mean
    func_integrals = np.zeros(ndofs)
    for elem in elems:
        for i in range(len(elem.elem_connect)):
            ifunc_g = elem.elem_connect[i]
            func_integrals[ifunc_g] += elem.integrate_func(i)
    return func_integrals

# square integration, that's all
def integrate_funcs_sq( elems, ndofs ):
    func_integrals = np.zeros(ndofs)
    for elem in elems:
        for i in range(len(elem.elem_connect)):
            ifunc_g = elem.elem_connect[i]
            func_integrals[ifunc_g] += elem.integrate_func_sq(i)
    return func_integrals

# no idea
def subinterval_op( min_val, max_val, degree ):
    M = np.zeros((degree + 1, degree + 1))
    for j in range(0,degree+1):
        for k in range(0,degree+1):
            for i in range(max(0, j+k-degree), min(j,k)+1):
                M[j,k] += BernsteinBasis(j)(max_val, i+1) * BernsteinBasis(degree - j)(min_val, k - i + 1)
    return M.T

# Get the matrix that performs fine-coarse projection, using sympy
class fine_coarse_projection_op(sympy.Function):
    nargs = 1
    @classmethod
    def eval(cls, p):
        x = sympy.Symbol('x')
        A = sympy.Matrix(2*p+2,2*p+2, lambda i,j:0)
        for i in range(2*p+2):
            for j in range(2*p+2):
                if i <= p and j <= p:
                    A[i,j] = sympy.integrate(bernstein_basis_poly(p, i, (2*x+1)) * bernstein_basis_poly(p, j, (2*x+1)), (x, -1, 0))
                elif i > p and j > p:
                    A[i,j] = A[i,j] + sympy.integrate(bernstein_basis_poly(p, i - p - 1, (2*x-1)) * bernstein_basis_poly(p, j - p - 1, (2*x-1)), (x, 0, 1))
        # if i==p and j == p:
        #     A[i,j] = 2*A[i,j]
        M = sympy.Matrix(p+1, 2 * p + 2, lambda i,j:0)
        for i in range(p+1):
            v = sympy.Matrix(2*p+2, 1, lambda i,j:0)
            for j in range(0, 2*p+2):
                if j <= p:
                    v[j,0] = sympy.integrate(bernstein_basis_poly(p, i, x) * bernstein_basis_poly(p, j, (2*x+1)), (x, -1, 0))
                elif j > p:
                    v[j,0] =  sympy.integrate(bernstein_basis_poly(p, i, x) * bernstein_basis_poly(p, j - p - 1, (2*x-1)), (x, 0, 1))
            # if j == p:
            #     v[j,0] = sympy.integrate(bernstein_basis_poly(p, i, x) * bernstein_basis_poly(p, p, (2*x+1)), (x, -1, 0)) + sympy.integrate(bernstein_basis_poly(p, i, x) * bernstein_basis_poly(p, j-p, (2*x-1)), (x, 0, 1))
            for j in range(0, 2*p+2):
                M[i,j] = (A.inv()*v)[j,0]
        return M

# This one doesn't change the dimension of the space so what's it doing
class coarse_projection_op(sympy.Function):
    nargs = 3
    @classmethod
    def eval(cls, p, local_min, local_max):
        x = sympy.Symbol('x')
        A = sympy.Matrix(p+1,p+1, lambda i,j:0)

        for i in range(p+1):
            for j in range(p+1):
                if i <= p and j <= p:
                    A[i,j] = sympy.integrate(bernstein_basis_poly(p, i, x) * bernstein_basis_poly(p, j, x), (x, local_min, local_max))
        M = sympy.Matrix(p+1, p + 1, lambda i,j:0)
        for i in range(p+1):
            v = sympy.Matrix(p+1, 1, lambda i,j:0)
            for j in range(0, p+1):
                if j <= p:
                    v[j,0] = sympy.integrate(bernstein_basis_poly(p, i, x) * bernstein_basis_poly(p, j, ( 2 * (x - local_min) / (local_max - local_min) - 1 )), (x, local_min, local_max))
            for j in range(0, p+1):
                M[i,j] = (A.inv()*v)[j,0]
        return np.array(np.array(M), np.float)

# Take a knot vector and get the bounds for all the elements. simple enough
def get_elem_bounds( kv ):
    ukv = get_unique_knots( kv )
    elem_bounds = []
    for i in range(len(ukv)-1):
        elem_bounds.append([ukv[i], ukv[i+1]])
    return elem_bounds

# Bounds is an (a,b) type thing. We're getting all the indices for the elements that
# intersect (a,b)
def get_isecting_elems( bounds, elems_bounds ):
    inds = []
    for test_bounds, i in zip(elems_bounds, itertools.count()):
        if (bounds[0] > test_bounds[0] and bounds[0] < test_bounds[1]) or equal( bounds[0], test_bounds[0]):
            inds.append(i)
        if (bounds[1] > test_bounds[0] and test_bounds[1] > bounds[0]) or equal( bounds[1], test_bounds[1]):
            inds.append(i)
        if test_bounds[0] > bounds[1]:
            break
    return list(set(inds))

# Get index-based maps from the elements of one to the other. 
# Knot vectors
def calc_elem_elem_maps( kv1, kv2 ):
    ukv1 = get_unique_knots( kv1 )
    ukv2 = get_unique_knots( kv2 )
    elem_elem_map = []
    rev_elem_elem_map = []
    elem_bounds1 = get_elem_bounds( ukv1 )
    elem_bounds2 = get_elem_bounds( ukv2 )
    for elem in elem_bounds1:
        elem_elem_map.append(get_isecting_elems(elem, elem_bounds2))
    for elem in elem_bounds2:
        rev_elem_elem_map.append(get_isecting_elems(elem, elem_bounds1))
    # rev is for reverse
    return elem_elem_map, rev_elem_elem_map

# Project the curve defined by knot_vector1 and degree1 and pts to the same curve, but defined
# by knot_vector2 and degree2 so we need the new pts
def project_pts( knot_vector1, degree1, knot_vector2, degree2, pts, local_pts = False, elem_elem_map = None ):

    new_pts = []
    elems1 = get_elems( knot_vector1, degree1 )
    elems2 = get_elems( knot_vector2, degree2 )
    num_elems1 = len(elems1)
    num_elems2 = len(elems2)
    refine = False
    coarsen = False
    if degree1 != degree2:
        M = degree_projection_op(degree1, degree2)
    else:
        M = np.eye(degree1 + 1)
    if elem_elem_map == None:
        if num_elems1 == num_elems2:
            if degree1 == degree2:
                M = np.eye(degree1 + 1)
            elem_elem_map = range(num_elems1)
        elif num_elems1 == 2 * num_elems2:
            coarsen = True
            if elem_elem_map == None:
                elem_elem_map = []
                for i in range(num_elems1):
                    elem_elem_map.append( (i+1) % 2 + (i+1) / 2 - 1 )
        elif num_elems2 == 2 * num_elems1:
            refine = True
            if elem_elem_map == None:
                elem_elem_map = []
                for i in range(num_elems1):
                    local_map = []
                    for j in range(2):
                        local_map.append(2 * i + j)
                    elem_elem_map.append(local_map)
    rev_elem_elem_map = [[] for i in range(num_elems2)]
    for i in range(num_elems1):
        try:
            for j in elem_elem_map[i]:
                rev_elem_elem_map[j] = i
        except TypeError:
            j = elem_elem_map[i]
            rev_elem_elem_map[j].append(i)
            if len(rev_elem_elem_map[j]) > 1:
                coarsen = True
    if not coarsen:
        for i in range(len(elems1)):
            if degree1 == degree2:
                M = np.eye(degree1 + 1)
            if local_pts:
                lpts = pts[i]
            else:
                lpts = [pts[j] for j in elems1[i].elem_connect]
            if degree1 != degree2 or not refine:
                proj_op = np.dot( np.dot( np.linalg.inv(elems2[i].elem_ext).T, M.T ), elems1[i].elem_ext .T )
                new_pts.append( np.dot( proj_op, lpts ) )
            if refine:
                for j in elem_elem_map[i]:
                    local_min = elems1[i].local_coord(elems2[j].x_min)
                    local_max = elems1[i].local_coord(elems2[j].x_max)
                    M = subinterval_op( local_min, local_max, degree1 )
                    proj_op = np.dot( np.dot( np.linalg.inv(elems2[j].elem_ext).T, M.T ), elems1[i].elem_ext .T )
                    new_pts.append( np.dot( proj_op, lpts ) )
    else:
        for i in range(num_elems2):
            # make empty arrays for the combined extraction op and projection op
            Ccomb = None
            A = np.array([])
            lpts = []
            for j in rev_elem_elem_map[i]:
                if local_pts:
                    lpts = lpts + pts[i]
                else:
                    lpts = lpts + [pts[k] for k in elems1[j].elem_connect]
                local_min = elems2[i].local_coord(elems1[j].x_min)
                local_max = elems2[i].local_coord(elems1[j].x_max)
                Asub = subinterval_op( local_min, local_max, degree1 )
                if Ccomb == None:
                    A = Asub
                    Ccomb = elems1[j].elem_ext
                else:
                    A = np.hstack([A,Asub])
                    Ccomb = la.block_diag(Ccomb, elems1[j].elem_ext)
            M = np.dot(A.T, np.linalg.inv(np.dot(A,A.T)))
            proj_op = np.dot( np.dot( np.linalg.inv(elems2[i].elem_ext).T, M.T ), Ccomb.T )
            new_pts.append( np.dot( proj_op, lpts ) )

    return new_pts

# % overlap between two elements
def overlap_frac( elem1, elem2 ):
    x_min = max(elem1.x_min, elem2.x_min)
    x_max = min(elem1.x_max, elem2.x_max)
    return (x_max - x_min) / float(elem2.x_max - elem2.x_min)

# Get the gramian, symbolically is all
def sympy_gramian(degree):
    class bernstein_basis(sympy.Function):
        nargs = 3
        @classmethod
        def eval(cls, n, i, x):
            return sympy.Piecewise( (0, x<-1), (0,x>1), (sympy.binomial(n,i) * ((x+1)/2) ** i * ((1-x)/2) **(n-i),True))
    class _gramian(sympy.Function):
        nargs = 1
        @classmethod
        def eval(cls, n):
            x = sympy.Symbol('x')
            A = sympy.Matrix(n+1,n+1, lambda i,j:0)
            for i in range(n+1):
                for j in range(n+1):
                    A[i,j] = sympy.integrate( bernstein_basis(n, i, x) * bernstein_basis(n, j, x), (x, -1, 1))
            return A
    return _gramian(degree)

# Get the gramian for a bernstein basis of given degree
def gramian(degree):

    class bernstein_basis(sympy.Function):
        nargs = 3
        @classmethod
        def eval(cls, n, i, x):
            return sympy.Piecewise( (0, x<-1), (0,x>1), (sympy.binomial(n,i) * ((x+1)/2) ** i * ((1-x)/2) **(n-i),True))
    class _gramian(sympy.Function):
        nargs = 1
        @classmethod
        def eval(cls, n):
            x = sympy.Symbol('x')
            A = sympy.Matrix(n+1,n+1, lambda i,j:0)
            for i in range(n+1):
                for j in range(n+1):
                    A[i,j] = sympy.integrate( bernstein_basis(n, i, x) * bernstein_basis(n, j, x), (x, -1, 1))
            return A
    return np.array(_gramian(degree)).astype('float')

# What is binom? Where is it coming from? I think it's scipy.stats.binom
# if it breaks you should def look into this
def inv_gramian(degree):
    A = np.zeros((degree+1,degree+1))
    for i in range(degree+1):
        for j in range(degree+1):
            bin_pi = binom(degree, i)
            bin_pj = binom(degree, j)
            for k in range(min(i,j)+1):
                bin_pkpmi = binom(degree+k+1, degree-i)
                bin_pkpmj = binom(degree+k+1, degree - j)
                bin_pmkpmi  = binom(degree-k, degree-i)
                bin_pmkpmj  = binom(degree-k, degree-j)
                A[i,j]+=(-1)**(i+j) / (bin_pi*bin_pj)*(2*k+1)*bin_pkpmi*bin_pkpmj*bin_pmkpmi*bin_pmkpmj
    # multiply by 0.5 because the above expression came from Juttler 1998 
    # on the dual basis over the unit interval and we're using the biunit interval
    return A*0.5

# Get the new points after projecting. Can handle degree projection and all of the 
# knot insertion, basis roughening and smoothing, etc.
# THIS FUNCTION IS THE WORKHORSE
def new_project_pts( knot_vector1, degree1, knot_vector2, degree2, pts, local_pts = False, cell_cell = False, coarsening = 'gramian' ):
    new_pts = []
    elems1 = get_elems( knot_vector1, degree1 )
    elems2 = get_elems( knot_vector2, degree2 )
    num_elems1 = len(elems1)
    num_elems2 = len(elems2)
    if degree1 != degree2:
        M = degree_projection_op(degree1, degree2)
        cell_cell = True
    else:
        M = np.eye(degree1 + 1)
        G = np.matrix(gramian(degree1))
        Ginv=np.linalg.inv(G)
    if cell_cell:
        elem_elem_map = []
        rev_elem_elem_map = []
        for i in range(len(elems1)):
            elem_elem_map.append([i])
            rev_elem_elem_map.append([i])
    else:
        elem_elem_map, rev_elem_elem_map = calc_elem_elem_maps(knot_vector1, knot_vector2)
    for i in range(num_elems2):
        # make empty arrays for the combined extraction op and projection op
        C = None
        A = np.array([])
        lpts = []
        for j in rev_elem_elem_map[i]:
            if local_pts:
                # Check whether this should be i or j
                lpts.append(pts[j])
            else:
                lpts = lpts + [pts[k] for k in elems1[j].elem_connect]
            print lpts
            if degree1 == degree2:
                # Need check for k-refinement
                # these are the [a,b] (or [\tilde{a},\tilde{a}])
                Ainter=np.matrix(np.eye(degree1+1))
                if len(rev_elem_elem_map[i])>1 and coarsening=='gramian':
                    x_min=max(elems1[j].x_min, elems2[i].x_min)
                    x_max=min(elems1[j].x_max, elems2[i].x_max)
                    local_min=elems2[i].local_coord(x_min)
                    local_max=elems2[i].local_coord(x_max)
                    inter_min=elems1[j].local_coord(x_min)
                    inter_max=elems1[j].local_coord(x_max)
                    Ainter=np.matrix(subinterval_op(inter_min, inter_max, degree1)).T
                else:
                    local_min = elems2[i].local_coord(elems1[j].x_min)
                    local_max = elems2[i].local_coord(elems1[j].x_max)
                # local_min = elems1[j].local_coord(elems2[i].x_min)
                # local_max = elems1[j].local_coord(elems2[i].x_max)
                # elem_frac is phi_i in the paper.
                if cell_cell:
                    elem_frac = 1.0
                else:
                    elem_frac = overlap_frac( elems1[j], elems2[i] )
                Asub = subinterval_op( local_min, local_max, degree1 )
                if coarsening=='least_sq':
                    if C == None:
                        A = Asub * elem_frac
                        C = elems1[j].elem_ext * elem_frac
                    else:
                        A = np.hstack([A,Asub * elem_frac])
                        C = la.block_diag(C, elems1[j].elem_ext * elem_frac)
                    if not A.shape[0] == A.shape[1]:
                        M = np.dot(A.T, np.linalg.inv(np.dot(A,A.T)))
                    else:
                        M = np.linalg.inv(A)
                elif coarsening=='avg':
                    if C == None:
                        Mt = np.linalg.inv(Asub).T
                        C = elems1[j].elem_ext * elem_frac
                    else:
                        Mt = np.hstack([Mt, np.linalg.inv(Asub).T])
                        C = la.block_diag(C, elems1[j].elem_ext * elem_frac)
                elif coarsening=='gramian':
                    if len(rev_elem_elem_map[i])>1:
                        Ai=np.matrix((Asub))
                        if C == None:
                            Mt = np.array(Ginv * Ai * G * Ainter)
                            C = elems1[j].elem_ext * elem_frac
                        else:
                            Mt = np.hstack([Mt, np.array(Ginv*Ai*G*Ainter)])
                            C = la.block_diag(C, elems1[j].elem_ext * elem_frac)
                    else:
                        Ai=np.linalg.inv(Asub).T
                        if C == None:
                            Mt = Ai
                            C = elems1[j].elem_ext * elem_frac
                        else:
                            Mt = np.hstack([Mt, Ai])
                            C = la.block_diag(C, elems1[j].elem_ext * elem_frac)
                M = Mt.T
            else:
                C = elems1[j].elem_ext
    #print "C.T"
    #print C.T
    #print M.T
    #print np.linalg.inv(elems2[i].elem_ext).T
        proj_op = np.dot( np.dot( np.linalg.inv(elems2[i].elem_ext).T, M.T ), C.T )
        new_pts.append( np.dot( proj_op, lpts ) )
    #print( proj_op )
    #print lpts
    #print new_pts
    return new_pts

# Pretty self-explanatory?
def compute_func_weights( ext_info ):
    wts = np.zeros(ext_info['nodeN'])
    for ifunc in range(ext_info['nodeN']):
        ielems = ext_info['func_map'][ifunc]['elems']
        ifuncs_l = ext_info['func_map'][ifunc]['local_ifuncs']
        for ielem, ifunc_l in zip(ielems, ifuncs_l):
            elem = BElem(ext_info['spatial_dim'], ext_info['param_dim'], ext_info['elems'][ielem])
            wts[ifunc] += elem.param_vol() * np.sum(elem.elem_info['ext_op'][ifunc_l])
    return wts

# like overlap_frac, but can handle more than one dimension
def new_overlap_frac(elem1, elem2):
    inter_ranges = []
    overlap_vol = 1
    for i in range(elem1.param_dim):
        s_min = max(elem1.bounds(i)[0], elem2.bounds(i)[0])
        s_max = min(elem1.bounds(i)[1], elem2.bounds(i)[1])
        overlap_vol = s_max - s_min
    return overlap_vol / elem2.param_vol()

# Convenience function? I'm assuming. Not sure what the point of it is.
def update_elem_bounds(ext_info):
    for i in range(len(ext_info['elems'])):
        elem = BElem(ext_info['spatial_dim'], ext_info['param_dim'], ext_info['elems'][i])
        ext_info['elems'][i]['bounds'] = elem.bounds()

# get_isecting_elems, equipped for multiple dimensions
def new_get_isecting_elems( bounds, elem_infos ):
    inds = []
    for elem_info, i in zip(elem_infos, itertools.count()):
        overlap = True
        for idim in range(len(elem_info['bounds'])):
            test_bounds = bounds[idim]
            dim_bounds = elem_info['bounds'][idim]
            if equal(test_bounds[0], dim_bounds[1]) or test_bounds[0] > dim_bounds[1]:
                overlap = False
                break
            elif equal(test_bounds[1], dim_bounds[0]) or test_bounds[1] < dim_bounds[0]:
                overlap = False
                break
        if overlap:
            inds.append(i)
    return list(set(inds))

# More class-based than calc_elem_elem_maps. 
def new_calc_elem_elem_maps(ext_info1, ext_info2):
    update_elem_bounds(ext_info1)
    update_elem_bounds(ext_info2)
    elem_elem_map = []
    rev_elem_elem_map = []
    for elem_info in ext_info1['elems']:
        elem_elem_map.append(new_get_isecting_elems(elem_info['bounds'], ext_info2['elems']))
    for elem_info in ext_info2['elems']:
        rev_elem_elem_map.append(new_get_isecting_elems(elem_info['bounds'], ext_info1['elems']))
    return elem_elem_map, rev_elem_elem_map

# Kronecker product of input matrices (tensor product matrix)
def kron_all(*mats):
    """
    Take the Kronecker product of all the input matrices.  The input order
    corresponds to the order in which the product will be computed.
    """

    ret_mat = np.ones(1)
    mat_list = list(mats)
    mat_list.reverse()
    for M in mat_list:
        ret_mat = np.kron(M,ret_mat)
    return ret_mat

# Class-based, multidimensional version of new_project_pts
def multid_project_pts( ext_info1, ext_info2, pts, local_pts = False, elem_elem_map = None, cell_cell=False ):

    new_pts = []
    M_list = []
    for d1, d2 in zip(ext_info1['elems'][0]['degrees'],ext_info2['elems'][0]['degrees']):
        if d1 != d2:
            M_list.append(degree_projection_op(d1, d2))
            elevate=True
        else:
            M_list.append(np.eye(d1 + 1))
            elevate=False
    M_list.reverse()
    M = kron_all(*M_list)
    if cell_cell:
        elem_elem_map = []
        rev_elem_elem_map = []
        for i in range(len(ext_info1["elems"])):
            elem_elem_map.append([i])
            rev_elem_elem_map.append([i])
    else:
        elem_elem_map, rev_elem_elem_map = new_calc_elem_elem_maps(ext_info1, ext_info2)
    for i in range(len(ext_info2["elems"])):
        elem2_info = ext_info2['elems'][i]
        elem2 = BElem(ext_info2['spatial_dim'], ext_info2['param_dim'], elem2_info)
        new_bpts=[]
        for j in rev_elem_elem_map[i]:
            elem1_info = ext_info1['elems'][j]
            elem1 = BElem(ext_info1['spatial_dim'], ext_info1['param_dim'], elem1_info)
            Ainter_list = []
            local_min=np.zeros(ext_info1['param_dim'])
            local_max=np.zeros(ext_info1['param_dim'])
            for idim in range(ext_info1['param_dim']):
                degree1 = elem1_info['degrees'][idim]
                Ainter_list.append(np.matrix(np.eye(degree1+1)))
                x_min=max(elem1.bounds(idim)[0], elem2.bounds(idim)[0])
                x_max=min(elem1.bounds(idim)[1], elem2.bounds(idim)[1])
                pt=np.zeros(ext_info1['param_dim'])
                pt[idim]=x_min
                local_min[idim]=elem2.local_coord(pt)[idim]
                inter_min=elem1.local_coord(pt)[idim]
                pt[idim]=x_max
                print x_min, x_max
                local_max[idim]=elem2.local_coord(pt)[idim]
                inter_max=elem1.local_coord(pt)[idim]
                Ainter_list[idim]=np.matrix(subinterval_op(inter_min, inter_max, degree1)).T
            Ainter_list.reverse()
            Ainter=kron_all(*Ainter_list)
            print i,j
            print local_min
            print local_max
            if local_pts:
                lpts = pts[i]
            else:
                lpts = [pts[k] for k in elem1_info['elem_connect']]
            bpts=np.dot(elem1.ext_op().T, lpts)
            # localize points to the overlap
            bpts=np.dot(Ainter, bpts)
            if elevate:
                new_bpts=np.dot(M.T, bpts)
            else:
                if cell_cell:
                    elem_frac = 1.0
                else:
                    elem_frac = new_overlap_frac( elem1, elem2 )
                gramian_list=[]
                inv_gramian_list=[]
                A_list=[]
                A_inv_list=[]
                for idim in range(ext_info2['param_dim']):
                    degree1 = elem1_info['degrees'][idim]
                    A_sub = subinterval_op( local_min[idim], local_max[idim], degree1 )
                    A_list.append(A_sub)
                    A_inv_list.append(np.linalg.inv(A_sub))
                    gramian_list.append(gramian(degree1))
                    inv_gramian_list.append(inv_gramian(degree1))
                A_list.reverse()
                A_inv_list.reverse()
                gramian_list.reverse()
                inv_gramian_list.reverse()
                G=np.matrix(kron_all(*gramian_list))
                Ginv=np.matrix(kron_all(*inv_gramian_list))
                A=np.matrix(kron_all(*A_list))
                Ainv=np.matrix(kron_all(*A_inv_list))
                temp_bpts=np.dot(elem_frac*Ginv*Ainv.T*G, bpts)
                if len(new_bpts)==0:
                    new_bpts=temp_bpts[:]
                else:
                    new_bpts=[new_bpts[i]+temp_bpts[i] for i in range(len(new_bpts))]
            new_pts.append(np.dot(np.linalg.inv(elem2.ext_op()).T, new_bpts))
    return new_pts

# Use the weighting scheme to smooth out the projected pieces
def smooth_pts( knot_vector, degree, pts, sq_weight = False ):
    elems = get_elems( knot_vector, degree )
    ndofs = len(knot_vector) - degree - 1
    if sq_weight:
        function_integrals = integrate_funcs_sq( elems, ndofs )
    else:
        function_integrals = integrate_funcs( elems, ndofs )
    new_pts = []
    w_check = np.zeros(ndofs)\

    for i in range(len(elems)):
        for a in range(elems[i].globalFuncN):
            A = elems[i].elem_connect[a]
            if sq_weight:
                coeff_weight = elems[i].integrate_func_sq(a) / function_integrals_sq[A]
            else:
                coeff_weight = elems[i].integrate_func(a) / function_integrals[A]
            try:
                new_pts[A] += pts[i][a] * coeff_weight
            except IndexError:
                new_pts.append(pts[i][a] * coeff_weight)
            w_check[A] += coeff_weight
    return new_pts

# Something about picking maximum weights
def select_pts( knot_vector, degree, pts, sq_weight = False ):
    elems = get_elems( knot_vector, degree )
    ndofs = len(knot_vector) - degree - 1
    if sq_weight:
        function_integrals = integrate_funcs_sq( elems, ndofs )
    else:
        function_integrals = integrate_funcs( elems, ndofs )
    new_pts = [None for A in range(ndofs)]
    w_check = np.zeros(ndofs)

    max_weights = np.zeros(ndofs)
    max_elems = np.zeros(ndofs, dtype = int)

    for i in range(len(elems)):
        for a in range(elems[i].globalFuncN):
            A = elems[i].elem_connect[a]
            if sq_weight:
                coeff_weight = elems[i].integrate_func_sq(a) / function_integrals_sq[A]
            else:
                coeff_weight = elems[i].integrate_func(a) / function_integrals[A]
            if coeff_weight > max_weights[A]:
                max_weights[A] = coeff_weight
                max_elems[A] = i
    for i in range(len(elems)):
        for a in range(elems[i].globalFuncN):
            A = elems[i].elem_connect[a]
            if max_elems[A] == i:
                new_pts[A] = pts[i][a]
    return new_pts

# smooth_pts, equipped for multiple dimensions
def multid_smooth_pts(ext_info, pts):
    new_pts = {}
    wts = compute_func_weights(ext_info)
    for elem_info,i in zip(ext_info['elems'],itertools.count()):
        epts = pts[i]
        elem = BElem(ext_info['spatial_dim'], ext_info['param_dim'], elem_info)
        for a in range(elem_info['global_func_n']):
            A = elem_info['elem_connect'][a]
            smooth_weight = sum(elem_info['ext_op'][a]) * elem.param_vol() / wts[A]
            try:
                new_pts[A] += epts[a] * smooth_weight
            except KeyError:
                new_pts[A] = epts[a] * smooth_weight
    return new_pts.values()

def legendre_project(x_min, x_max, degree, func):
    gauss_rule = degree + 1
    gpts, wts = gauss_integrate.get1dGaussRule( gauss_rule )
    def local_func(s):
        x = ( s + 1 )* (x_max - x_min) / 2.0 + x_min
        return func(x)
    leg_coeffs = [(2*n+1)/float(2) for n in range(degree+1)]
    Pm2 = np.ones(gpts.shape)
    f = local_func(gpts)
    leg_coeffs[0] *= np.dot(Pm2*wts,f)
    Pm1 = gpts
    leg_coeffs[1] *= np.dot(Pm1*wts,f)
    for n in range(2,degree+1):
        P = ((2*n - 1)*gpts*Pm1-(n-1)*Pm2)/float(n)
        Pm2 = Pm1
        Pm1 = P
        leg_coeffs[n] *= np.dot(P*wts,f)
    return leg_coeffs
    
def legendre_bernstein_mat(degree):

    M = np.zeros((degree+1,degree+1))
    for j in range(degree+1):
        bin_nj = scipy.special.binom(degree,j)
        for k in range(degree+1):
            i_min = max(0, j + k - degree)
            i_max = min(j,k)
            for i in range(i_min, i_max+1):
                bin_ki = scipy.special.binom(k,i)
                bin_nk_ji = scipy.special.binom(degree-k,j-i)
                M[j,k] += 1 / bin_nj * (-1)**(k+i) * bin_ki ** 2 * bin_nk_ji
    return M

def bernstein_legendre_mat(degree):

    M = np.zeros((degree+1,degree+1))
    for j in range(degree+1):
        bin_npj = scipy.special.binom(degree+j,j)
        den = bin_npj * (degree+j+1)
        sqrj = np.sqrt(2*j+1)
        for k in range(degree+1):
            i_min = max(0, j + k - degree)
            i_max = min(j,k)
            for i in range(0, j+1):
                bin_kpi = scipy.special.binom(k+i,k)
                bin_ji = scipy.special.binom(j,i)
                bin_nkji = scipy.special.binom(degree-k+j-i,degree-k)
                M[j,k] += (2*j+1) / den * (-1)**(j+i) * bin_ji * bin_kpi * bin_nkji
    return M

def legendre_bernstein(leg_coeffs):

    degree = len(leg_coeffs) - 1
    M = legendre_bernstein_mat(degree)
    return np.dot(M,leg_coeffs)

def approx_l2_project(func, domain, degree, num_elems=10, knot_vector=None, sq_weight=False, smoothing_method='weight', use_legendre=False):

    if knot_vector == None:
        knots = np.linspace(domain[0], domain[1], num_elems + 1)
        knot_vector = np.array(dup(knots[0], degree) + list(knots) + dup(knots[-1], degree))
    else:
        knots = get_unique_knots(knot_vector)
    
    ext_ops = computeExtractionOperator(knot_vector, degree)
    ndofs = len(knot_vector) - degree - 1
    d = np.zeros(ndofs)
    w_check = np.zeros(ndofs)
    num_shared = np.ones(ndofs) * min(degree + 1, num_elems)
    
    for i in range(min(degree+1,num_elems)):
        num_shared[i] = i + 1
    for i in range(min(degree+1,num_elems)):
        num_shared[-i-1] = i + 1
    
    elems = get_elems( knot_vector, degree )
    gauss_rule = degree + 1
    gpts, wts = gauss_integrate.get1dGaussRule( gauss_rule )
    d_local = []
    
    if sq_weight:
        function_integrals = integrate_funcs_sq( elems, ndofs )
    else:
        function_integrals = integrate_funcs( elems, ndofs )
    
    leg_bern = legendre_bernstein_mat(degree)
    bern_leg = bernstein_legendre_mat(degree)

    for i in range(len(ext_ops)):
        x_min = knots[i]
        x_max = knots[i+1]

        if use_legendre:
            leg_coeffs = legendre_project(x_min, x_max, degree, func)
            bez_coeffs = np.dot(leg_bern, leg_coeffs)
            local_coeffs = np.linalg.solve(ext_ops[i].T,bez_coeffs)

        else:
            temp_elem = Elem( x_min, x_max, range(i, i + degree+1), np.eye(ext_ops[i].shape[0]))
            ke, fe = LocalAssemblyL2(temp_elem, func, gpts, wts )
            bez_coeffs = np.linalg.solve(ke,fe.T)
            local_coeffs = np.linalg.solve(ext_ops[i].T,bez_coeffs)
            # Additional code from alternate version of function
            # leg_coeffs = legendre_project(x_min, x_max, degree, func)
            # bez_coeffs = legendre_bernstein(leg_coeffs) #np.dot(leg_bern, leg_coeffs)
            # local_coeffs = np.linalg.solve(ext_ops[i].T,bez_coeffs)
            # local_coeffs = np.linalg.solve(ext_ops[i].T,bez_coeffs)

        d_local.append(local_coeffs)
        # More additional code from alternate version of function
        # for a in range(elems[i].globalFuncN):
        #     A = elems[i].elem_connect[a]
        #     if sq_weight:
        #         coeff_weight = elems[i].integrate_func_sq(a) / function_integrals_sq[A]
        #     else:
        #         coeff_weight = elems[i].integrate_func(a) / function_integrals[A]
        #     d[A] += local_coeffs[a] * coeff_weight
        #     w_check[A] += coeff_weight

    if smoothing_method == 'weight':
        d = smooth_pts(knot_vector, degree, d_local)
    elif smoothing_method == 'select':
        d = select_pts(knot_vector, degree, d_local)
    else:
        raise(ValueError('Invalid input value for approx_l2_project'))
    return elems, d, d_local

def multid_local_l2_project(func, ext_info):
    d_local = []
    func_cache = {}
    for elem_info in ext_info['elems']:
        elem = BElem(ext_info['spatial_dim'], ext_info['param_dim'], elem_info, func_cache = func_cache)
        elem.elem_info['ext_op'] = np.eye(elem.elem_info['ext_op'].shape[0])
        k, f = LocalAssemblyL2New( elem, func )
        Q_local = np.linalg.solve(k,f)
        P_local = np.dot( np.linalg.inv( elem_info['ext_op'].T ), Q_local)
        d_local.append(P_local)
    d = multid_smooth_pts(ext_info, d_local)
    return d, d_local

def new_color_cycle( num_colors, cmap = None, clims = (-0.1,1.1) ):

    if cmap == None:
        
        cmap = matplotlib.cm.ScalarMappable(cmap = 'spectral')
    cmap.set_clim(-0.1,1.05)
    color_list = list(cmap.to_rgba(np.linspace(0,1,num_colors)))
    color_list.reverse()
    color_cycle = itertools.cycle(color_list)
    return color_cycle

def convergence_compare(filename, f = None, min_elems = 2, max_elems = 10, elem_step = 4, min_degree = 2, max_degree = 5, smoothing_method = 'weight' ):
    
    out_file = h5py.File(filename)
    if f == None:
        def f(x):
            return np.sin(2*np.pi*x)
        # return 1 - 4*(0.5 - x)**2.0 #


    degs = range(min_degree,max_degree + 1)
    num_degs = len(degs)
    num_elems = np.array(range(min_elems, max_elems,elem_step))

    elem_size = 1 / num_elems.astype(float)
    out_file['degs'] = degs
    for deg in degs:
        group_name = "p=" + str(deg)
        l2_error = np.zeros(len(num_elems))
        approx_l2_error = np.zeros(len(num_elems))
        gauss_rule = deg + 1
        for i in range(len(num_elems)):
            n = num_elems[ i ]
            elems, d = l2_project(f, [0,1], deg, n)
            l2_error[i] = L2Error( elems, d, f, gauss_rule)
            elems_approx, d_approx, d_local = approx_l2_project(f, [0,1], deg, n, smoothing_method = smoothing_method)
            approx_l2_error[i] = L2Error( elems_approx, d_approx, f, gauss_rule)
        out_file[group_name + "/elem_size"] = elem_size
        out_file[group_name + "/l2_error"] = l2_error
        out_file[group_name + "/approx_l2_error"] = approx_l2_error

def extract_cell(knot_vector, degree, icell):
    cell_ext_op = []
    unique_knots = get_unique_knots(knot_vector)
    knot_mults = get_knot_multiplicities(knot_vector)
    iknot = sum(knot_mults[:icell])
    for i in range(degree + 1):
        start_ind = iknot - degree + i
        stop_ind = start_ind + degree + 2
        local_kv = knot_vector[start_ind:stop_ind]
        seg_rows = local_knots_to_segment_rows(local_kv, filter_zlen = False)
        cell_ext_op.append(seg_rows[degree - i])
    return cell_ext_op

def cpt_bezier_plot(elems, d_x_local, d_y_local, ax, color_cycle = None, num_pts = 100, markersize = 10, **kwargs):

    num_elems = len(elems)
    if color_cycle == None:
        color_cycle = new_color_cycle( num_elems )
    else:
        color_cycle = color_cycle
    for i in range(num_elems):
        elem = elems[i]
        elem_color = color_cycle.next()

        b_kwargs = kwargs.copy()
        b_kwargs["linewidth"] = 1
        b_kwargs["c"] = "k"
        b_kwargs["markerfacecolor"] = elem_color
        b_kwargs["markersize"] = markersize
        elem_ext = elem.elem_ext
        bpt_x = np.dot(elem_ext.T, d_x_local[i])
        bpt_y = np.dot(elem_ext.T, d_y_local[i])
        ax.plot(bpt_x, bpt_y, **b_kwargs)

def extend_knot_vector(knot_vector, step, new_mult = 1):
    unique_knots = get_unique_knots(knot_vector)
    knot_mults = get_knot_multiplicities(knot_vector)
    new_vector = []
    for i in range(len(unique_knots) - 1):
        for m in range(knot_mults[i]):
            new_vector.append(unique_knots[i])
    for m in range(new_mult):
        new_vector.append(unique_knots[-1])
    for m in range(knot_mults[-1]):
        new_vector.append(unique_knots[-1] + step)
    return new_vector

def bezier_seg_plot(elems, d_x_local, d_y_local, ax, color_cycle = None, num_pts = 100, markersize = 10, **kwargs):


    num_elems = len(elems)
    if color_cycle == None:
        color_cycle = new_color_cycle( num_elems )
    else:
        color_cycle = color_cycle
    for i in range(num_elems):
        elem = elems[i]
        elem_color = color_cycle.next()
        t, x_sol = elem.eval_sol(d_x_local[i], num_pts = num_pts)
        t, y_sol = elem.eval_sol(d_y_local[i], num_pts = num_pts)
        a_kwargs = {"c" : elem_color}
        a_kwargs.update(kwargs)
        a_kwargs["linewidth"] = 2
        ax.plot(x_sol, y_sol, **a_kwargs)

def cpt_local_plot(elems, d_x_local, d_y_local, ax, color_cycle = None, markersize = 10, **kwargs):


    num_elems = len(d_x_local)
    if color_cycle == None:
        color_cycle = new_color_cycle( num_elems )
    else:
        color_cycle = color_cycle
    for i in range(num_elems):
        elem_color = color_cycle.next()
        b_kwargs = kwargs.copy()
        b_kwargs["linewidth"] = 1
        if not(kwargs.has_key("marker")):
            b_kwargs["marker"] = "o"
        b_kwargs["c"] = "k"
        b_kwargs["markerfacecolor"] = elem_color
        b_kwargs["markersize"] = markersize
        ax.plot(d_x_local[i], d_y_local[i], **b_kwargs)

def cpt_global_plot_ratio(elems, d_x, d_y, ax, color_cycle = None, markersize = 40, linewidth = 1, sq_weight = False, **kwargs):

    num_elems = len(elems)
    if color_cycle == None:
        color_cycle = new_color_cycle( num_elems )
    else:
        color_cycle = color_cycle

    num_funcs = elems[-1].elem_connect[-1] + 1
    if sq_weight:
        function_integrals = integrate_funcs_sq( elems, num_funcs )
    else:
        function_integrals = integrate_funcs( elems, num_funcs )
    func_ratios = [[] for i in range(num_funcs)]
    markers = [[] for i in range(num_funcs)]
    func_elem_colors = [[] for i in range(num_funcs)]
    for i in range(num_elems):
        elem = elems[i]
        elem_color = color_cycle.next()
        for a in range(elem.globalFuncN):
            A = elem.elem_connect[a]
            if sq_weight:
                func_ratios[A].append(elem.integrate_func_sq(a) / function_integrals_sq[A])
            else:
                func_ratios[A].append(elem.integrate_func(a) / function_integrals[A])
            func_elem_colors[A].append(elem_color)
    ax.plot(d_x, d_y, linewidth = linewidth, color = "k", zorder = 0, **kwargs)
    for i in range(num_funcs):
        markers[i] = ratio_markers( func_ratios[i] )
        for j in range(len(func_ratios[i])):
            if len(func_ratios[i])==1:
                lw=0
            else:
                lw=0.3
            ax.scatter(d_x[i], d_y[i], marker = markers[i][j], color = func_elem_colors[i][j], s = markersize, zorder = 20, c = "k", edgecolor='w', linewidths = lw, **kwargs)

def cpt_global_plot(d_x = None, d_y = None, pts = None, ax = None, markersize = 10, marker = 'o', color_cycle = None, **kwargs):
    
    if ax == None:
        ax = matplotlib.pyplot.gca()
    if pts == None and ( d_x == None or d_y == None ):
        raise ValueError( 'Either pts or d_x and d_y must not be equal to None')
    if pts != None:
        d_x = np.array(pts)[:,0]
        d_y = np.array(pts)[:,1]
    p, = ax.plot(d_x, d_y, **kwargs )
    zorder = p.get_zorder()
    kwargs['zorder'] = zorder + 1
    if color_cycle == None:
        try:
            colors = kwargs['color']
        except AttributeError:
            colors = 'k'
    else:
        colors = []
        try:
            kwargs.pop('color')
            kwargs.pop('c')
        except:
            pass
        for x in d_x:
            colors.append(color_cycle.next())
    for x, y, i in zip(d_x, d_y, itertools.count()):
        try:
            c = colors[i]
        except:
            c = colors
        ax.plot(x, y, marker, markersize = markersize, c = c, marker = marker, markeredgecolor = 'None', **kwargs)

def global_seg_plot(elems, d_x, d_y, ax = None, color_cycle = None, linewidth = 1, num_pts = 100, **kwargs):
    if ax == None:
        ax = plt.gca()
    num_elems = len(elems)
    if color_cycle == None:
        color_cycle = new_color_cycle( num_elems )
    else:
        color_cycle = color_cycle

    for i in range(num_elems):
        elem = elems[i]
        elem_color = color_cycle.next()
        nkwargs = kwargs.copy()
        if not(nkwargs.has_key("color") or nkwargs.has_key("c")):
            nkwargs["color"] = elem_color
        d_x_elem = elem.localize(d_x)
        d_y_elem = elem.localize(d_y)
        t, x_sol = elem.eval_sol(d_x_elem, num_pts = num_pts)
        t, y_sol = elem.eval_sol(d_y_elem, num_pts = num_pts)
        ax.plot(x_sol, y_sol, linewidth = linewidth, **nkwargs)

def cpt_example(x, y, degree, num_elems, ax):
    elems_approx, d_x_approx, d_x_local = approx_l2_project(x, [0,1], degree, num_elems)
    elems_approx, d_y_approx, d_y_local = approx_l2_project(y, [0,1], degree, num_elems)
    # elems, d_x = l2_project(x, [0,1], degree, num_elems)
    # elems, d_y = l2_project(y, [0,1], degree, num_elems)
    # global_seg_plot( elems, d_x, d_y, ax, linewidth = 2 )
    # ax.plot(d_x_approx, d_y_approx)
    # ax.scatter(d_x_approx, d_y_approx)
    color_cycle = itertools.cycle(['b', 'orange', 'g'])
    bezier_seg_plot( elems_approx, d_x_local, d_y_local, ax, linestyle = "--", color_cycle = color_cycle )
    # cpt_local_plot( elems_approx, d_x_local, d_y_local, ax )
    global_seg_plot( elems_approx, d_x_approx, d_y_approx, ax, linewidth = 2, color_cycle = color_cycle )
    cpt_global_plot_ratio( elems_approx, d_x_approx, d_y_approx, ax, color_cycle = color_cycle )

def cpt_plot(w = 4, h = 1, base_name = "../pics/local_proj_steps", filetype = "pdf"):

    def x(t):
        return 2 *np.sqrt( t ) ** 3

    def y(t):
        return 0.1 * np.sin(3 * np.pi * t)

    # fig, axs = plt.subplots(5,1,sharey=True, sharex=True)
    color_cycle = itertools.cycle(color_blind_list()[0:3])
    num_figs = 5
    figs = []
    axs = []

    for i in range(num_figs):
        fig = plt.figure(frameon=False)
        fig.set_size_inches(w,h)
        if len(axs) > 0:
            share = axs[-1]
        else:
            share = None
        ax = plt.Axes(fig, [0., 0., 1., 1.], sharex = share, sharey = share)
        ax.set_axis_off()

        fig.add_axes(ax)
        figs.append(fig)
        axs.append(ax)
    for ax in axs:
        t = np.linspace(0,1)
        ax.plot(x(t),y(t), '--k')
        # ax.axis("off")
    axs[0].plot(x(t),y(t), '-k')
    degree = 2
    num_elems = 3
    kv=[0,0,0,0.5,0.75,1,1,1]
    kv=[0,0,0,1./3,2./3,1,1,1]
    elems_approx, d_x_approx, d_x_local = approx_l2_project(x, [0,1], degree, knot_vector=kv)
    elems_approx, d_y_approx, d_y_local = approx_l2_project(y, [0,1], degree, knot_vector=kv)
    bezier_seg_plot( elems_approx, d_x_local, d_y_local, axs[1], linestyle = "-", color_cycle = color_cycle )
    cpt_bezier_plot( elems_approx, d_x_local, d_y_local, axs[1], linestyle = "-", color_cycle = color_cycle, marker = 's', markersize = 5 )
    cpt_local_plot( elems_approx, d_x_local, d_y_local, axs[2], color_cycle = color_cycle, marker = "v", markersize = 5 )
    bezier_seg_plot( elems_approx, d_x_local, d_y_local, axs[2], linestyle = "-", color_cycle = color_cycle )
    global_seg_plot( elems_approx, d_x_approx, d_y_approx, axs[3], linewidth = 2, color_cycle = color_cycle )
    cpt_global_plot_ratio( elems_approx, d_x_approx, d_y_approx, axs[3], color_cycle = color_cycle )
    global_seg_plot( elems_approx, d_x_approx, d_y_approx, axs[4], linewidth = 2, color = color_blind_list()[3], zorder = -10 )

    elems, d_x = l2_project(x, [0,1], degree, num_elems)
    elems, d_y = l2_project(y, [0,1], degree, num_elems)
    # global_seg_plot( elems, d_x, d_y, axs[4], linewidth = 2, color = "MediumPurple" )
    # cpt_global_plot( elems, d_x, d_y, axs[4], color_cycle = color_cycle )
    fig_num = 1
    x_pad = 0.05
    y_pad = 0.05
    axs[0].set_xlim(min(np.amin(d_x), np.amin(np.array(d_x_local))) - x_pad, max(np.amax(d_x), np.amax(np.array(d_x_local))) + x_pad)
    axs[0].set_ylim(min(np.amin(d_y), np.amin(np.array(d_y_local))) - y_pad, max(np.amax(d_y), np.amax(np.array(d_y_local))) + y_pad)
    # axs[0].axis("equal")

    for fig in figs:
        fig.savefig(base_name + str(fig_num) + "." + filetype)
        fig.show()
        fig_num += 1

def eval_curve( kv, degree, pts, num_pts = 100 ):
    funcs = spline_basis_funcs( degree, kv )
    t = np.linspace(kv[0],kv[-1], num_pts)
    w = np.zeros(t.shape)
    if len(pts[0])==2:
        w += 1;
    x = np.zeros(t.shape)
    y = np.zeros(t.shape)
    for pt, func in zip(pts, funcs):
        f_vals = func(t)
        # ax.plot(t,f_vals)
        if len(pts[0]) == 3:
            x += pt[0] * f_vals * pt[2]
            y += pt[1] * f_vals * pt[2]
            w += pt[2] * f_vals
        else:
            x += pt[0] * f_vals
            y += pt[1] * f_vals
    return t,x,y,w

def deriv_pts(kv, degree, pts):
    dpts = []
    for i in range(len(pts)-1):
        dpts.append(degree / (kv[i+degree+1] - kv[i+1]) * (np.array(pts[i+1]) - np.array(pts[i])))
    dkv = list(kv)
    dkv.pop(0)
    dkv.pop(-1)
    return dkv, dpts

def spline_func(kv,degree,pts):
    funcs = spline_basis_funcs( degree, kv )
    def eval(t):
        x = np.zeros_like(t)
        y = np.zeros_like(t)
        for pt, func in zip(pts, funcs):
            f_vals = func(t)
            x += pt[0] * f_vals
            y += pt[1] * f_vals
        return x,y
    return eval

def tangent_func(kv,degree,pts):
    dkv, dpts = deriv_pts(kv, degree, pts)
    dfuncs = spline_basis_funcs( degree - 1, dkv )
    def eval(t):
        u = np.zeros_like(t)
        v = np.zeros_like(t)
        for dpt, dfunc in zip(dpts, dfuncs):
            df_vals = dfunc(t)
            u += dpt[0] * df_vals
            v += dpt[1] * df_vals
        return u,v
    return eval

def unit_tangent_func(kv,degree,pts):
    tfunc = tangent_func(kv,degree,pts)
    def eval(t):
        u,v = tfunc(t)
        mag = np.sqrt(u**2+v**2)
        return u / mag, v / mag
    return eval

def unit_normal_func(kv,degree,pts):
    utfunc = unit_tangent_func(kv,degree,pts)
    def eval(t):
        u,v = utfunc(t)
        return -v, u
    return eval

def eval_tangent(kv, degree, pts, num_pts=10):
    t = np.linspace(kv[0],kv[-1], num_pts)
    func = spline_func(kv,degree,pts)
    x,y = func(t)
    tfunc = tangent_func(kv,degree,pts)
    u,v = tfunc(t)
    return t,x,y, u, v

def eval_unit_tangent(kv, degree, pts, num_pts=10):
    t = np.linspace(kv[0],kv[-1], num_pts)
    func = spline_func(kv,degree,pts)
    x,y = func(t)
    utfunc = unit_tangent_func(kv,degree,pts)
    u,v = utfunc(t)
    return t,x,y, u, v

def eval_unit_normal(kv, degree, pts, num_pts=10):
    t = np.linspace(kv[0],kv[-1], num_pts)
    func = spline_func(kv,degree,pts)
    x,y = func(t)
    unfunc = unit_normal_func(kv,degree,pts)
    u,v = unfunc(t)
    return t,x,y, u, v

def plot_curve( kv, degree, pts, num_pts = 100, t_range=None, ax = None, **kwargs ):
    if ax == None:
        ax = plt.gca()
    t,x,y,w = eval_curve(kv, degree, pts, num_pts)
    if t_range==None:
        ax.plot(x/w,y/w, **kwargs)
    else:
        mask=np.where((t > t_range[0]) & (t < t_range[1]))
        xp = x[mask]
        yp = y[mask]
        wp = w[mask]
        ax.plot(xp/wp,yp/wp, **kwargs)
    ax.axis('equal')

def plot_curve_2d( skv, tkv, sdeg, tdeg, pts, num_pts = 100, **kwargs ):
    # skv is the knot vector in the s direction
    # tkv is the knot vector in the t direction
    # sdeg is the degree of the polynomial in the s direction
    # tdeg is the degree of the polynomial in the t direction
    # pts is the control points in the s-fastest then t direction
    xgrid,ygrid,wgrid = eval_curve_2d(skv, tkv, sdeg, tdeg, pts, num_pts)
    
    fig = plt.figure()
    ax = axes3d.Axes3D(fig)
    ax.plot_wireframe(xgrid/wgrid, ygrid/wgrid, 1, rstride=2, cstride=2, **kwargs)
    plt.axis('equal')        
    plt.show()            

def eval_curve_2d( skv, tkv, sdeg, tdeg, pts, num_pts = 100 ):
    sfuncs = spline_basis_funcs( sdeg, skv )
    tfuncs = spline_basis_funcs( tdeg, tkv )
    s = np.linspace(skv[0],skv[-1], num_pts)
    t = np.linspace(tkv[0],tkv[-1], num_pts)
    sknotlen = len(skv) # length of the knot vector in the s direction
    tknotlen = len(tkv) # length of the knot vector in the t direction
    scontpoint = sknotlen-sdeg-1 # number of control points in the s direction
    tcontpoint = tknotlen-tdeg-1 # number of control points in the t direction
    wgrid = np.zeros([num_pts,num_pts]) # grid of w values
    xgrid = np.zeros(wgrid.shape) # grid of x values
    ygrid = np.zeros(wgrid.shape) # grid of y values
    # xvec = np.zeros(num_pts**2) # values of x in a vector
    # yvec = np.zeros(xvec.shape) # values of y in a vector
    # wvec = np.zeros(xvec.shape) # values of w in a vector
    for scurve in range(1,tcontpoint+1): # iterate through the s basis functions (in the t direction)
        for tcurve in range(1,scontpoint+1):
            currentpt = (tcurve-1) + (scurve-1)*scontpoint # find the point is overlapped by both basis functions
            xcurrent = pts[currentpt][0] # x value of current point
            ycurrent = pts[currentpt][1] # y value of current point
            wcurrent = pts[currentpt][2] # w value of current point
            xgrid += np.outer( sfuncs[tcurve-1](s), tfuncs[scurve-1](t) * xcurrent * wcurrent )
            ygrid += np.outer( sfuncs[tcurve-1](s), tfuncs[scurve-1](t) * ycurrent * wcurrent )
            wgrid += np.outer( sfuncs[tcurve-1](s), tfuncs[scurve-1](t) * wcurrent )
    # for spoint in range(0,num_pts):
    #     for tpoint in range(0,num_pts):
    #         xvec[spoint+tpoint*num_pts] = xgrid[spoint,tpoint]
    #         yvec[spoint+tpoint*num_pts] = ygrid[spoint,tpoint]
    #         wvec[spoint+tpoint*num_pts] = wgrid[spoint,tpoint]
    # return xvec, yvec, wvec
    return xgrid, ygrid, wgrid

def refine_knot_vector( knot_vector ):
    new_knot_vector = []
    knot_mults=get_knot_multiplicities(knot_vector)
    knots=get_unique_knots(knot_vector)
    for i, k in enumerate(knots[1:]):
        ki = knots[i]
        for m in range(knot_mults[i]):
            new_knot_vector.append(ki)
        new_knot_vector.append((ki+k)/2.0)
    for m in range(knot_mults[-1]):
        new_knot_vector.append(knots[-1])
    return new_knot_vector

def coarsen_knot_vector( vector ):
    new_knot_vector = []
    knot_mults=get_knot_multiplicities(knot_vector)
    knots=get_unique_knots(knot_vector)
    for k, mults in zip(knots[::2], knot_mults[::2]):
        for m in range(mults):
            new_knot_vector.append(k)
    return new_knot_vector
    
def elevate_knot_vector( knot_vector ):
    new_knot_vector = []
    knot_mults=get_knot_multiplicities(knot_vector)
    knots=get_unique_knots(knot_vector)
    for i, k in enumerate(knots):
        for m in range(knot_mults[i]+1):
            new_knot_vector.append(k)
    return new_knot_vector
    
def reduce_knot_vector( knot_vector ):
    new_knot_vector = []
    knot_mults=get_knot_multiplicities(knot_vector)
    knots=get_unique_knots(knot_vector)
    for i, k in enumerate(knots):
        for m in range(max(1,knot_mults[i]-1)):
            new_knot_vector.append(k)
    return new_knot_vector
    
def roughen_knot_vector( knot_vector ):
    new_knot_vector = []
    i = 0
    unique_knots = get_unique_knots(knot_vector)
    knot_mults = get_knot_multiplicities(knot_vector)
    for i in range(len(unique_knots)):
        if knot_mults[i] == knot_mults[0]:
            for j in range(knot_mults[i]):
                new_knot_vector.append(unique_knots[i])
        else:
            for j in range(knot_mults[i]+1):
                new_knot_vector.append(unique_knots[i])
    return new_knot_vector
    
def smooth_knot_vector( knot_vector ):
    new_knot_vector = []
    unique_knots = get_unique_knots(knot_vector)
    knot_mults = get_knot_multiplicities(knot_vector)
    for i in range(len(unique_knots)):
        knot = unique_knots[i]
        if knot == unique_knots[0] or knot == unique_knots[-1]:
            for j in range(knot_mults[i]):
                new_knot_vector.append(knot)
        elif knot_mults[i] == 1:
            new_knot_vector.append(knot)
        else:
            for j in range(knot_mults[i]-1):
                new_knot_vector.append(knot)
    return new_knot_vector
    
def shift_interior_knots( knot_vector, shift = 0.2 ):
    new_knot_vector = []
    knot_mults=get_knot_multiplicities(knot_vector)
    knots=get_unique_knots(knot_vector)
    for i, k in enumerate(knots):
        for m in range(knot_mults[i]):
            if i == 0 or i == len(knot_mults)-1:
                new_knot_vector.append(k)
            else:
                new_knot_vector.append(k + shift)
    return new_knot_vector

def format_knot_label(k_val, mult = 1, val_sep = ",", group_begin = "(", group_end = ")", **kwargs):
    label = ''
    if mult == 1:
        group_begin = ""
        group_end = ""
    try:
        if k_val.denominator == 1:
            label = label + "$" + group_begin + val_sep.join([str(k_val.numerator) for m in range(mult)]) + group_end + "$"
        else:
            label = label + r"$" + group_begin  + val_sep.join([r"\frac{" + str(k_val.numerator) + "}{" + str(k_val.denominator) + "}" for m in range(mult)]) + group_end + "$"
    except AttributeError:
        numer_denom = k_val.as_numer_denom()
        if numer_denom[1] == 1:
            label = label + "$ + group_begin " + val_sep.join([str(numer_denom[0]) for m in range(mult)]) + group_end + "$"
        else:
            label = label + r"$" + group_begin  + val_sep.join([r"\frac{" + str(numer_denom[0]) + "}{" + str(numer_denom[1]) + "}" for m in range(mult)]) + group_end + "$"
    except AttributeError:
        label = label + r"$" + group_begin  + val_sep.join([str(k_val) for m in range(mult)]) + group_end + "$"
    return label

def get_knot_vector_labels( knot_vector, label_repeated = True ):

    unique_knots = get_unique_knots( knot_vector )
    try:
        unique_knots = [Fraction(knot) for knot in unique_knots]
    except TypeError:
        numer, denom = knot.as_numer_denom()
        unique_knots = [Fraction(int(numer), int(denom)) for knot in unique_knots]
    knot_mults = get_knot_multiplicities( knot_vector )
    if label_repeated:
        knot_strings = [format_knot_label(k,mult = m) for k,m in zip(unique_knots,knot_mults)]
    else:
        knot_strings = [format_knot_label(k) for k in unique_knots]

    print knot_strings

def plot_basis(knot_vector, degree, num_pts = 100, ax = None, color_cycle = None, **kwargs):

    if ax == None:
        ax = plt.gca()

    nkwargs = kwargs.copy()
    for local_kv in window(knot_vector, degree + 2):
        if color_cycle != None:
            nkwargs['color'] = color_cycle.next()
        func = spline_basis(local_kv)
        func_bounds = [local_kv[0],local_kv[-1]]
        t = np.linspace(func_bounds[0],func_bounds[1],num_pts)
        ax.plot(t, func(t), **nkwargs)

def plot_curve2( kv, degree, ws, num_pts = 100, yp = 0, ax = None, **kwargs ):

    funcs = spline_basis_funcs( degree, kv )
    t = np.linspace(kv[0],kv[-1], num_pts)
    y = np.zeros(t.shape)
    if ax == None:
        ax = plt.gca()
    for w, func in zip(ws, funcs):
        f_vals = func(t)
        # ax.plot(t,f_vals)
        y += w * f_vals
    y += yp
    ax.plot(t,y, **kwargs)

####################### TEST CODE ################################

def maintest():
    weights = get_elem_func_weights(2, [0,0,0,1,2,3,4,4,4])
    # b = BFunc([0,1,2,3,4],[0,1,2,3])
    

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # x = np.linspace(0,4)
    # y = np.linspace(0,3)
    # X, Y = np.meshgrid(x,y)
    # Z = b(X,Y)
    # ax.plot_wireframe(X, Y, Z, rstride=5, cstride=5)

    # fig1 = plt.figure()
    # ax = fig1.add_subplot(111)
    # ax.plot(x,b.internal_func(x, 0))
    # ax.plot(y,b.internal_func(y, 1))
    # plt.show()

    
    kv1 = [0,0,0,0.5,1,1,1]
    degree1 = 2
    pts1 = [[1,2],[3,5],[4,5],[5,2]]
    kv2 = [0,0,0,0.7,1,1,1]
    degree2 = 2
    local_pts2 = new_project_pts( kv1, degree1, kv2, degree2, pts1)
    pts2 = smooth_pts( kv2, degree2, local_pts2 )
    fig = plt.figure()
    fig.add_subplot(111)
    plot_curve(kv1,degree1,pts1, linestyle = '-', c = 'b')
    plot_curve(kv2,degree2,pts2, linestyle = '--', c = 'r')
    #fig.show()
    plt.show()

    degree1 = 2

    pts1 = copy.copy(pts2)
    kv1 = [0,0,0,0.7,1,1,1]
    kv2 = [0,0,0,0.5,1,1,1]
    # pts1[3][0] = 4
    # kv1 = [0,0,0,0.3,0.5,0.7,1,1,1]
    # degree2 = 2
    # kv2 = [0,0,0,1,1,1]
    local_pts2 = new_project_pts( kv1, degree1, kv2, degree2, pts1)
    pts2 = smooth_pts( kv2, degree2, local_pts2 )
    # fig = plt.figure()
    # fig.add_subplot(111)
    # plot_curve(kv1,degree1,pts1, linestyle = ':', c = 'g')
    plot_curve(kv2,degree2,pts2, linestyle = '--', c = 'orange')
    #fig.show()
    plt.show()

    # kv1 = [0,0,0,0.5,1,1,1]
    # degree1 = 2
    # pts1 = [[1,2],[3,6], [4,7],[5,6], [6,1]]
    # kv2 = [0,0,0,0.5,0.5, 1,1,1]
    # degree2 = 2
    # local_pts2 = new_project_pts( kv1, degree1, kv2, degree2, pts1)
    # pts2 = smooth_pts( kv2, degree2, local_pts2 )
    # fig = plt.figure()
    # fig.add_subplot(111)
    # plot_curve(kv1,degree1,pts1)
    # plot_curve(kv2,degree2,pts2)
    # fig.show()
    # # # fig1, fig2 = convergence_compare()
    # # # fig2.show()
    # # cpt_plot()

maintest()

if __name__ == "__main__":
    weights = get_elem_func_weights(2, [0,0,0,1,2,3,4,4,4])
    # b = BFunc([0,1,2,3,4],[0,1,2,3])
    

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # x = np.linspace(0,4)
    # y = np.linspace(0,3)
    # X, Y = np.meshgrid(x,y)
    # Z = b(X,Y)
    # ax.plot_wireframe(X, Y, Z, rstride=5, cstride=5)

    # fig1 = plt.figure()
    # ax = fig1.add_subplot(111)
    # ax.plot(x,b.internal_func(x, 0))
    # ax.plot(y,b.internal_func(y, 1))
    # plt.show()

    
    kv1 = [0,0,0,0.5,1,1,1]
    degree1 = 2
    pts1 = [[1,2],[3,5],[4,5],[5,2]]
    kv2 = [0,0,0,0.7,1,1,1]
    degree2 = 2
    local_pts2 = new_project_pts( kv1, degree1, kv2, degree2, pts1)
    pts2 = smooth_pts( kv2, degree2, local_pts2 )
    fig = plt.figure()
    fig.add_subplot(111)
    plot_curve(kv1,degree1,pts1, linestyle = '-', c = 'b')
    plot_curve(kv2,degree2,pts2, linestyle = '--', c = 'r')
    fig.show()

    degree1 = 2

    pts1 = copy.copy(pts2)
    kv1 = [0,0,0,0.7,1,1,1]
    kv2 = [0,0,0,0.5,1,1,1]
    # pts1[3][0] = 4
    # kv1 = [0,0,0,0.3,0.5,0.7,1,1,1]
    # degree2 = 2
    # kv2 = [0,0,0,1,1,1]
    local_pts2 = new_project_pts( kv1, degree1, kv2, degree2, pts1)
    pts2 = smooth_pts( kv2, degree2, local_pts2 )
    # fig = plt.figure()
    # fig.add_subplot(111)
    # plot_curve(kv1,degree1,pts1, linestyle = ':', c = 'g')
    plot_curve(kv2,degree2,pts2, linestyle = '--', c = 'orange')
    fig.show()
    
    # kv1 = [0,0,0,0.5,1,1,1]
    # degree1 = 2
    # pts1 = [[1,2],[3,6], [4,7],[5,6], [6,1]]
    # kv2 = [0,0,0,0.5,0.5, 1,1,1]
    # degree2 = 2
    # local_pts2 = new_project_pts( kv1, degree1, kv2, degree2, pts1)
    # pts2 = smooth_pts( kv2, degree2, local_pts2 )
    # fig = plt.figure()
    # fig.add_subplot(111)
    # plot_curve(kv1,degree1,pts1)
    # plot_curve(kv2,degree2,pts2)
    # fig.show()
    # # # fig1, fig2 = convergence_compare()
    # # # fig2.show()
    # # cpt_plot()
