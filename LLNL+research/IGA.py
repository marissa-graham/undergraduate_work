import numpy as np
import h5py
import copy
import sympy
import itertools
import operator
import string
import gauss_integrate
from fractions import Fraction
from mayavi import mlab

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
import local_projection as lp
import multigrid as mg

def notes():
	"""
	- Test code goes in a separate file
	- WATCH OUT FOR NAMESPACE ERRORS
	- The helper functions will go in their own files later, but for now they can stay here.
		- relax_code.py
		- residual_code.py
		- misc_IGA_helpers.py (most projection helpers should come from local_projection.py)


	Cases:
	- Simplest: -u'' + sig*u = f, one dimension, affine basis, 
	Start with the affine functions (take multigrid code from before), write it in terms
	of those, uniform grid spacing, -u'' + sig*u = f, easy as possible
	- Model problems
	- # of dimensions
	- Relaxation methods (dependent on model problem)
	- Projection methods
	- Bezier segment to actual splines
	- Interior knot multiplicity
		
	
	MAIN THING I AM STUCK ON:

	Normally, we have the grid, and the function is defined on the grid. With the splines, 
	we have the knot vector (which is analogous to the grid for fine/coarseness), then the 
	control points, which define the curve, and then the ODE/PDE we want to solve. Trying
	to figure out how to visualize this. I'm assuming it still turns into the form Av=f, 
	where A is a problem-dependent matrix with a structure that is taken advantage of.
	How does this work when the problem is defined on a curve? I need pictures.
	Is flat grid a special case of splines, and if so, how?

	Once I know this I will be able to better figure out what arguments each function needs.
	"""
	pass

####### HELPER CODE TO OFFLOAD LATER ################

# Discretization stuff based on a knot vector and basis type (calculate h vector)

# Going to need something for spline arc lengths. Might be in local_projection

################ MAIN FUNCTIONS #####################

# These will mostly call subroutines, like, 1d, 2d, flat vs spline, etc.

# Basis functions and control points
	# basis function is entirely dependent on the knot vector, so everything needs it
	# What depends on the discretization?
		# Relaxation and residual. Other than that, if it's in prob_params, it's fine

# What goes in prob_params?
	# name (of problem)
	# knot_vector
	# basis type (affine or spline)
	# projection type
	# problem-dependent parameters
	# dimension, I think

# What goes in relaxm?
	# relaxation method
	# number of sweeps
	# problem dependent parameters (e.g. w)

def relax(v, f, prob_params, relaxm, sweeps=1):
	pass

def residual(v, f, prob_params):
	pass

# Needs the knot vector, degree? control pts? prob_params?
def project(v, refine, method):
	# refine is a bool. true = refine, false = coarsen
	# method is a dictionary which includes the method and any parameters relevant to it

	# Flat grid case, as before. Call subfunctions for interpolation, restriction, etc.

	# Spline case
		# Create new knot vector (not changing degree)
		# Call new_project_pts, I think
	
	# When it turns into the spline case, the only extra argument is the knot vector
	pass

# TODO for V-cycle and FMG: figure out how initial guess will be adjusted
# For splines, the knot vector will be size 2**L + 2*degree

def V_cycle(vh, fh, L, l, v1, v2, prob_params, relaxm, projectm):
	"""
	Perform a V-cycle of height L using the given parameters.

	INPUTS
	------
	vh - Initial guess.
	fh - Right hand vector.
	L  - Finest level of grid spacing (excluding end knots).
	l  - Coarsest level of grid spacing (Note l <= L).
	v1 - Number of relaxation sweeps to perform on the way down.
	v2 - Number of relaxation sweeps to perform on the way up.
	prob_params - Relevant problem information to pass to helper functions.
	relaxm - Relaxation type to pass to helper functions.
	projectm - Projection method.

	OUTPUTS
	-------
	v - Updated initial guess.
	"""

	def V_h(v,f,k):
		"""Recursive V-cycling, starting on level k."""
		# Relax on prob, f v1 times with initial guess v
		# Grid size will be 2^k - 1
		v = relax(v, f, prob_params, relaxm, sweeps=v1)

		# Base case
		if k == l:
			# Relax v2 times on prob,f with initial guess v
			v = relax(v, f, prob_params, relaxm, sweeps=v2)
			return v

		else:
			r = residual(v, f, prob_params)
			f2h = project(r, False, projectm)

			# TODO: Get initial guess based on problem type
			if probl == 1:
				v2h = np.zeros(2**(k-1)-1)
			elif probl == 2:
				v2h = np.zeros((2**(k-1)-1,2**(k-1)-1))

			vh = V_h(v2h, f2h, k-1)
			
			v += project(vh, True, projectm)
			v = relax(v, f, prob_params, relaxm, sweeps=v2)
			
		return v

	return V_h(vh,fh,L)

# L and l still works with the splines, right? Right now we're assuming no interior knot
# multiplicity greater than one, and cells are doubled or halved.
def FMG(fh, L, l, v0, v1, v2, prob_params, relaxm, projectm):
	"""
	Perform a full multigrid V-cycle using the given parameters.

	INPUTS
	------
	fh - Right hand vector on the finest grid spacing.
	L  - Finest level of grid spacing (finest grid is size 2**L).
	l  - Coarsest level of grid spacing (Note l <= L).
	v0 - Number of times to run the V-cycle on each level.
	v1 - Number of relaxation sweeps to perform on the way down.
	v2 - Number of relaxation sweeps to perform on the way up.
	prob_params - Relevant problem information to pass to helper functions.
	relaxm - Relaxation type, to pass to helper functions.

	OUTPUTS
	-------
	v - Updated initial guess.
	"""
	f = {'f{0}h'.format(L):fh}
	for k in xrange(L-1,l-1,-1):
		# 'refine' argument is False because restriction is coarsening.
		f['f{0}h'.format(k)] = project(f['f{0}h'.format(k+1)], False, projectm)
		
	# TODO: adjust initial guess and where dimension is coming from
	if len(fh.shape) == 1:
		initial = np.zeros(2**l-1)
	elif len(fh.shape) == 2:
		initial = np.zeros((2**l-1,2**l-1))

	v = relax(initial, f['f{0}h'.format(l)], problem, relaxm, sweeps=v0)
	v = project(v, True, projectm)

	# Start running V-cycles; will start at l+1, go all the way to level L
	for k in xrange(l+1, L):
		for j in xrange(v0):
			v = V_cycle(v, f['f{0}h'.format(k)], k, l, v1, v2, prob_params, relaxm)
		v = project(v, True, projectm)

	return v
