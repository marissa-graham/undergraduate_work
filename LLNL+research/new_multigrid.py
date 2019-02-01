import numpy as np 
from scipy import linalg as la
from matplotlib import pyplot as plt
#from mayavi import mlab

# Extensions: nontrivial boundary conditions
#			  1d model problems besides poisson

########## HELPER FUNCTIONS ###############

# TODO: sig =/= 0 case
def jacobi_relax1d(v, f, prob_params, relaxm):
	k = v.shape[0]
	n = k+1

	if prob_params['prob'] == 'poisson1d':

		for j in xrange(prob_params['sweeps']):
			w = relaxm['w']

			# v_j = (1-w)*v_j + w*(f + v_j-1*n^2 + v_j+1*n^2)/(2n^2 + sig)

			# If sig == 0, we can avoid introducing some instability?
			if prob_params['sig'] == 0:
				v += -w*v + w/2.*(f/n**2 + np.append(0.0,v[:-1]) + np.append(v[1:], 0.0))

			else:
				# Check and see later if it gets screwy
				v += n**2/(2.0*n**2 + sig)*(f/n**2 + np.append(0.0,v[:-1]) + np.append(v[1:], 0.0))

		return v

def jacobi_relax2d(v, f, prob_params, relaxm):
	pass

# TODO: sig =/= 0 case
def GS_relax1d(v, f, prob_params, relaxm):
	k = v.shape[0]
	n = k+1

	if prob_params['prob'] == 'poisson1d':

		for j in xrange(prob_params['sweeps']):
			if prob_params['sig'] == 0:
				# TODO: Should add 0.5*f/n**2 outside the looping? Faster, right?

				v[0] = 0.5*(v[1] + f[0]/n**2)

				for i in xrange(1, n-2):
					v[i] = 0.5*(v[i-1] + v[i+1] + f[i]/n**2)

				v[n-2] = 0.5*(v[n-3] + f[n-2]/n**2)

		return v

def GS_relax12d(v, f, prob_params, relaxm):
	pass

def poisson1d_residual():
	pass

########## MAIN FUNCTIONS #################

def relax(v, f, prob_params, relaxm):
	"""
	Relax on a given model problem using a given method.

	INPUTS:
	v - Initial guess. Array of same dimension as problem.
	f - Right hand vector. Same dimension as problem.
	TODO: if 2d, guaranteed square?
	prob_params - Dictionary with problem parameters. 
		prob - model problem (poisson1d, poisson2d, etc.)
		dim - int, dimension of problem
		sweeps - int, number of relaxations to perform
		sig - float, parameter for Poisson problem
			If sig is None, just put sig=0
		eps - float, optional convection term for 2d Poisson
	relaxm - Relaxation method to use
		method - jacobi, GS, RBGS, SOR, etc.
		w - Weighting factor for Jacobi and SOR

	RETURNS:
	v - Updated solution vector.
	"""

	if relaxm['method'] == 'jacobi':
		if prob_params['dim'] == 1:
			v = jacobi_relax1d(v, f, prob_params, relaxm)
		elif prob_params['dim'] == 2:
			v = jacobi_relax2d(v, f, prob_params, relaxm)

	if relaxm['method'] == 'GS':
		if prob_params['dim'] == 1:
			v = GS_relax1d(v, f, prob_params, relaxm)
		elif prob_params['dim'] == 2:
			v = GS_relax2d(v, f, prob_params, relaxm)

	return v
			
def residual(v, f, prob_params):
	"""
	Find f-Av, where A is the matrix corresponding to the
	model problem given in prob_params.
	"""
	pass

def interpolate():
	pass

def restrict():
	pass

def V_cycle():
	pass

def FMG():
	pass

########## TESTING #######################

def test_relax(cases):
	# Cases is a vector of bools that decides which
	# test cases to run when test_jacobi is called

	# Plot Jacobi iteration on Au=0 for various wavenumbers, sig=0, w=2/3.
	if cases[0]:
		# Produces plot that looks like Fig 2.3(a)
		# Error terms look the same as from old multigrid.py
		n = 63
		f = np.zeros(n)
		domain = np.arange(100)

		prob_params = {'prob':'poisson1d', 'dim':1, 'sweeps':1, 'sig':0.0}
		relaxm = {'method':'jacobi', 'w':2/3.}

		# Reproduce plot from multigrid book
		plt.subplot(121)
		for k in [1,3,6]:
			v = np.sin(k*np.pi*np.linspace(0,1,n))
			error = np.zeros(100)
			error[0] = la.norm(-1*v,ord=np.inf)
			for i in xrange(1,100):
				v = relax(v, f, prob_params, relaxm)
				error[i] = la.norm(-1*v,np.inf)
			print error[-1]
			plt.plot(domain, error)
			
		# Whole bunch of wavenumbers

		plt.subplot(122)
		for k in xrange(1,150,5):
			v = np.sin(k*np.pi*np.linspace(0,1,n))
			error = np.zeros(100)
			error[0] = la.norm(-1*v,ord=np.inf)
			for i in xrange(1,100):
				v = relax(v, f, prob_params, relaxm)
				error[i] = la.norm(-1*v,np.inf)
			plt.plot(domain, error)
			
		plt.show()

	# Plot Gauss-Siedel iteration on Au=0 for various wavenumbers, sig=0
	if cases[1]:
		# Produces plot that looks like Fig 2.3(b)
		# Convergence is significantly better than Jacobi; see error[-1]
			# k = 1: 0.91875 vs 0.77819
			# k = 3: 0.48561 vs 0.14436
			# k = 6: 0.07938 vs 0.02097
		# Error terms still the same as in previous implementation

		n = 63
		f = np.zeros(n)
		domain = np.arange(100)

		prob_params = {'prob':'poisson1d', 'dim':1, 'sweeps':1, 'sig':0.0}
		relaxm = {'method':'GS', 'w':2/3.}

		# Reproduce plot from multigrid book
		plt.subplot(121)
		for k in [1,3,6]:
			v = np.sin(k*np.pi*np.linspace(0,1,n))
			error = np.zeros(100)
			error[0] = la.norm(-1*v,ord=np.inf)
			for i in xrange(1,100):
				v = relax(v, f, prob_params, relaxm)
				error[i] = la.norm(-1*v,np.inf)
			print error[-1]
			plt.plot(domain, error)
			
		# Whole bunch of wavenumbers

		plt.subplot(122)
		for k in xrange(1,150,5):
			v = np.sin(k*np.pi*np.linspace(0,1,n))
			error = np.zeros(100)
			error[0] = la.norm(-1*v,ord=np.inf)
			for i in xrange(1,100):
				v = relax(v, f, prob_params, relaxm)
				error[i] = la.norm(-1*v,np.inf)
			plt.plot(domain, error)
			
		plt.show()

	# Homogenous Jacobi iteration
	if cases[2]:
		pass

test_relax([1,1])

