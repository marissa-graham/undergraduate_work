import numpy as np
from scipy import linalg as la
from matplotlib import pyplot as plt
#from mayavi import mlab


def relax(v, f, prob, method, sweeps=1, sig=0., e=1, w=2/3.):
	"""Relax on the 1d or 2d model problem using a given method.
	1-d problem: -u'' + sig*u = f
	2-d problem: -u_xx - u_yy + sig*u = f
	If e != 1: -e(u_xx + u_yy) + sig*u_x = f

	INPUTS:
	v - Initial guess. If prob is 2-d, will be a 2-d array.
	f - Right hand vector. If prob is 2-d, will be a 2-d array.
	prob - Model problem to use. 1 is the 1-d problem, 2 is the
		standard 2-d problem (convection optional).
	sweeps - Number of relaxations to perform.
	method - Relaxation method to use; can be 'jacobi', 'GS',
		'RBGS', or 'SOR'.
	sig - Parameter in model problem. A float.
	e - Parameter in model problem 2c
	w - Weighting factor for Jacobi method and SOR.

	RETURNS:
	v - The updated solution vector
	"""
	k = v.shape[0]
	n = k+1

	# 1-d case. v is a vector
	if prob == 1:
		c = 2.0+sig/n**2

		if method == 'jacobi':
			for j in xrange(sweeps):
				vprev = np.append(0.0,v[:-1])
				vplus = np.append(v[1:],0.0)
				v += -w*v + w*(vprev + vplus + f/n**2)/c

		elif method == 'GS':
			for j in xrange(sweeps):
				v[0] = (v[1] + f[0]/n**2)/c
				for i in xrange(1,n-2):
					v[i] = (v[i-1] + v[i+1] + f[i]/n**2)/c
				v[n-2] = (v[n-3] + f[n-2]/n**2)/c
				#v /= c

		elif method == 'RBGS':
			for j in xrange(sweeps):
				v[0] = (v[1] + f[0]/n**2)/c
				# Evens; start at 1 since v_2 is v[1]
				for i in xrange(1,n-2,2):
					v[i] = (v[i-1] + v[i+1] + f[i]/n**2)/c
				# Odds
				for i in xrange(2,n-2,2):
					v[i] = (v[i-1] + v[i+1] + f[i]/n**2)/c
				v[n-2] = (v[n-3] + f[n-2]/n**2)/c
				#v /= c

		else:
			raise ValueError('Relaxation method not recognized.')

		return v

	# 2-d case. v is a 2-d array. Is it square? Yes, for now
	elif prob == 2:
		const = 4.0*e + 1.0*sig/n**2
		print const
		if method == 'jacobi':
			for j in xrange(sweeps):
				viprev = np.hstack((np.zeros((k,1)), v[:,:-1]))
				viplus = np.hstack((v[:,1:], np.zeros((k,1))))
				vjprev = np.vstack((np.zeros((1,k)), v[:-1,:]))
				vjplus = np.vstack((v[1:,:], np.zeros((1,k))))
				eterms = viprev + viplus + vjprev + vjplus
				v += -w*v 
				print f, eterms
				v += w*(f/n**2 + eterms)/const

		elif method == 'GS':
			for a in xrange(sweeps):
				for i in xrange(k):
					for j in xrange(k):
						if i == 0:
							if j == 0:
								eterms = e*(v[i,j+1] + v[i+1,j])
							elif j == k-1:
								eterms = e*(v[i,j-1] + v[i+1,j])
							else:
								eterms = e*(v[i,j-1] + v[i,j+1] + v[i-1,j] + v[i+1,j])
						elif i == k-1:
							if j == 0:
								eterms = e*(v[i,j+1] + v[i-1,j])
							elif j == k-1:
								eterms = e*(v[i,j-1] + v[i-1,j])
							else:
								eterms = e*(v[i,j-1] + v[i,j+1] + v[i-1,j])
						else:
							if j == 0:
								eterms = e*(v[i,j+1] + v[i-1,j] + v[i+1,j])
							elif j == k-1:
								eterms = e*(v[i,j-1] + v[i-1,j] + v[i+1,j])
							else:
								eterms = e*(v[i,j-1] + v[i,j+1] + v[i-1,j] + v[i+1,j])
					
						v[i,j] = (f[i,j]/n**2 + eterms)/const

		elif method == 'RBGS':
			for z in xrange(sweeps):
				# Red points, i.e. even indices
				for i in xrange(k):
					for j in xrange(i % 2, k, 2):
						if i == 0:
							if j == 0:
								eterms = e*(v[i,j+1] + v[i+1,j])
							elif j == k-1:
								eterms = e*(v[i,j-1] + v[i+1,j])
							else:
								eterms = e*(v[i,j-1] + v[i,j+1] + v[i-1,j] + v[i+1,j])
						elif i == k-1:
							if j == 0:
								eterms = e*(v[i,j+1] + v[i-1,j])
							elif j == k-1:
								eterms = e*(v[i,j-1] + v[i-1,j])
							else:
								eterms = e*(v[i,j-1] + v[i,j+1] + v[i-1,j])
						else:
							if j == 0:
								eterms = e*(v[i,j+1] + v[i-1,j] + v[i+1,j])
							elif j == k-1:
								eterms = e*(v[i,j-1] + v[i-1,j] + v[i+1,j])
							else:
								eterms = e*(v[i,j-1] + v[i,j+1] + v[i-1,j] + v[i+1,j])
					
						v[i,j] = (f[i,j]/n**2 + eterms)/const

				# Black points, i.e. odd indices
				for i in xrange(k):
					for j in xrange((i+1) % 2, k, 2):
						if i == 0:
							if j == 0:
								eterms = e*(v[i,j+1] + v[i+1,j])
							elif j == k-1:
								eterms = e*(v[i,j-1] + v[i,j+1] + v[i-1,j] + v[i+1,j])
							else:
								eterms = e*(v[i,j-1] + v[i,j+1] + v[i-1,j] + v[i+1,j])
						elif i == k-1:
							if j == 0:
								eterms = e*(v[i,j+1] + v[i-1,j])
							elif j == k-1:
								eterms = e*(v[i,j-1] + v[i-1,j])
							else:
								eterms = e*(v[i,j-1] + v[i,j+1] + v[i-1,j])
						else:
							if j == 0:
								eterms = e*(v[i,j+1] + v[i-1,j] + v[i+1,j])
							elif j == k-1:
								eterms = e*(v[i,j-1] + v[i-1,j] + v[i+1,j])
							else:
								eterms = e*(v[i,j-1] + v[i,j+1] + v[i-1,j] + v[i+1,j])
						
						v[i,j] = (f[i,j]/n**2 + eterms)/const

		else:
			raise ValueError('Relaxation method not recognized.')

		return v

def residual(v, f, prob, sig=0, e=1):
	"""Find f-Av, where A is the matrix corresponding to prob.
	Prob can be 1, 2, '2c'.
	"""
	k = v.shape[0]
	n = k + 1
	if prob == 1:
		u2prime = (-np.append(0,v[:-1]) + 2*v - np.append(v[1:],0))*n**2
		return f - u2prime - sig*v

	if prob == 2:
		viprev = np.hstack((np.zeros((k,1)), v[:,:-1]))
		viplus = np.hstack((v[:,1:], np.zeros((k,1))))
		vjprev = np.vstack((np.zeros((1,k)), v[:-1,:]))
		vjplus = np.vstack((v[1:,:], np.zeros((1,k))))
		eterms = viprev + viplus + vjprev + vjplus
		return f + eterms*n**2 - 4*v*n**2 - sig*v

def interpolate(v2h, method='linear'):
	"""Take a vector from a 2h grid spacing to h grid spacing; that
	is, we have a linear operator from R^{n/2-1} to R^{n-1}.

	INPUTS:
	v2h - Input vector. If d==2, will be a 2-d array (since we 
		don't actually use matrix multiplication to relax).
	d - Number of dimensions (1 or 2).
	method - Whether to use linear interpolation or something else.

	RETURNS:
	vh - Transformed vector.
	"""
	d = len(v2h.shape)
	if method == 'linear':
		if d == 1:

			k = v2h.shape[0]
			vh = np.zeros(2*k+1)

			# v_j and v_j+1
			vh[::2] += np.append(0,v2h)
			vh[::2] += np.append(v2h,0)
			vh[1::2] += 2*v2h
			vh /= 2

		if d == 2:
			k = v2h.shape[0]
			khzeros = np.zeros((1,k))
			kvzeros = khzeros.T
			vh = np.zeros((2*k+1, 2*k+1))

			#print 'Interpolating...'
			#X, Y = np.meshgrid(np.linspace(0,1,k), np.linspace(0,1,k))
			#mlab.surf(v2h, warp_scale='auto')
			#mlab.show()

			# v_2i,2j case
			vh[1::2,1::2] += v2h

			# v_2i+1,2j case
			vh[::2,1::2] += 0.5*(np.vstack((khzeros,v2h)) + np.vstack((v2h,khzeros)))
			
			# v_2i,2j+1 case
			vh[1::2,::2] += 0.5*(np.hstack((v2h,kvzeros)) + np.hstack((kvzeros,v2h)))
			
			# v_2i+1,2j+1 case
			b = np.zeros((k+1, k+1))
			b1, b2, b3, b4 = b, b, b, b

			b1[1:,1:], b2[:-1,1:], b3[1:,:-1], b4[:-1,:-1] = v2h, v2h, v2h, v2h
			vh[::2,::2] += 0.25*(b1 + b2 + b3 + b4)

			#dom = np.linspace(0,1,2*k+1)
			#X, Y = np.meshgrid(dom, dom)
			#mlab.surf(vh, warp_scale='auto')
			#mlab.show()

		return vh

	else:
		raise ValueError('Interpolation method not recognized.')

def restrict(vh, method):
	"""Take a vector from h grid spacing to 2h grid spacing; that
	is, we have a linear operatore from R^{n-1} to R^{n/2-1}.

	INPUTS:
	vh - Input vector. If d==2, will be a 2-d array (since we
		don't actually use matrix multiplication to relax)
	d - Number of dimensions (1 or 2).
	method - What sort of restriction to use. Options are 'full' for
		full weighting, 'injection' for injection, and 'half' for
		half-injection.

	RETURNS:
	v2h - Transformed vector.
	"""
	d = len(vh.shape)

	if d == 1:
		if method == 'full':
			v2h = 0.5*vh[1::2] + 0.25*(vh[:-1:2]+vh[2::2])

		elif method == 'injection':
			v2h = vh[1::2]

		elif method == 'half':
			v2h = 0.5*vh[1::2]

		else:
			raise ValueError('Injection method not recognized.')

		return v2h

	elif d == 2:
		k = vh.shape[0]

		#print 'Restricting...'
		#X, Y = np.meshgrid(np.linspace(0,1,k), np.linspace(0,1,k))
		#mlab.surf(vh, warp_scale='auto')
		#mlab.show()
		#print vh
	
		if method == 'full':
			singles = vh[:-1:2,:-1:2] + vh[:-1:2,2::2] + vh[2::2,:-1:2] + vh[2::2,2::2]
			doubles = vh[:-1:2,1::2] + vh[1::2,:-1:2] + vh[1::2,2::2] + vh[2::2,1::2]
			v2h = (singles + 2*doubles + 4*vh[1::2,1::2])/16.0

		elif method == 'injection':
			v2h = vh[1::2,1::2]

		elif method == 'half':
			v2h = 0.5*vh[1::2,1::2]

		else:
			raise ValueError('Injection method not recognized.')

		#dom = np.linspace(0,1,v2h.shape[0])
		#X, Y = np.meshgrid(dom, dom)
		#mlab.surf(v2h, warp_scale='auto')
		#mlab.show()

		#return v2h

	else:
		raise ValueError('Dimension not recognized.')

def V_cycle(vh, fh, L, l, probl, relaxm, restrictm, v1, v2, sigma=0, eps=1, om=2/3.):
	"""
	Perform a V-cycle, using the parameters of the model problem
	as given, with v1 sweeps on the way down, v2 sweeps on the way
	up, starting on a grid of size 2^L, going down to a grid of size
	2^l as the coarsest level.

	INPUTS:
	vh - Initial guess.
	fh - Right hand vector.
	L - Finest level of the grid; we start out with size 2^L.
	l - Coarsest level of the grid; we stop at size 2^l.
	prob - Model problem whose solution we are computing. 1 is the 
		standard 1-d problem, 2 is the standard 2-d problem, and
		2c is the 2-d problem with a convection term.
	relaxm - Relaxation method to use; can be 'jacobi', 'GS', or 'RBGS'.
	restrictm - Restriction method to use. Can be 'full', 'injection',
		or 'half'.
	v1 - Number of relaxation sweeps to perform on the way down.
	v2 - Number of relaxation sweeps to perform on the way up.
	sig - Parameter in model problem
	e - Parameter in model problem 2c
	w - Weighting factor for Jacobi method

	OUTPUT:
	v - Updated initial guess.
	"""

	def V_h(v,f,k):
		"""Recursive V-cycling, starting on level k."""
		# Relax on prob, f v1 times with initial guess v
		# Grid size will be 2^k - 1
		v = relax(v, f, probl, relaxm, sweeps=v1, sig=sigma, e=eps, w=om)

		# Base case
		if k == l:
			# TODO: shouldn't the base case be solving?
			
			# Relax v2 times on prob,f with initial guess v
			v = relax(v, f, probl, relaxm, sweeps=v2, sig=sigma, e=eps, w=om)
			return v

		else:
			r = residual(v, f, prob=probl, sig=sigma, e=eps)
			f2h = restrict(r, restrictm)
			if probl == 1:
				v2h = np.zeros(2**(k-1)-1)
			elif probl == 2:
				v2h = np.zeros((2**(k-1)-1,2**(k-1)-1))
			vh = V_h(v2h, f2h, k-1)
			
			v += interpolate(vh)
			v = relax(v, f, probl, relaxm, sweeps=v2, sig=sigma, e=eps, w=om)
			
		return v

	return V_h(vh,fh,L)

def FMG(fh, L, l, problem, relaxm, restrictm, v1, v2, v0=1, sig=0, e=1, w=2/3.):
	"""Perform a full multigrid V-cycle, using the parameters of the 
	model problem as given, with v1 sweeps on the way down, v2 sweeps 
	on the way up, initially solving on a grid of size 2^l (the coarsest level), and working
	to a grid of size 2^L. v0 V-cycles on each level.

	INPUTS: 
	fh - Right hand vector. Given in terms of the finest grid spacing.
	L - Finest level of the grid; we stop at size 2^L.
	l - Coarsest level of the grid; we start out at size 2^l. Note
		l must be <= L
	prob - Model problem whose solution we are computing. 1 is the 
		standard 1-d problem, 2 is the standard 2-d problem, and
		2c is the 2-d problem with a convection term.
	relaxm - Relaxation method to use. 'jacobi', 'GS', or 'RBGS'.
	restrictm - Restriction method. 'full', 'injection', or 'half'.
	v0 - Number of times to run the V-cycle on each level. 
		Defaults to 1.
	v1 - Number of relaxation sweeps to perform on the way down.
	v2 - Number of relaxation sweeps to perform on the way up.
	sig - Parameter in model problem
	e - Parameter in model problem 2c
	w - Weighting factor for Jacobi method

	OUTPUT:
	v - Updated initial guess.
	"""
	f = {'f{0}h'.format(L):fh}
	for k in xrange(L-1,l-1,-1):
		f['f{0}h'.format(k)] = restrict(f['f{0}h'.format(k+1)], restrictm)
	# Solve on coarsest grid
	if len(fh.shape) == 1:
		initial = np.zeros(2**l-1)
	elif len(fh.shape) == 2:
		initial = np.zeros((2**l-1,2**l-1))
	v = relax(initial, f['f{0}h'.format(l)], problem, relaxm, sweeps=v0)
	v = interpolate(v)

	# Start running V-cycles; will start at l+1, go all the way to level L
	for k in xrange(l+1, L):
		for j in xrange(v0):
			v = V_cycle(v,f['f{0}h'.format(k)], k, l, problem, relaxm, restrictm, v1, v2)
		v = interpolate(v)

	return v


########## TESTING ##########

# Tests relaxation in 1-d ONLY
def ch2ex20():
	n = 63
	f = np.zeros(n)
	domain = np.arange(100)

	# Weighted Jacobi iteration error plot
	plt.subplot(121)
	for k in [1,3,6]:
		v = np.sin(k*np.pi*np.linspace(0,1,n))
		error = np.zeros(100)
		error[0] = la.norm(-1*v,ord=np.inf)
		for i in xrange(1,100):
			v = relax(v, f, 1, method='jacobi')
			error[i] = la.norm(-1*v,np.inf)
		plt.plot(domain, error)
		plt.tick_params(axis='x',which='both',labelbottom='off',labelleft='off',bottom='off',left='off')
		plt.tick_params(axis='y',which='both',labelbottom='off',labelleft='off',bottom='off',left='off')
	

	plt.subplot(122)
	for k in xrange(1,150,5):
		v = np.sin(k*np.pi*np.linspace(0,1,n))
		error = np.zeros(100)
		error[0] = la.norm(-1*v,ord=np.inf)
		for i in xrange(1,100):
			v = relax(v, f, 1, method='jacobi')
			error[i] = la.norm(-1*v,np.inf)
		plt.plot(domain, error)
		print error[-1]
		plt.tick_params(axis='x',which='both',labelbottom='off',labelleft='off',bottom='off',left='off')
		plt.tick_params(axis='y',which='both',labelbottom='off',labelleft='off',bottom='off',left='off')
	

	
	# Gauss-Siedel iteration error plot
	plt.subplot(323)
	for k in [1,3,6]:
		v = np.sin(k*np.pi*np.linspace(0,1,n))
		error = np.zeros(100)
		error[0] = la.norm(-1*v,ord=np.inf)
		for i in xrange(1,100):
			v = relax(v, f, 1, method='GS')
			error[i] = la.norm(-1*v,np.inf)
		print error[-1]
		plt.plot(domain, error)

	plt.subplot(324)
	for k in xrange(1,150,5):
		v = np.sin(k*np.pi*np.linspace(0,1,n))
		error = np.zeros(100)
		error[0] = la.norm(-1*v,ord=np.inf)
		for i in xrange(1,100):
			v = relax(v, f, 1, method='GS')
			error[i] = la.norm(-1*v,np.inf)
		plt.plot(domain, error)


	"""
	# Red-black Gauss-Siedel
	plt.subplot(325)
	for k in [1,3,6]:
		v = np.sin(k*np.pi*np.linspace(0,1,n))
		error = np.zeros(100)
		error[0] = la.norm(-1*v,ord=np.inf)
		for i in xrange(1,100):
			v = relax(v, f, 1, method='RBGS')
			error[i] = la.norm(-1*v,np.inf)
		plt.plot(domain, error)
	
	plt.subplot(326)
	for k in xrange(1,150,5):
		v = np.sin(k*np.pi*np.linspace(0,1,n))
		error = np.zeros(100)
		error[0] = la.norm(-1*v,ord=np.inf)
		for i in xrange(1,100):
			v = relax(v, f, 1, method='RBGS')
			error[i] = la.norm(-1*v,np.inf)
		plt.plot(domain, error)
	""" 	
	plt.show()
	#plt.savefig("510projectJacobi.png",dpi=300)

ch2ex20()

def test_relax():
	m = 100
	k = 4
	L = 8
	C = 20.0
	n = 2**L
	pi = np.pi
	relaxm = 'jacobi'

	dom = np.linspace(0,1,n+1)[1:-1]
	v = (np.sin(pi*dom) + np.sin(6*pi*dom) + np.sin(32*pi*dom))/3.0
	f = C*np.sin(k*np.pi*dom)#np.linspace(0,1,n+1))[1:-1]
	u = C/(np.pi*np.pi*k*k)*np.sin(k*np.pi*dom)#np.linspace(0,1,n+1))[1:-1]
	plt.plot(dom, u, dom, v)

	error = np.zeros(m+2)
	error[0] = la.norm(-1*v,ord=np.inf)
	for i in xrange(m):
		v = relax(v, f, 1, method=relaxm)
		error[i] = la.norm(u - v,np.inf)
	plt.plot(dom, v)
	plt.plot(dom,u)
	plt.show()
	plt.plot(np.arange(m+2), error)
	plt.show()

	v = (np.sin(pi*dom) + np.sin(6*pi*dom) + np.sin(32*pi*dom))/3.0
	error = np.zeros(m+2)
	error[0] = la.norm(-1*v,ord=np.inf)
	for i in xrange(m):
		v = relax(v, f, 1, method='GS')
		error[i] = la.norm(u - v,np.inf)
	plt.plot(dom, u, dom, v)
	plt.show()
	plt.plot(np.arange(m+2), error)
	plt.show()

	v = (np.sin(pi*dom) + np.sin(6*pi*dom) + np.sin(32*pi*dom))/3.0
	error = np.zeros(m+2)
	error[0] = la.norm(-1*v,ord=np.inf)
	for i in xrange(m):
		v = relax(v, f, 1, method='RBGS')
		error[i] = la.norm(u - v,np.inf)
	plt.plot(dom, u, dom, v)
	plt.show()
	plt.plot(np.arange(m+2), error)
	plt.show()

# Test V-cycles for all combinations of methods, in 1-d
# Let's leave that alone for now
def testV():
	C = 2.0
	k = 4
	m = 10
	errors = np.zeros(9)
	for L in xrange(3,12):
		n = 2**L
		print '\n', 'L=',L, 'so n=',n, '\n'
		pi = np.pi
		f = C*np.sin(k*np.pi*np.linspace(0,1,n+1))[1:-1]
		#f = np.zeros(n-1)
		d = np.linspace(0,1,n+1)[1:-1]
		v = (np.sin(pi*d)+np.sin(6*pi*d)+np.sin(32*pi*d))/3.0
		u = C/(np.pi**2*k**2)*np.sin(k*np.pi*np.linspace(0,1,n+1))[1:-1]

		error = np.zeros(m)
		for i in xrange(m):
			v = V_cycle(v, f, L, 2, 1, 'jacobi', 'full', 2, 1)
			error[i] = la.norm(u-v)/2**(L/2.0)
		errors[L-3] = error[-1]
		print 'jacobi, full weighting:\t', error[-1]

		"""
		error = np.zeros(m)
		v = (np.sin(pi*d)+np.sin(6*pi*d)+np.sin(32*pi*d))/3.0
		for i in xrange(m):
			v = V_cycle(v, f, L, 2, 1, 'jacobi', 'injection', 2, 1)
			error[i] = la.norm(u-v)/2**(L/2.0)
		print 'jacobi, injection\t', error[-1]

		error = np.zeros(m)
		v = (np.sin(pi*d)+np.sin(6*pi*d)+np.sin(32*pi*d))/3.0
		for i in xrange(m):
			v = V_cycle(v, f, L, 2, 1, 'jacobi', 'half', 2, 1)
			error[i] = la.norm(u-v)/2**(L/2.0)
		print 'jacobi, half injection:\t', error[-1]
		
		error = np.zeros(m)
		v = (np.sin(pi*d)+np.sin(6*pi*d)+np.sin(32*pi*d))/3.0
		for i in xrange(m):
			v = V_cycle(v, f, L, 2, 1, 'GS', 'full', 2, 1)
			error[i] = la.norm(u-v)/2**(L/2.0)
		print 'GS, full weighting:\t', error[-1]
		
		error = np.zeros(m)
		v = (np.sin(pi*d)+np.sin(6*pi*d)+np.sin(32*pi*d))/3.0
		for i in xrange(m):
			v = V_cycle(v, f, L, 2, 1, 'GS', 'injection', 2, 1)
			error[i] = la.norm(u-v)/2**(L/2.0)
		print 'GS, injection:\t\t', error[-1]

		error = np.zeros(m)
		v = (np.sin(pi*d)+np.sin(6*pi*d)+np.sin(32*pi*d))/3.0
		for i in xrange(m):
			v = V_cycle(v, f, L, 2, 1, 'GS', 'half', 2, 1)
			error[i] = la.norm(u-v)/2**(L/2.0)
		print 'GS, half injection:\t', error[-1]
		
		error = np.zeros(m)
		v = (np.sin(pi*d)+np.sin(6*pi*d)+np.sin(32*pi*d))/3.0
		for i in xrange(m):
			v = V_cycle(v, f, L, 2, 1, 'RBGS', 'full', 2, 1)
			error[i] = la.norm(u-v)/2**(L/2.0)
		print 'RGBS, full weighting:\t', error[-1]

		error = np.zeros(m)
		v = (np.sin(pi*d)+np.sin(6*pi*d)+np.sin(32*pi*d))/3.0
		for i in xrange(m):
			v = V_cycle(v, f, L, 2, 1, 'RBGS', 'injection', 2, 1)
			error[i] = la.norm(u-v)/2**(L/2.0)
		print 'RGBS, injection:\t', error[-1]
		
		error = np.zeros(m)
		v = (np.sin(pi*d)+np.sin(6*pi*d)+np.sin(32*pi*d))/3.0
		for i in xrange(m):
			v = V_cycle(v, f, L, 2, 1, 'RBGS', 'half', 2, 1)
			error[i] = la.norm(u-v)/2**(L/2.0)
		print 'RGBS, half injection:\t', error[-1]
		"""
	plt.plot(np.arange(3,12),errors)
	plt.savefig('510projectVcycle.png',dpi=300)
	plt.show()

def testFMG():
	for L in xrange(3,12):
		n = 2**L
		print '\nn = ', n, '\n'
		C = 2.0
		k = 4
		f = C*np.sin(k*np.pi*np.linspace(0,1,n+1))[1:-1]
		u = C/(np.pi**2*k**2)*np.sin(k*np.pi*np.linspace(0,1,n+1))[1:-1]

		v = FMG(f, L, 2, 1, 'jacobi', 'full', 2, 1)
		print 'jacobi, full:\t\t', la.norm(u-v)/2**(L/2.0)

		v = FMG(f, L, 2, 1, 'jacobi', 'injection', 2, 1)
		print 'jacobi, injection:\t', la.norm(u-v)/2**(L/2.0)

		v = FMG(f, L, 2, 1, 'jacobi', 'half', 2, 1)
		print 'jacobi, half:\t\t',la.norm(u-v)/2**(L/2.0)

		v = FMG(f, L, 2, 1, 'GS', 'full', 2, 1)
		print 'GS, full:\t\t', la.norm(u-v)/2**(L/2.0)

		v = FMG(f, L, 2, 1, 'GS', 'injection', 2, 1)
		print 'GS, injection:\t\t', la.norm(u-v)/2**(L/2.0)

		v = FMG(f, L, 2, 1, 'GS', 'half', 2, 1)
		print 'GS, half:\t\t', la.norm(u-v)/2**(L/2.0)

		v = FMG(f, L, 2, 1, 'RBGS', 'full', 2, 1)
		print 'RBGS, full:\t\t', la.norm(u-v)/2**(L/2.0)

		v = FMG(f, L, 2, 1, 'RBGS', 'injection', 2, 1)
		print 'RBGS, injection:\t', la.norm(u-v)/2**(L/2.0)

		v = FMG(f, L, 2, 1, 'RBGS', 'half', 2, 1)
		print 'RBGS, half:\t\t', la.norm(u-v)/2**(L/2.0)
def test2():
	C = 2.0
	k = 4
	m = 10
	sigma = 0.1
	domain = np.arange(m)
	errors2 = np.zeros(m)
	for L in xrange(3,2+m):
		n = 2**L
		f = C*np.sin(k*np.pi*np.linspace(0,1,n+1))[1:-1]
		v = FMG(f, L, 2, 1, 'jacobi', 'full', 2, 1, sig=sigma)
		u = C/(np.pi**2*k**2+ sigma)*np.sin(k*np.pi*np.linspace(0,1,n+1))[1:-1]
		errors2[L-3] = la.norm(u-v)/np.sqrt(n)
	print errors
	plt.plot(domain, errors2)
	plt.savefig('510projectFMG.png',dpi=300)
	plt.plot(domain, errors2)
	plt.show()

# Test relaxation alone in 2d with a grid of size L
def relax2(L):
	C = 2.0
	k = 4

def mytest2d():
	k = 30
	l = 10
	m = 1

	L = 6
	n = 2**L

	X, Y = np.meshgrid(np.linspace(0,1,n+1)[1:-1],np.linspace(0,1,n+1)[1:-1])
	f = 2*np.pi*np.pi*np.multiply(np.sin(k*np.pi*X), np.sin(l*np.pi*Y))
	u = 0.5*f/(np.pi*np.pi)
	v = f/(2*np.pi*np.pi)

	#v = np.zeros((n-1,n-1))
	#f = np.zeros((n-1,n-1))
	error = np.zeros(m)
	for i in xrange(m):
		v = V_cycle(v, f, L, 2, 2, 'jacobi', 'full', 2, 1)
		error[i] = la.norm(u-v)/2**(L/2.0)
	#plt.plot(np.arange(m), error)
	#plt.show()
	print error

	v = FMG(f, L, 2, 2, 'jacobi', 'full', 2, 1)
	print la.norm(u-v)/2**(L/2.0)

def test2d():
	C = 2.0
	k = 30
	l = 10
	sig = 0
	m = 10

	L = 4
	n = 2**L

	X, Y = np.meshgrid(np.linspace(0,1,n+1)[1:-1],np.linspace(0,1,n+1)[1:-1])
	f = C*np.multiply(np.sin(k*np.pi*X), np.sin(l*np.pi*Y))
	u = f/(np.pi**2*k**2 + np.pi**2*l**2 + sig)
	v = f/C

	#f = np.zeros((n-1,n-1))
	error = np.zeros(m)
	for i in xrange(m):
		v = V_cycle(v, f, L, 2, 2, 'jacobi', 'full', 2, 2)
		error[i] = la.norm(u-v)/2**(L/2.0)
	plt.plot(np.arange(m), error)
	plt.show()
	print error

	v = FMG(f, L, 2, 2, 'jacobi', 'full', 2, 1)
	print la.norm(u-v)/2**(L/2.0)


