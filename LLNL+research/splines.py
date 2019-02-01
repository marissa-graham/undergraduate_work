import numpy as np
from scipy import linalg as la
#from mayavi import mlab
from matplotlib import pyplot as plt

###### BEZIER SURFACES ########

def horner(a,u0):
    """Compute point u0 on power basis curve with coefficients in a."""
    n = a.shape[0]
    C = 0
    for i in xrange(n-1,-1,-1):
        C = C*u0 + a[i]
    return C

def bernstein(i,n,u):
    """Compute the value of the B_{i,n} at the point u."""
    temp = np.zeros(n+1)
    temp[n-i] = 1.0
    u1 = 1-u
    for k in xrange(1,n+1):
        for j in xrange(n,k-1,-1):
            temp[j] = u*temp[j] + u1*temp[j-1]
    return temp[-1]

def all_bernstein(n,u):
    """Compute value of all n-th degree Bernstein polynomials at u."""
    B = np.empty(n+1)
    B[0] = 1.0
    u1 = 1.0-u
    for j in xrange(1,n+1):
        saved = 0
        for k in xrange(j):
            temp = B[k]
            B[k] = saved+u1*temp
            saved = u*temp
        B[j] = saved
    return B

def bezier_point(P,n,u):
    """Compute the point on an nth degree Bezier curve defined by the array P at u."""
    B = all_bernstein(n,u)
    C = 0.0
    for k in xrange(n+1):
        C += B[k]*P[k]
    return C

def deCasteljua(P,n,u):
    """Better because doesn't require calculating Bernstein schtuff."""
    Q = np.copy(P)
    # P is an array of points; each column is a point
    for k in xrange(1,n+1):
        for i in xrange(n-k+1):
            Q[:,i] = (1.0-u)*Q[:,i] + u*Q[:,i+1]
    return Q[:,0]

def horner2(A,m,n,u0,v0):
    """Compute point on a power basis surface; matrix A, vectors u0,v0."""
    b = np.empty(n+1)
    for i in xrange(n+1):
        b[i] = horner(A[i,:],m,v0)
    return horner(b,n,u0)

def deCasteljua2(P,n,m,u0,v0):
    """n rows by m columns. MUST BE N ROWS BY M COLUMNS."""
    if n <= m:
        # 3 because 3d points
        Q = np.zeros((3,n+1))
        for j in xrange(n+1):
            # Go through the rows
            Q[:,j] = deCasteljua(P[j,:],m,u0)
        return deCasteljua(Q,m,v0)
    else:
        Q = np.zeros((3,m+1))
        for i in xrange(m+1):
            # Go through the columns
            Q[:,i] = deCasteljua(P[:,i].T,n,v0)
        return deCasteljua(Q,m,u0)

def bezier_surface(A,save_img=False):
    """
    Plot the nonrational bezier surface defined by the (n+1)x(m+1)
    net of points defined by the (n+1)x(m+1)x3 matrix A.
    """
    n, m, z = A.shape
    n, m = n-1, m-1
    res = 25
    B = np.zeros((res,res,3))
    u = np.linspace(0,1,res)
    v = np.linspace(0,1,res)
    for i in xrange(res):
        for j in xrange(res):
            B[i,j,:] = deCasteljua2(A,n,m,u[i],v[j])    
    mlab.mesh(B[:,:,0],B[:,:,1],B[:,:,2])
    mlab.points3d(A[:,:,0],A[:,:,1],A[:,:,2],scale_factor=0.5)
    mlab.mesh(A[:,:,0],A[:,:,1],A[:,:,2],color=(0,0,0),line_width=0.2,representation='wireframe')
    if save_img:
    	mlab.savefig('bezier_surface_demo.png',size=((1280,720)))
    mlab.show()

########### NURBS #############

# A2.1: Find the knot span index
def find_span(n,p,u,U):
    """
    Find the knot span for a point u.

    p is the degree of the basis functions.
    U is the knot vector.
    
    I'm confused about what n is. 

    m + 1: total number of knots
    n = m - p - 1; One less than total number of basis
        functions since we're indexing by 0.
        n = len(U) - p - 2

    total number of basis functions: m - p. Why? Think the
    triangle table. Start with m + 1 knots, at each step of
    computation the number decreases by 1, over p steps.

    Illustrate with stuff on p. 54?
    """
    if u == U[n+1]: # guessing we're trying to get the endpoint knot here
        return n
    low, high = p, n+1 # why is p the lowest span it could be in?
        # is it to do with the fact that we always have a lot of
        # repeat knots at the endpoints?
    mid = int(0.5*(low+high))
    while (u < U[mid] or u >= U[mid+1]):
        if u < U[mid]:
            high = mid
        else:
            low = mid
        mid = int(0.5*(low+high)) # is it possible for mid to be a non-integer?
    return mid

# A2.2: Compute nonvanishing basis functions
def basis_funcs(i,u,p,U):
    """
    Compute the nonvanishing basis functions.

    i is the knot span of the point U, with the knot vector U,
    and the degree is p. So there are p+1 repeated knots on 
    each end of U
    """
    # N is an array of length p+1 that stores the output
    # left and right are also arrays of length p+1, I think
    N = np.zeros(p+1)
    left = np.zeros(p+1)
    right = np.zeros(p+1)
    N[0] = 1.0
    for j in xrange(1,p+1):
        left[j] = u - U[i+1-j]
        right[j] = U[i+j] - u
        saved = 0.0
        for r in xrange(j):
            temp = N[r]/(right[r+1] + left[j-r])
            N[r] = saved + right[r+1]*temp
            saved = left[j-r]*temp
        N[j] = saved

    return N

# A2.4: Compute the basis function N_ip
def one_basis_func(p,m,U,i,u):
    """Compute the basis function N_ip.

    p is the degree. U is the knot vector. m is the index of the furthest knot,
    so m = len(U) - 1, getting the ith basis function, at the point u.
    """

    # Special cases
    if (i == 0 and u == U[0]) or (i == m-p-1 and u == U[m]): 
        return 1.0

    # Local property
    if u < U[i] or u >= U[i+p+1]:
        return 0.0

    # Initialize zeroth degree functions
    N = np.zeros(p+1)
    for j in xrange(p+1):
        if u >= U[i+j] and u < U[i+j+1]:
            N[j] = 1.0
        else:
            N[j] = 0.0
    for k in xrange(1,p+1):
        if N[0] == 0.0:
            saved = 0.0
        else:
            saved = ((u-U[i])*N[0])/(U[i+k]-U[i])
        for j in xrange(p-k+1):
            Uleft = U[i+j+1]
            Uright = U[i+j+k+1]
            if N[j+1] == 0.0:
                N[j] = saved
                saved = 0.0
            else:
                temp = N[j+1]/(Uright-Uleft)
                N[j] = saved + (Uright-u)*temp
                saved = (u-Uleft)*temp
    return N[0]

# A3.1: Compute point on a B-spline curve for fixed u value
def curve_point(n,p,U,P,u):
    """P is the control points. One point per column."""
    span = find_span(n,p,u,U)
    N = basis_funcs(span,u,p,U)
    C = np.zeros(2)
    for i in xrange(p+1):
        C += N[i]*P[:,span-p+i]
    return C

# Plot a B-spline curve
def plot_spline(n,p,U,P,numpts,save_img=False,filename=None):
    # P has one column per point
    u = np.linspace(U[0],U[-1],numpts)
    out = np.zeros((2,numpts))
    for i in xrange(numpts):
        out[:,i] = curve_point(n,p,U,P,u[i])

    plt.scatter(P[0,:],P[1,:])
    plt.plot(out[0,:],out[1,:])
    if save_img:
    	plt.savefig(filename,dpi=400)
    plt.show()

# A3.5: Get a surface point
def surface_point(n,m,p,q,U,V,u,v,P):
    """
    Generate a point on a B-spline surface.

    n: Highest index of basis functions in the u-direction.
    m: Highest index of basis functions in the v-direction.
    p: Basis function degree in the u-direction.
    q: Basis function degree in the v-direction.

    U: Knot vector in the u-direction.
    V: Knot vector in the v-direction.
    u: u-coordinate of desired point.
    v: v-coordinate of desired point.
    P: Net of control points. The (i,j)th point is accessed
       like P[i,j,:]. Get x-coordinates like [:,:,0], etc. 
    """
    uspan = find_span(n,p,u,U)
    Nu = basis_funcs(uspan,u,p,U)

    vspan = find_span(m,q,v,V)
    Nv = basis_funcs(vspan,v,q,V)
    #print "v", v
    #print "vspan: ", vspan
    #print "q: ", q

    uind = uspan-p
    S = np.zeros(3)
    for l in xrange(q+1):
        temp = np.zeros(3)
        vind = vspan - q + l
        #print vind
        for k in xrange(p+1):
            temp += Nu[k]*P[uind+k,vind,:]
        S += Nv[l]*temp
    return S    

# Plot a B-spline surface
def plot_surface(n,m,p,q,U,V,P,res,save_img=False):
    B = np.zeros((res,res,3))
    u = np.linspace(U[0],U[n+1],res)
    v = np.linspace(V[0],V[m+1],res)
    for i in xrange(res):
        for j in xrange(res):
            B[i,j,:] = surface_point(n,m,p,q,U,V,u[i],v[j],P)    
    mlab.mesh(B[:,:,0],B[:,:,1],B[:,:,2])
    mlab.points3d(P[:,:,0],P[:,:,1],P[:,:,2],scale_factor=0.25)
    mlab.mesh(P[:,:,0],P[:,:,1],P[:,:,2],color=(0,0,0),line_width=0.2,representation='wireframe')
    if save_img:
    	mlab.savefig('b_spline_surface_demo.png',size=(1280,720))
    mlab.show()

# A4.1: Compute point on rational B-spline curve
def weighted_curve_point(n,p,U,Pw,u):
    span = find_span(n,p,u,U)
    N = basis_funcs(span,u,p,U)
    # Pw should be a n+1 array 
        # Number of control points is same as number of
        # basis functions you end up with.
    Cw = np.zeros(3)
    for i in xrange(p+1):
        Cw += N[i]*Pw[span-p+i]
    C = Cw/Cw[-1]
    return C[:-1]

# Plot a NURBS curve
def plot_NURBS_curve(n,p,U,Pw,res):
    # P has one column per point
    u = np.linspace(U[0],U[-1],numpts)
    out = np.zeros((2,numpts))
    for i in xrange(numpts):
        out[:,i] = weighted_curve_point(n,p,U,Pw,u[i])

    plt.scatter(P[0,:],P[1,:])
    plt.plot(out[0,:],out[1,:])
    plt.show()

# A4.3: Compute point on rational B-spline surface
def weighted_surface_point(n,m,p,q,U,V,u,v,Pw):
    """
    Generate a point on a B-spline surface.

    n: Highest index of basis functions in the u-direction.
    m: Highest index of basis functions in the v-direction.
    p: Basis function degree in the u-direction.
    q: Basis function degree in the v-direction.

    U: Knot vector in the u-direction.
    V: Knot vector in the v-direction.
    u: u-coordinate of desired point.
    v: v-coordinate of desired point.
    Pw: Net of weighted control points (n+1 by m+1).
    """
    uspan = find_span(n,p,u,U)
    Nu = basis_funcs(uspan,u,p,U)

    vspan = find_span(m,q,v,V)
    Nv = basis_funcs(vspan,v,q,V)

    uind = uspan-p
    Sw = np.zeros(4)

    for l in xrange(q+1):
        temp = 0.0
        vind = vspan - q + l
        for k in xrange(p+1):
            temp += Nu[k]*Pw[uind+k,vind]
        Sw += Nv[l]*temp

    S = Sw/Sw[-1]
    return S[:-1]

# Plot a NURBS surface
def plot_NURBS(n,m,p,q,U,V,Pw,res):
    B = np.zeros((res,res,3))
    u = np.linspace(U[0],U[n+1],res)
    v = np.linspace(V[0],V[n+1],res)
    for i in xrange(res):
        for j in xrange(res):
            B[i,j,:] = weighted_surface_point(n,m,p,q,U,V,u[i],v[j],Pw)    
    mlab.mesh(B[:,:,0],B[:,:,1],B[:,:,2])
    mlab.points3d(P[:,:,0],P[:,:,1],P[:,:,2],scale_factor=0.5)
    mlab.show()

# Get localized extraction operators corresponding to a certain knot vector
def extraction_operators(p,U):
    """
    Inputs
    ------
    p: curve degree
    U: knot vector

    Outputs
    -------
    nb: number of elements
    Ce = list of element extraction operators
    """
    """
    ########### STUPID WAY (not working) ##########

    # Helper function for figuring out alphas that I don't feel like vectorizing
    # ZERO INDEXING
    def alpha(i,u,U,k):
        if i < k - p:
            return 1
        else if i < k:
            return (u - U[i]) / (U[i+p] - U[i])
        else:
            return 0

    # Figure out how many and which knots to add to the knot vector
    Ubar = []
    curr = U[0]
    numcurr = 1
    for i in xrange(1,len(U)-1):
        # Still on the same knot
        if U[i] == curr:
            numcurr += 1
        else:
            while p - numcurr > 0:
                Ubar.append((curr,i))
                numcurr += 1
            curr = U[i]
            numcurr = 1
    m = len(Ubar)

    # Build global extraction operator
    n = len(U) - p - 2
    C = np.eye(n)
    for j in xrange(m):
        # Get alpha(1) through alpha(n+j)
        # Need i, u, U, k
        Cj = np.zeros(n+j,n+j+1)
        for i in xrange(n+j):
            alph = alpha(i,Ubar[i][0],U,Ubar[i][1])
            if i < n+j-1:
                Cj[i,i] = alph
            if i > 0:
                Cj[i-1,i] = 1 - alph
        C = np.dot(Cj.T,C.T)

    k = 0
    Ce = []
    while k + p + 1 <= n:
        Ce.append(C[k:k+p+1,k*p:k*p+p+1])
        k += 1

    return k+1, Ce
    """

    ########## Algorithm 1 ##########

    a = p+1
    b = a+1
    nb = 1
    m = len(U)

    
    C = np.eye(m) # HOW BIG

    while b < m:

        # WAIT THESE ALL NEED TO GO IN A LIST
        Cnext = np.eye(m) # HOW BIG
        i = b

        # Count multiplicity of the knot at location b
        while b < m and U[b+1] == U[b]:
            b += 1
        mult = b - i + 1

        alphas = np.zeros(p)
        if mult < p:

            # Use (10) to compute the alphas
            numer = U[b] - float(U[a])
            for j in xrange(p,mult+1,-1):
                alphas[j-mult] = numer / (U[a+j] - U[a])

            r = p - mult
            # Update the matrix coefficients for r new knots
            for j in xrange(1,r+1):
                save = r - j + 1
                s = mult + j

                for k in xrange(p+1,s+2):
                    alpha = alphas[k-s]
                    
                    # WHAT IS HAPPENING
                    C[:,k] = alpha*C[:,k] + (1.0-alpha)*C[:,k-1]

                if b < m:
                    # WHAT IS HAPPENING HERE
                    Cnext[save:j+save,save] = C[p-j+1:p+1,p+1]

            nb += 1

            if b < m:
                a = b
                b = b+1

    # Ce needs to be a list of all the things
    return nb, Ce


########## GEOMETRIC ALGORITHMS ########################

# A5.1: Knot insertion
def curve_knot_ins(n,p,UP,Pw,u,k,s,r):
    """
    Compute new curve using knot insertion

    INPUT
    -----
    n: Index of last control point
        DON'T CALL IT NP. NAMESPACE ERRORS.
    p: Degree of the curve
    UP: Knot vector before insertion
    Pw: Control points before insertion. 
        Each point is a column, so the i-th point is Pw[:,i]
    u: Knot to be inserted
    k: Index of knot to insert u after
    s: Initial multiplicity of u
    r: Number of times u is to be inserted

    OUTPUT
    ------
    nq: Index of last control point after insertion
    UQ: New knot vector
    Qw: New control points
    """
    d = Pw.shape[0]

    # Last index of knot vector before insertion
    mp = n + p + 1
    nq = n + r

    # Create new knot vector
    UQ = np.zeros(mp + r + 1)
    for i in xrange(k+1):
        UQ[i] = UP[i]
    for i in xrange(1,r+1):
        UQ[k+i] = u
    for i in xrange(k+1,mp + 1):
        UQ[i+r] = UP[i]

    # Copy over unaltered control points
    Qw = np.zeros((d,n+2)) # YOU JUST CHANGED THIS LINE JUST NOW
    Rw = np.zeros((d,p+1))
    for i in xrange(k-p+1):
        Qw[:,i] = Pw[:,i]
    for i in xrange(k-s,n+1):
        Qw[:,i+r] = Pw[:,i]
    for i in xrange(p-s+1):
        Rw[:,i] = Pw[:,k-p+i]

    # Insert the knot r times
    for j in xrange(1,r+1):
        L = k - p + j
        for i in xrange(p-j-s+1):
            alpha = (u - UP[L+i])/(UP[i+k+1] - UP[L+i])
            Rw[:,i] = alpha*Rw[:,i+1] + (1.0 - alpha)*Rw[:,i]
        Qw[:,L] = Rw[:,0]
        Qw[:,k+r-j-s] = Rw[:,p-j-s]

    # Copy over remaining control points
    for i in xrange(L+1, k-s):
        Qw[:,i] = Rw[:,i-L]

    return nq, UQ, Qw

# A5.8: Curve knot removal (DO)
def remove_curve_knot(n,p,U,Pw,u,r,s,num,TOL=1e-10):
    """
    Attempt to remove the knot u at index r num times.
    Note that index r represents the LAST instance of the knot.
    So if we have the knot vector {0,0,0,1,1,3,3,3}, and we want
    to remove 1, r = 4, not 3.

    INPUT
    -----
    n: Index of furthest control point
    p: Degree of curve
    U: Knot vector of length n + p + 2
    Pw: d x (n+1) matrix of control points
        i-th control point is Pw[:,i]
    u: Knot to be removed
    r: Index of knot to be removed
    s: Multiplicity of knot to be removed
    num: Number of knot removals to be attempted
    TOL: Tolerance level for coincident knots

    OUTPUT
    ------
    t: Number of successful knot removals
    U: Knot vector after removals
    Pw: Control points after removals 
    """
    d = Pw.shape[0]
    m = n + p + 1 #furthest knot vector index
    ord = p + 1

    fout = (2*r - s - p)/2. # first control point out
    last = r - s
    first = r - p
    temp = np.zeros((d,2*p+1))

    for t in xrange(num):
        
        off = first - 1 
        temp[:,0] = Pw[:,off]
        temp[:,last+1-off] = Pw[last+1]

        i, j = first, last
        ii, jj = 1, last - off

        remflag = 0

        while j - i > t:
            # Compute control points for one removal step
            alf_i = (u - U[i]) / (U[i+ord+t] - U[i])
            alf_j = (u - U[j-t]) / (U[j+ord] - U[j-t])

            temp[:,ii] = (Pw[:,i] - (1.0-alf_i)*temp[:,ii-1]) / alf_i
            temp[:,jj] = (Pw[:,j] - alf_j*temp[:,jj+1]) / (1.0-alf_j)

            i, ii = i+1, ii+1
            j, jj = j-1, jj-1

        # Check if knot is removable
        if j - i < t:
            if la.norm(temp[:,ii-1] - temp[:,jj+1]) <= TOL:
                remflag = 1
        else:
            alf_i = (u - U[i]) / (U[i+ord+t] - U[i])
            if la.norm(Pw[:,i] - (alf_i*temp[:,ii+t+1] + (1.0-alf_i)*temp[:,ii-1])) <= TOL:
                remflag = 1

        if remflag == 0: # Cannot remove any more knots
            break # Get out of the for loop
        else:
            # Successful removal. Save new control points.
            i, j = first, last
            while j - i > t:
                Pw[:,i] = temp[:,i-off]
                Pw[:,j] = temp[:,j-off]

                i += 1
                j -= 1

        first -= 1
        last += 1

    if t == 0:
        return t, U, Pw

    for k in xrange(r+1,m+1):
        U[:,k-t] = U[:,k] # Shift knots

    j = fout
    i = j # Going to overwrite Pj through Pi

    for k in xrange(1,t):
        if k % 2 == 1:
            i += 1
        else:
            j -= 1

    for k in xrange(i+1,n+1): # Shift
        Pw[:,j] = Pw[:,k]
        j += 1

    return t, U, Pw

################## TEST CODE ############################

# Test the bezier surface code
def prob4():
    u = np.linspace(0,1,400)
    plt.scatter(-1-u+2*u**2, -2*u+u**2,s=1)
    plt.show()
    plt.scatter(-1+2*u-u**2, -2-u+2*u**2,s=1)
    plt.show()
def test_bernstein(n):
    domain = np.linspace(0,1,100)
    vals = np.zeros((n+1,100))
    for i in xrange(100):
        vals[:,i] = all_bernstein(n, domain[i])
    for i in xrange(n+1):
        plt.plot(domain, vals[i,:])
    plt.savefig('bernstein_demo.png',dpi=400)
    plt.show()
def test_bezier():
    n = 3
    P = np.random.rand(2,n+1)
    #P = np.sort(P,axis=1)
    
    y = np.zeros((2,100))
    x = np.zeros((2,100))
    for i in xrange(100):
        x[:,i] = deCasteljua(P,n,0.01*i)
        #y[:,i] = bezier_point(P,n,0.02*i)
    plt.plot(x[0,:],x[1,:])
    plt.plot(P[0,:],P[1,:])
    #plt.plot(y[0,:],y[1,:])
    #plt.show()
    plt.savefig('bezier_demo.png',dpi=400)
    plt.show()
#test_bernstein(3)
#print "Testing Bezier"
#test_bezier()

def test_surface(): 
    A = np.array([[[0,0,0],[0,2,2],[0,4,0]], 
                  [[3,0,3],[3,2,6],[3,4,3]],
                  [[6,0,3],[6,2,5],[6,4,3]],
                  [[9,0,0],[9,2,2],[9,4,0]]],dtype=np.float);  
    bezier_surface(A,save_img=True)
#test_surface()
def test_rational1d():
    P = np.array([[1,1,0],[0,1,2],[1,1,2]],dtype=np.float)
    y = np.zeros((3,100))
    for i in xrange(100):
        y[:,i] = deCasteljua(P,2,0.01*i)
    y[0,:] /= y[2,:]
    y[1,:] /= y[2,:]
    plt.plot(y[0,:],y[1,:])
    plt.show()
    """ Note: Projecting is not difficult. You just do 
    everything the normal way and at the end divide
    all of the points in the array by their last entry.
    Doesn't require a special function. Gee that must be
    why they do it this way. 3-d shouldn't be any harder
    because the deCasteljua algorithm doesn't care about
    the dimension of the points it's given; the only reason
    surfaces are harder is because it's a NET of points."""


def test_span_finder():
    # n = # of knots - p+1 (takes care of endpoints)
    #     - 1 (takes care of the zero indexing)
    # find_span(n,p,u,U)

    # p = 2, m = 5, n = 5 - 2 - 1 = 2
    u1 = np.array([0.,0,0,1,1,1])

    # p = 3, m = 10, n = 10 - 3 - 1 = 6
    u2 = np.array([0,0,0,0,0.4,0.6,0.6,1,1,1,1])

    # p = 2, m = 10, n = 10 - 2 - 1 = 7
    u3 = np.array([0,0,0,1,2,3,3,4,5,5,5])

   
    print "Testing knot vector", u1
    for u in [0,0.5,1.]:
        print "u = ", u, " gives ", find_span(2,2,u,u1)

    print "\nTesting knot vector", u2
    for u in [0,0.2,0.4,0.5,0.6,0.75,1.]:
        print "u = ", u, " gives ", find_span(6,3,u,u2)

    print "\nTesting knot vector", u3
    for u in [0.,0.5,1,1.5,2.,2.5,3.,3.5,4,4.5,5]:
        print "u = ", u, " gives ", find_span(7,2,u,u3)
def test_curve_plot():
    # def plot_spline(n,p,U,P,numpts)

    # p = 3, m+1 = 12, m = 11, n = m - p - 1 = 7
    U = np.array([0,0,0,0,0.25,0.5,0.75,0.75,1,1,1,1])

    P = np.array([[0,-2,3,4,4,7,9,6],[0,2,4,3,0,0,2,5]])

    plot_spline(7,3,U,P,50,save_img=True)

    # p = 2, m+1 = 11, m = 10, n = 10 - 2 - 1 = 7
    U = np.array([0,0,0,1,2,3,4,5,5,5],dtype=np.float)
    P = np.array([[-6,-5,-3,-1,0,3,2.5,1],[-1,2,3,2,0,1,3,5]])
    #plot_spline(7,2,U,P,100)

    U = np.array([0,0,0,1,1.25,2,3,4,5,5,5],dtype=np.float)
    P = np.array([[-6,-5,-30/8.,-22/8.,-1,0,3,2.5,1],[-1,2,21/8.,23/8.,2,0,1,3,5]])
    #plot_spline(8,2,U,P,100)

def get_figures():
	# Figure 3.13
	U = np.array([0,0,0,0.25,0.5,0.75,1,1,1])
	P = np.array([[0,1,2,2,3,4],[4,0,3,3,0.5,5]])
	plot_spline(5,2,U,P,50,save_img=True,filename='cusp_spline.png')

	# Figure 3.8
	U1 = np.array([0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1])
	U2 = np.array([0,0,0,1/8.,0.25,3/8.,0.5,5/8.,0.75,7/8.,1,1,1])
	P = np.array([[0,2,3,6,10,12,16,18,22,24],[4,6,4,9,8.5,-1,4,2.5,7.5,6]])

	plot_spline(9,9,U1,P,50,save_img=True,filename='ninth_degree_bezier.png')
	plot_spline(9,2,U2,P,50,save_img=True,filename='quadratic_bezier.png')
#get_figures()

def test_surface_plot():
    
    # def plot_surface(n,m,p,q,U,V,P,res)

    # p = 2, m+1 = 8, n = 7 - p - 1 = 4
    U = np.array([0,0,0,0.4,0.6,1,1,1])

    # q = 2, m+1 = 11, m = 10 - p - 1 = 7
    V = np.array([0,0,0,0.2,0.4,0.6,0.7,0.8,1,1,1])

    # Need a 5 * 8 control point net.
    P = np.array([[[0,1,0],[1,1,3],[2,1,3],[3,1,5],[4,1,4],[5,1,2],[6,1,3],[7,1,0]],
                  [[0,2,1],[1,2,2],[2,2,5],[3,2,4],[4,2,5],[5,2,6],[6,2,5],[7,2,2]],
                  [[0,3,2],[1,3,3],[2,3,5],[3,3,3],[4,3,5],[5,3,6],[6,3,5],[7,3,4]],
                  [[0,4,1],[1,4,2],[2,4,4],[3,4,3],[4,4,4],[5,4,5],[6,4,6],[7,4,3]],
                  [[0,5,0],[1,5,1],[2,5,3],[3,5,5],[4,5,4],[5,5,8],[6,5,5],[7,5,2]]],dtype=np.float);
    plot_surface(4,7,2,2,U,V,P,50,save_img=True)
    P = P.transpose((1,0,2))
    #plot_surface(7,4,2,2,V,U,P,50)

    '''
    # This is working

    # If points is a 4x4, and degree is 2, n = 3 + 1 + p = 6
    # p = 2, m+1 = 8, n = 7 - p - 1
    # n = len(U) - p - 2
    # So size of grid = len(U) - p - 1
    # Thus len(U) = size of grid + p + 1
    U = np.array([0,0,0,0.5,1,1,1])
    V = U

    P = np.array([[[0,0,1],[1,0,1],[2,0,1],[3,0,1]],
                  [[0,1,2],[1,1,3],[2,1,3],[3,1,1]],
                  [[0,2,1],[1,2,3],[2,2,3],[3,2,2]],
                  [[0,3,1],[1,3,1],[2,3,1],[3,3,1]]],dtype=np.float);

    plot_surface(3,3,2,2,U,V,P,50)
    '''
#test_surface_plot()


def test_knot_insertion():
    # Create curve


    # p = 3, m+1 = 12, m = 11, n = m - p - 1 = 7
    p = 3
    m = 11
    n = 7
    U = np.array([0,0,0,0,0.25,0.5,0.75,0.75,1,1,1,1])
    P = np.array([[0,-2,3,4,4,7,9,6],[0,2,4,3,0,0,2,5]])
    print U.shape
    print P.shape
    print P
    plot_spline(n,p,U,P,50,save_img=False)
    
    u = 0.35
    k = 4
    s = 0
    r = 1
    nq, UQ, Q = curve_knot_ins(n,p,U,P,u,k,s,r)
    print nq
    print UQ.shape
    print UQ
    print Q.shape
    print Q
    plot_spline(nq,p,UQ,Q,50,save_img=False)

    r = 5
    s = 1
    num = 1
    t, U, P = remove_curve_knot(nq,p,UQ,Q,u,r,s,num)
    print t
    print U
    print P
    plot_spline(n,p,U,P,50,save_img=False)

    #def curve_knot_ins(n,p,UP,Pw,u,k,s,r):
    """
    Compute new curve using knot insertion

    INPUT
    -----
    n: Index of last control point
        DON'T CALL IT NP. NAMESPACE ERRORS.
    p: Degree of the curve
    UP: Knot vector before insertion
    Pw: Control points before insertion. 
        Each point is a column, so the i-th point is Pw[:,i]
    u: Knot to be inserted
    k: Index of knot to insert u after
    s: Initial multiplicity of u
    r: Number of times u is to be inserted

    OUTPUT
    ------
    nq: Index of last control point after insertion
    UQ: New knot vector
    Qw: New control points
    """

    #def remove_curve_knot(n,p,U,Pw,u,r,s,num,TOL=1e-10):
    """
    Attempt to remove the knot u at index r num times.
    Note that index r represents the LAST instance of the knot.
    So if we have the knot vector {0,0,0,1,1,3,3,3}, and we want
    to remove 1, r = 4, not 3.

    INPUT
    -----
    n: Index of furthest control point
    p: Degree of curve
    U: Knot vector of length n + p + 2
    Pw: d x (n+1) matrix of control points
        i-th control point is Pw[:,i]
    u: Knot to be removed
    r: Index of knot to be removed
    s: Multiplicity of knot to be removed
    num: Number of knot removals to be attempted
    TOL: Tolerance level for coincident knots

    OUTPUT
    ------
    t: Number of successful knot removals
    U: Knot vector after removals
    Pw: Control points after removals
    """
    # Plot
    # Remove the inserted knot
    # Plot
    # Insert a repeated knot
    # Plot
    # Remove the repeated knot
    # Plot
#test_knot_insertion()