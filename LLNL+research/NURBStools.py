import numpy
def floor_index(a,x):
    ind = abs(a-x).argmin()
    if a[ind]<=x:
        return ind
    else:
        return ind-1
    
class RecursiveNURBSBasis():
    def __init__(self, degree, knotVector):
        self.degree = degree
        self.knotVector = knotVector
        self.numBasisFuncs = len(self.knotVector)-self.degree-1

    def preVecFunc(self, t, index):
        if self.degree==0:
            if t>= self.knotVector[index] and t<self.knotVector[index+1]:
                return 1.0
            elif t==self.knotVector[index+1] and self.knotVector[index+1]==self.knotVector[-1]:
                return 1.0
            else:
                return 0.0
        else:
            tempDegree = self.degree-1
            tempNURBSBasis = RecursiveNURBSBasis(tempDegree, self.knotVector)
            if self.knotVector[self.degree+index]-self.knotVector[index]==0:
                term1 = 0.0
            else:
                term1 = (t-self.knotVector[index]) / (self.knotVector[self.degree+index]-self.knotVector[index]) * tempNURBSBasis(t, index)

            if self.knotVector[self.degree+index+1]-self.knotVector[index+1]==0:
                term2 = 0.0
            else:
                term2 = (self.knotVector[self.degree+index+1]-t) / (self.knotVector[self.degree+index+1]-self.knotVector[index+1]) * tempNURBSBasis(t, index+1)

            return term1 + term2

    def __call__(self, t, index):
        from numpy import vectorize
        tempFunc = numpy.vectorize(self.preVecFunc)
        return tempFunc(t, index)

class BernsteinBasis():
    def __init__(self, degree):
        self.degree = degree
    def __call__(self, t, i):
        from scipy.misc import comb
        if i>0:
            return 1.0 / 2**self.degree * comb(self.degree, i-1)*(1-t)**(self.degree-(i-1)) * (1+t)**(i-1)
        else:
            return 0.0

def fraction_eye(n):
    import numpy
    import fractions
    mat_list = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append(1)
            else:
                row.append(0)
        mat_list.append(row)
    return numpy.array(mat_list, dtype = fractions.Fraction)

def computeSympyExtractionOperator(U,p):
    from numpy import eye, zeros
    import sympy
    """Compute the extraction operator relating a NURBS basis to the Bernstein polynomials over a single patch.
    
    Arguments:
    - `U`: Knot vector
    - `p`: Curve degree
    """
    def equals(a,b,tol=1e-10):
        return abs(a-b) < tol
    m = len(U)
    a = p + 1
    b = a + 1
    nb = 0
    C = []
    C.append(sympy.eye(p+1))
    while b<m:
        C.append(sympy.eye(p+1))
        i = b
        mult = 1
        while b<m and equals(U[b],U[b-1]):
            mult = mult+1
            b = b+1
        # mult = b - i + 1
        if mult < p:
            numer = U[b-1]-U[a-1]
            alphas = zeros((p+1), dtype = sympy.Rational)
            for j in reversed(range(mult,p+2)):
                alphas[j-mult-1] = sympy.Rational(numer, ( U[a + j -1] - U[a-1] ))
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

# def computeSympyExtractionOperator(U,p):
#     from numpy import eye, zeros
#     import sympy
#     """Compute the extraction operator relating a NURBS basis to the Bernstein polynomials over a single patch.
    
#     Arguments:
#     - `U`: Knot vector
#     - `p`: Curve degree
#     """
#     m = len(U)
#     a = p + 1
#     b = a + 1
#     nb = 1
#     C = []
#     C.append(sympy.eye(p+1))
#     if p == 1:
#         return C
#     while b<m:
#         i = b

#         while b<m and U[b]==U[b-1]: b = b+1
#         mult = b - i + 1
#         if mult < p:
#             C.append(sympy.eye(p+1))
#             numer = U[b-1]-U[a-1]
#             alphas = zeros((p+1), dtype = sympy.Rational)
#             for j in reversed(range(mult,p+1)):
#                 alphas[j-mult-1] = sympy.Rational(numer, ( U[a + j -1] - U[a-1] ))
#             r = p-mult
#             for j in range(1,r+1):
#                 save = r - j + 1
#                 s = mult + j
#                 for k in reversed(range(s+1,p+2)):
#                     alpha = alphas[k-s-1]
#                     C[nb-1][:,k-1] = alpha * C[nb-1][:,k-1] + (1 - alpha) * C[nb-1][:, k-2]
#                 if b<m:
#                     C[nb][save-1:j+save,save-1] = C[nb-1][p-j:p+1,p]
#             nb = nb + 1
#         if b < m:
#             a = b
#             b = b + 1
#     return C

def computeExtractionOperator(U,m,p):
    from numpy import eye, zeros
    """Compute the extraction operator relating a NURBS basis to the Bernstein polynomials over a single patch.
    
    Arguments:
    - `U`: Knot vector
    - `m`: number of knots
    - `p`: Curve degree
    """
    a = p + 1
    b = a + 1
    nb = 0
    C = []
    C.append(eye(p+1))
    def equals(a,b,tol=1e-10):
        return abs(a-b) < tol

    while b<m:
        C.append(eye(p+1))
        i = b
        mult = 1
        while b<m and equals(U[b],U[b-1]):
            mult = mult+1
            b = b+1
        # mult = b - i + 1
        if mult < p:
            numer = U[b-1]-U[a-1]
            alphas = zeros((p+1))
            for j in reversed(range(mult,p+2)):
                alphas[j-mult-1] = numer /float( U[a + j -1] - U[a-1] )
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
    return C

class NURBSElementBasis():
    def __init__(self, elementMin, elementMax, extractionOperator, weights=None):
        from numpy import shape, ones
        self.degree = shape(extractionOperator)[0]-1
        self.extractionOperator = extractionOperator
        self.min = elementMin
        self.max = elementMax
        self.B = BernsteinBasis(self.degree)
        self.dBBasis = BernsteinBasis(self.degree-1)
        self.ddBBasis = BernsteinBasis(self.degree-2)

        if weights==None:
            self.weights = ones(shape(extractionOperator)[0]) # need number of weights
            self.Bspline = True
        else:
            self.weights = weights
            self.Bspline = False
    def N(self,t,index):
        output = 0.0
        for i in range(1,self.degree+2):
            output = output + self.extractionOperator[index-1,i-1] * self.B(2 * (t-self.min) / (self.max-self.min) - 1, i)
        return output

    def W(self,t):
        output = 0.0
        for i in range(1, len(self.weights)+1):
            output = output + self.weights[i-1] * self.N(t, i)
        return output
    
    def dB(self,x,i):
        return self.degree * (self.dBBasis(x,i-1) - self.dBBasis(x, i))

    def dN(self,t,index):
        output = 0.0
        for i in range(1,self.degree+2):
            output = output + 1.0/(self.max-self.min) * self.extractionOperator[index-1,i-1] * self.dB(2 * (t-self.min) / (self.max-self.min) - 1, i)
        return output

    def dW(self,t):
        output = 0.0
        for i in range(1, len(self.weights)+1):
            output = output + self.weights[i-1] * self.dN(t, i)
        return output

    def ddB(self,x,i):
        return self.degree * (self.degree - 1) * ( ( self.ddBBasis(x,i-2) - self.ddBBasis(x,i-1)) - ( self.ddBBasis(x, i-1) - self.ddBBasis(x, i) ) )

    def ddN(self, t, index):
        output = 0.0
        for i in range(1,self.degree+2):
            output = output + 1.0/((self.max-self.min))**2 * self.extractionOperator[index-1,i-1] * self.ddB(2 * (t-self.min) / (self.max-self.min) - 1, i)
        return output
    
    def ddW(self, t):
        output = 0.0
        for i in range(1, len(self.weights)+1):
            output = output + self.weights[i-1] * self.ddN(t, i)
        return output
    
    def __call__(self,t,index):
        if self.Bspline:
            return self.N(t, index)
        else:
            return self.weights[index-1] * self.N(t, index) / self.W(t)

    def diff(self, t, index, diffOrder=1):
        if diffOrder==1:
            if self.Bspline:
                return self.dN(t, index)
            else:
                return self.weights[index-1] * ( self.dN(t, index) / self.W(t) - self.N(t, index) * self.dW(t) / self.W(t)**2)
        elif diffOrder==2:
            if self.Bspline:
                return self.ddN(t, index)
            else:
                return self.weights[index-1] * ( self.ddN(t, index) / self.W(t) - 2 * self.dN(t, index) * self.dW(t) / self.W(t)**2 + self.N(t, index) * (2 * self.dW(t)**2 - self.W(t) * self.ddW(t)) / self.W(t)**3 )

class NURBSBasis():
    def __init__(self, degree, knotVector, weights=None, dvalue = 0.0):
        from numpy import unique, shape, ones
        self.degree = degree
        self.knotVector = knotVector
        self.numBasisFunctions = len(self.knotVector)-self.degree-1
        if weights==None:
            self.weights = ones(self.numBasisFunctions)
        else:
            self.weights = weights
        self.dvalue = dvalue
        self.knotSets = []
        for i in range(len(self.knotVector)-self.degree-1):
            self.knotSets.append(self.knotVector[i:i+self.degree+2])
        self.numBasisFuncs = len(self.knotVector)-self.degree-1
        self.uniqueKnotVector = unique(self.knotVector)
        self.extractionOperator = computeExtractionOperator(self.knotVector,len(self.knotVector)-self.degree, self.degree)
    def __call__(self, t, index):
        import numpy as np
        from numpy import unique, vectorize
        knotSet = self.knotSets[index-1]
        uniqueKnotSet = unique(knotSet)
        numIntervalsSpanned = len(uniqueKnotSet)-1
        startKnot = knotSet[0]
        startUniqueKnotGlobalIndex = floor_index(self.uniqueKnotVector, startKnot)
        localKnotSet = uniqueKnotSet[0:numIntervalsSpanned+1]
        
        def tempfunc(x):
            indexOffset = floor_index(localKnotSet, x)
            elementIndex = startUniqueKnotGlobalIndex + indexOffset
            localWeightSet = self.weights[startUniqueKnotGlobalIndex:startUniqueKnotGlobalIndex+self.degree+1]

            if x==self.knotVector[-1]:
                indexOffset = indexOffset-1
                elementIndex = elementIndex-1
            if indexOffset>=0 and indexOffset<len(uniqueKnotSet)-1 and uniqueKnotSet[indexOffset]<=x and uniqueKnotSet[indexOffset+1]>=x:
                localExtractionOperator = self.extractionOperator[elementIndex]
                N = NURBSElementBasis(uniqueKnotSet[indexOffset], uniqueKnotSet[indexOffset+1], localExtractionOperator, weights=localWeightSet)
                return N(x, index-elementIndex)
            else:
                return self.dvalue
        return numpy.vectorize(tempfunc)(t)

    def diff(self, t, index, diffOrder=1):
        from numpy import unique, vectorize
        knotSet = self.knotSets[index-1]
        uniqueKnotSet = unique(knotSet)
        numIntervalsSpanned = len(uniqueKnotSet)-1
        startKnot = knotSet[0]
        startUniqueKnotGlobalIndex = floor_index(self.uniqueKnotVector, startKnot)
        localKnotSet = uniqueKnotSet[0:numIntervalsSpanned+1]
        def tempfunc(x):
            indexOffset = floor_index(localKnotSet, x)
            elementIndex = startUniqueKnotGlobalIndex + indexOffset
            if x==self.knotVector[-1]:
                indexOffset = indexOffset-1
                elementIndex = elementIndex-1
            if indexOffset>=0 and indexOffset<len(uniqueKnotSet)-1 and uniqueKnotSet[indexOffset]<=x and uniqueKnotSet[indexOffset+1]>=x:
                localExtractionOperator = self.extractionOperator[elementIndex]
                N = NURBSElementBasis(uniqueKnotSet[indexOffset], uniqueKnotSet[indexOffset+1], localExtractionOperator)
                return N.diff(x, index-elementIndex, diffOrder=diffOrder)
            else:
                return self.dvalue
        return numpy.vectorize(tempfunc)(t)

    def greville(self, index):
        knotSet = self.knotSets[index-1]
        greville_knots = knotSet[1:-1]
        return greville_knots.mean()

    def avg_knots(self, index):
        knot_set = self.knotSets[index-1]
        return knot_set.mean()

    
    def getRefineCoeffs(self, kvec, kcount):
        from scipy import zeros
        p = kcount - 2
        # see how many knots we need to add to the front.
        front_count = 0;
        first = kvec[ 0 ]

        for i in range(kcount):
            if kvec[ i ] == first:
                front_count += 1
            else:
                break

        back_count = 0
        last = kvec[ kcount - 1 ]
        for i in range(kcount - 1, -1, -1):
            if kvec[ i ] == last:
                back_count += 1
            else:
                break

        # add the needed knots to create the extended
        # knot vector
        front_add = p + 1 - front_count
        back_add = p + 1 - back_count
        U = []
        for i in range(front_add):
            U.append( kvec[ 0 ] )
        for i in range(kcount):
            U.append( kvec[ i ] )
        for i in range(back_add):
            U.append( kvec[ kcount - 1 ] );
        # create the knot insertion vector.
        X = []
        for i in range(len(U) - 1):
            current = U[ i ]
            next = U[ i + 1 ]
            if current != next:
                X.append( current + ( next - current ) / 2.0 )
        # set the number of non-zero knot intervals which exist in this basis function
        seg_n = len(X)
        # knot insertion parameters
        m = len(U) - 1
	n = m - p - 1
	a = p
	b = m - p
	r = len(X) - 1
        # create a nodes vector with a one in the slot
	# corresponding to the basis function associated with
	# the local knot vector
        Pw = zeros(n+1)
        Pw[ front_add ] = 1.0
        Qw = zeros( n + r + 2)
        Ubar = zeros( len(X) + len(U))
        for j in range(a - p + 1):
            Qw[ j ] = Pw[ j ]
        for j in range(b - 1, n + 1):
            Qw[ j + r + 1 ] = Pw[ j ]
        for j in range(a+1):
            Ubar[ j ] = U[ j ]
        for j in range(b, m+1):
            Ubar[ j + r + 1 ] = U[ j ]
        # generate the coefficients
        i = b + p - 1
        k = b + p + r
        def equals(a,b,tol):
            return abs(a-b) < tol
        for j in range(r, -1, -1):
            while X[ j ] <= U[ i ] and i > a:
                Qw[ k - p - 1 ] = Pw[ i - p - 1 ]
                Ubar[ k ] = U[ i ]
                k -= 1
                i -= 1

            Qw[ k - p - 1 ] = Qw[ k - p ]
            for l in range(1, p+1):
                ind = k - p + l
                alfa = Ubar[ k + l ] - X[ j ]
                if equals( abs( alfa ), 0.0, 1e-8 ):
                    Qw[ ind - 1 ] = Qw[ ind ]
                else:
                    alfa = alfa / ( Ubar[ k + l ] - U[ i - p + l ] )
                    Qw[ ind - 1 ] = alfa * Qw[ ind - 1 ] + ( 1.0 - alfa ) * Qw[ ind ]

            Ubar[ k ] = X[ j ];
            k -= 1
        # cut off the zero entries
        coeffs = []
        for i in range(len(Qw)):
            if equals( Qw[ i ], 0.0, 1e-8 ):
                continue
            coeffs.append( Qw[ i ] )
        return coeffs, seg_n
    
    def getFuncRefineCoeffs(self, index):
        coeffs, seg_n = self.getRefineCoeffs(self.knotVector[index:index + self.degree + 2], self.degree + 2)
        return coeffs
    
class NURBSPatch():
    def __init__(self,degree, controlPoints, knotVector, weights=None):
        self.degree = degree
        self.points = controlPoints
        self.knotVector = knotVector
        self.interval = [knotVector[0], knotVector[-1]]
        self.basis = NURBSBasis(self.degree, self.knotVector, weights = weights)
        if len(self.points)<self.basis.numBasisFuncs:
            print 'Insufficient number of points'
        
    def __call__(self, t):
        output = 0.0
        for i in range(0, len(self.points)):
            output = output + self.points[i] * self.basis(t, i+1)
        return output

    def plot(self, numPoints = 1000):
        from pylab import plot, show
        from numpy import linspace
        t = linspace(self.interval[0], self.interval[1], numPoints)
        plot(t, self(t))

def linear_ode(y0, t_range, num_cells, degree = 3, alpha = 1, num_pts = 200):
    from numpy import ones
    from pylab import gca
    from scipy.optimize import fsolve, newton
    knots = linspace(t_range[0], t_range[1], num_cells + 1)
    knots = list(knots)
    knots.reverse()
    for i in range(degree):
        knots.append(knots[-1])
    knots.reverse()
    for i in range(degree):
        knots.append(knots[-1])
    knotVector = array(knots)
    uniqueKnotVector = unique(knotVector)
    t = linspace(knotVector[0],knotVector[-1],num_pts) #
    weights = ones(len(knotVector) - degree - 1)
    N = NURBSBasis(degree, knotVector, weights = weights)
    C = zeros(N.numBasisFuncs)
    C[0] = y0
    A = zeros((degree, degree))
    b = zeros(degree)
    for i in range(degree):
        for j in range(degree):
            t_i = N.greville(i + 1)
            b[i] = C[0] * (alpha * N(t_i, 1) - N.diff(t_i, 1, 1))
            A[i,j] = (N.diff(t_i, j + 2, 1) - alpha * N(t_i, j + 2))

    x = dot(inv(A),b)
    C[1:degree + 1] = x
    for i in range(degree + 1, N.numBasisFuncs):
        t_next = N.greville(i)
        prev_terms = 0
        for j in range(i - degree, i):
            prev_terms = prev_terms - C[j] * (N.diff(t_next, j + 1, 1) - alpha * N(t_next, j + 1))
        curr_terms = (N.diff(t_next, i + 1, 1) - alpha * N(t_next, i + 1))
        def temp(Ci):
            return Ci * curr_terms - prev_terms
        C[i] = newton(temp, C[i-1])
    y = zeros(shape(t))
    for i in range(N.numBasisFuncs):
        y = y + C[i] * N(t, i + 1)
    plot(t,abs(y))
    gca().set_title(r'solution to $y^\prime = \alpha y$')
    gca().set_ylabel('y')
    gca().set_yscale('log')
    return t, y, N, C

def linear_ode_vec(func, y0, t_range, num_cells, degree = 3, num_pts = 200, tol = 1e-8):
    import numpy
    from numpy import ones, dot
    from numpy.linalg import inv
    from pylab import gca
    from scipy.optimize import fsolve
    knots = linspace(t_range[0], t_range[1], num_cells + 1)
    knots = list(knots)
    knots.reverse()
    for i in range(degree):
        knots.append(knots[-1])
    knots.reverse()
    for i in range(degree):
        knots.append(knots[-1])
    knotVector = array(knots)
    uniqueKnotVector = unique(knotVector)
    weights = ones(len(knotVector) - degree - 1)
    print "Creating basis"
    N = NURBSBasis(degree, knotVector, weights = weights, dvalue = 0.0)
    C = zeros((len(y0), N.numBasisFuncs))
    C[:,0] = y0
    t = []
    print "Building initialization matrix"
    Nmat = zeros((degree, degree))
    t_init = []
    t_inter = (N.knotVector[degree+1] - N.knotVector[degree])
    for icol in range(1, degree + 1):
        t_i = N.greville(icol)
        t_init.append(t_i)
        t.append(t_i)

        for ifunc in range(2, degree + 2):
            Nmat[icol - 1, ifunc - 2] = N.diff(t_i, ifunc, 1)
    Nmat_inv = inv(Nmat)
    t_init = array(t_init)
    
    y_init = odeint(func, y0, t_init)

    left_vec = []
    for i in range(0, degree):
        t_i = N.greville(i+1)
        left_vec.append( func(y_init[i], t_i) - C[:,0] * N.diff(t_i, 1, 1))
    left_vec = array(left_vec)
    
    C[:,1:degree+1] = dot(Nmat_inv, left_vec).T
    nfuncs = N.numBasisFuncs
    for i in range(degree + 1, nfuncs):
        t_next = N.greville(i)
        t.append(t_next)
        def temp(Ci):
            y_i = N(t_next, i + 1) * Ci
            dy_i = N.diff(t_next, i + 1, 1) * Ci
            for ifunc in range(i - degree, i):
                y_i = y_i + N(t_next, ifunc + 1) * C[:,ifunc]
                dy_i = dy_i + N.diff(t_next, ifunc + 1, 1) * C[:,ifunc]
            residual = dy_i - func(y_i, t_next)
            return residual
        temp, infodict, ier, mesg = fsolve(temp, C[:,i-1], full_output = 1, xtol = tol)
        if ier != 1:
            print mesg
        C[:,i] = temp
        if i%10 == 0:
            print "Completed step " + str(i) + " of " + str(nfuncs)
    tc = numpy.array(t)

    t = linspace(t_range[0], t_range[-1], num_pts)
    y = zeros((len(y0), len(t)))
    for j in range(N.numBasisFuncs):
        N_j =  N(t, j + 1)
        for i in range(len(y0)):
            y[i,:] = y[i,:] + C[i,j] * N_j
    
    return t, y, tc, N, C


def collocate_poisson(t, N):
    knotVector = N.knotVector
    collocationPoints = []
    for i in range(2,N.numBasisFuncs):
        collocationPoint = N.greville(i)
        collocationPoints.append(collocationPoint)

    collocationPoints = array(collocationPoints)

    K = zeros((N.numBasisFuncs-2,N.numBasisFuncs-2))
    for i in range(0, N.numBasisFuncs-2):
        for j in range(0, N.numBasisFuncs-2):
            K[i,j] = N.diff(collocationPoints[i - 1], j + 2, diffOrder=2)

    
    f = vectorize(lambda x:1)
    F = f(collocationPoints)
    d = dot(inv(K), F)
    x = zeros(N.numBasisFuncs)
    x[1:-1] = d
    temp = NURBSPatch(3, x, knotVector, weights = weights)
    figure()
    def exactSolution(x):
        return x * (x-knotVector[-1])/2.0
    subplot(211)
    temp.plot()
    plot(t, exactSolution(t))
    subplot(212)
    plot(t, temp(t)-exactSolution(t))
    

if __name__=="__main__":
    from numpy import linspace, array, sin, cos, pi, zeros, unique, vectorize, seterr, dot, shape, exp, ones
    import pylab
    from pylab import plot, show, subplot, figure, scatter
    from numpy.linalg import inv
    from scipy.integrate import odeint
    degree = 3
    knots = list(linspace(0,5,6))
    knots.reverse()
    for i in range(degree):
        knots.append(knots[-1])
    knots.reverse()
    for i in range(degree):
        knots.append(knots[-1])
    knots = array(knots)
    n = NURBSBasis(degree, knots)
    x = linspace(knots[0], knots[-1],300)
    subplot(111)
    for i in range(len(knots) - (degree + 1)):
        plot(x, n(x,i+1))
    show()

