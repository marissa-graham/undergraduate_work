# Functions from local_projection.py that do not call any other functions and are not called
# by any other functions
# Also functions that are duplicates of others and not used
from local_projection import *

from matplotlib.transforms import offset_copy 

# Symbolic representation of a bernstein basis
class bernstein_basis_sympy(sympy.Function):
    nargs = 3
    @classmethod
    def eval(cls, n, i, x):
        return sympy.binomial(n,i) * ((x+1)/2) ** i * ((1-x)/2) ** (n-i)

def compute_inv_ext_op_col(interval, support_kv, degree, tol=1e-10):
    extended_kv = []
    mult1=0
    mult2=0
    for k in support_kv:
        if abs(k-interval[0]) < tol:
            mult1+=1
        if abs(k-interval[1]) < tol:
            mult2+=1

def sympy_single_elevation_op(p):
    A = sympy.zeros((p+1,p+2))
    A[0,0] = 1
    A[-1,-1] = 1
    for i in range(1,p+1):
        A[i-1,i] = sympy.Rational(i, p + 1)
        A[i,i] = 1 - sympy.Rational(i, p + 1)
    return A


# elevation_op
def single_elevation_op(p):
    A = np.zeros((p+1,p+2))
    A[0,0] = 1
    A[-1,-1] = 1
    for i in range(1,p+1):
        A[i-1,i] = i / np.float(p + 1)
        A[i,i] = 1 - i / np.float(p + 1)
    return A

def elevation_op(p1, p2):
	if p1 == p2:
		return np.eye(p1+1)
	else:
		A = single_elevation_op(p1)
		for p in range(p1+1, p2):
			A = np.dot(A, single_elevation_op(p))
		return A


# plot_kv, label_kv
def compute_kv_plot_pts(kv, start = (0,0), orient = (1,0), split_repeated = True, split_space = 0.1):
    
    orient = np.array(orient)
    orient = orient / np.linalg.norm(orient)
    start_pt = np.array(start)
    pts = []
    unique_knots = get_unique_knots( kv )
    knot_mults = get_knot_multiplicities( kv )
    for i, k in enumerate(unique_knots):
        if split_repeated:
            mult = knot_mults[ i ]
            if i == 0:
                mult_start = - (mult) * split_space
            elif i == len(unique_knots) - 1:
                mult_start = k - split_space
            else:
                mult_start = k - (mult+1) / 2.0 * split_space
            if mult > 1:
                for m in range(mult):
                    pts.append(start_pt + ( mult_start + (m+1) * split_space ) * orient )
            else:
                pts.append(start_pt + k * orient)
        else:
            pts.append(start_pt + k * orient)
    return pts

def plot_kv(kv, start = (0,0), orient = (1,0), marker = "D", split_repeated = True, split_space = 0.1, **kwargs):
    
    pts=np.array(compute_kv_plot_pts(kv, start, orient, split_repeated = split_repeated, split_space = split_space))
    x=pts[:, 0]
    y=pts[:, 1]
    if kwargs.has_key("axes"):
        ax = kwargs["axes"]
    else:
        ax = plt.gca()
    ax.plot(x,y,marker=marker,**kwargs)
    return pts

def label_kv(kv, start = (0,0), orient = (1,0), offset = (0,-0.1), 
             decimal = True, frac_string = "frac",
             template = "${val}$", zero_based = False,
             label_repeated = True,
             **kwargs):
    """
    The template can use the formatting strings {val} to use the knot value or {ind} to use the integer index.
    """
    
    unique_knots=get_unique_knots(kv)
    knot_mults=get_knot_multiplicities(kv)
    pts=compute_kv_plot_pts(unique_knots, start, orient)
    
    label_pts = np.array([pt + np.array(offset) for pt in pts])
    x=label_pts[:, 0]
    y=label_pts[:, 1]
    
    if kwargs.has_key("axes"):
        ax = kwargs["axes"]
    else:
        ax = plt.gca()
    if zero_based:
        first_knot_ind = 0
    else:
        first_knot_ind = 1
    for i, k in enumerate(unique_knots):
        if decimal:
            knot_string = sympy.latex(str(k))
        else:
            knot_string = sympy.latex(sympy.Rational(str(k)).limit_denominator(1000))
        knot_string.replace("frac", frac_string)
        if zero_based:
            ind = i
        else:
            ind = i + 1
        
        if label_repeated and knot_mults[ i ] > 1:
            strip_temp = template.replace("$","")
            label_string = "$" + ",".join([strip_temp.format(ind = first_knot_ind + j, val = "{val}" ) for j in range(knot_mults[i])]) + "$"
            if label_string.count("{val}") > 0:
                label_string = label_string.format(val = knot_string)
            
        else:
            label_string = template.format(val = knot_string, ind = ind)
        ax.text(x[i], y[i], label_string, **kwargs)
        first_knot_ind += knot_mults[ i ]



def compute_kv_label_pts( pts, offset = (0,-0.1) ):
    
    return [pt + np.array(offset) for pt in pts]


# interval_bernstein_vector
class interval_bernstein_basis_poly(sympy.Function):
    nargs = 5
    def __init__(a, b):
        self.a = a
        self.b = b
    @classmethod
    def eval(cls, a,b, n, i, x):
        return sympy.binomial(n,i) * ((x-a)) ** i * (b-x) ** (n-i) / (b - a ) ** n

class interval_bernstein_vector(sympy.Function):
    nargs = 4
    @classmethod
    def eval(cls, a,b, p, x):
        return sympy.Matrix(p+1,1,lambda i,j: interval_bernstein_basis_poly(a,b,p,i, x))



def build_elem_global_local(connects):
	maps = []
	for c in connects:
		for i, j in enumerate(c):
			maps.append({})
			maps[-1][j] = i
	return maps


# levi_civita
def perm_parity(a):
    '''\
    Using algorithm from http://stackoverflow.com/questions/337664/counting-inversions-in-an-array/6424847#6424847
    But substituting Pythons in-built TimSort'''

    a = list(a)
    b = sorted(a)
    inversions = 0
    while a:
        first = a.pop(0)
        inversions += b.index(first)
        b.remove(first)
    return -1 if inversions % 2 else 1

def levi_civita(dim):
    
    e = np.zeros([dim for i in range(dim)])
    perms = itertools.permutations(range(dim))
    for p in perms:
        e[p] = perm_parity(p)
    return e


class TSplineFunc():
    
    def __init__(self, funcs, pts = None, homogeneous_pts = None, pts_wts = None):
        self.num_bfuncs = len(funcs)
        self.basis = funcs
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
            func_val = self.basis[ifunc](*param_pos)
            if self.wts == None:
                ret_val += func_val[..., np.newaxis] *  np.array(self.pts[ifunc])
            else:
                ret_val += self.wts[ifunc ] * func_val[..., np.newaxis] *  np.array(self.pts[ifunc])
                wt_val += self.wts[ifunc] * func_val

        if self.wts != None:
            ret_val = np.einsum('a...,a...->a...', ret_val, 1 / wt_val)
        if ret_val.shape[0] == 1:
            return ret_val[0]
        else:
            return ret_val


def get_func_start_elem( ifunc, ielem, elems ):
    if ielem == 0:
        return 0
    if elems[ielem].elem_connect.count(ifunc) == 0:
        return None
    while elems[ielem - 1].elem_connect.count(ifunc) > 0:
        ielem = ielem - 1
    return ielem

def get_func_elems( ifunc, elems, startelem = 0 ):
    func_elems = []
    for elem in elems[startelem:]:
        func_on_elem = elem.elem_connect.count(ifunc) > 0
        if not(func_on_elem):
            startelem += 1
        else:
            break
    for elem in elems[startelem:]:
        func_on_elem = elem.elem_connect.count(ifunc) > 0
        if func_on_elem:
            func_elems.append(elem)
        else:
            break
    return func_elems

def leg_eval(coeffs,x):
    
    ret = 0
    for n in range(len(coeffs)):
        ret += coeffs[n] * scipy.special.legendre(n)(x)
    return ret

def offset(ax, x, y):
    return offset_copy(ax.transData, x=x, y=y, units='dots')

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)

def set_ax_lims(pts, ax = None, bound_proportion = 0.05):
    if ax == None:
        ax = plt.gca()

    d_x = np.array(pts)[:,0]
    d_y = np.array(pts)[:,1]
    x_min = d_x.min()
    y_min = d_y.min()
    x_max = d_x.max()
    y_max = d_y.max()
    x_mid = (x_min + x_max) / 2.0
    y_mid = (y_min + y_max) / 2.0
    dx = x_max - x_min
    dy = y_max - y_min
    ax.set_xlim( x_mid - dx * ( 1 + bound_proportion ), x_mid + dx * ( 1 + bound_proportion ) )
    ax.set_ylim( y_mid - dy * ( 1 + bound_proportion ), y_mid + dy * ( 1 + bound_proportion ) )

def plot_single_basis(local_kv, num_pts = 100, ax = None, color_cycle = None, w = 1, **kwargs):

    if ax == None:
        ax = plt.gca()

    nkwargs = kwargs.copy()

    if color_cycle != None:
        nkwargs['color'] = color_cycle.next()
    func = spline_basis(local_kv)
    func_bounds = [local_kv[0],local_kv[-1]]
    t = np.linspace(func_bounds[0],func_bounds[1],num_pts)
    ax.plot(t, func(t) * w, **nkwargs)


class SegElem:
    def __init__(self):
        self.seg_list = []

    def add_seg( self, seg ):
        in_list = False
        for test,i in zip(self.seg_list, itertools.count()):
            if len(seg) != len(test):
                continue
            elif equal( np.array(seg), np.array(test) ):
                in_list = True
                break
        if not in_list:
            self.seg_list.append(seg)

    def add_func_segs( self, func ):
        for seg in func.extraction_rows:
            self.add_seg(seg)

    def get_seg_ind(self, seg):
        for test,i in zip(self.seg_list, itertools.count()):
            if len(seg) != len(test):
                continue
            elif equal( np.array(seg), np.array(test) ):
                return i
            else:
                return -1
    def segN(self):
        return len(self.seg_list)





def bernstein_horner(coeffs, n, t):
    t = t/2.0 + 0.5
    u = 1.0 - t
    bc = 1
    tn = 1
    tmp = coeffs[0]*u;
    for i in range(1, n):
        tn = tn*t
        bc = bc*(n-i+1)/i
        tmp = (tmp + tn*bc*coeffs[i])*u
    return (tmp + tn*t*coeffs[n]);

class ExtractionFunc():
	def __init__(self, degree, ext_rows, elem_lines, default_val = 0.0):
		
		self.degree = degree
		self.ext_rows = ext_rows
		self.elem_lines = elem_lines
		self.min = elem_lines[0]
		self.max = elem_lines[-1]
		def temp_func(s):
			if s > self.max or s < self.min:
				return default_val
			else:
				ielem = 0
				for i in range(len(self.elem_lines)-1):
					if s >= self.elem_lines[i] and s <= self.elem_lines[i+1]:
						ielem = i
						break
					else:
						continue
				t = (s - self.elem_lines[i]) / (self.elem_lines[i+1] - self.elem_lines[i])
				t = 2 * t - 1
				return bernstein_horner(self.ext_rows[ielem], self.degree, t)
		
		self._func = np.vectorize(temp_func)
	def __call__(self, s):
		return self._func(s)
	
	def plot(self, num_pts=100, color = 'r', **kwargs):
		
		s = np.linspace(-1,1,num_pts)
		for i in range(len(self.ext_rows)):
			x = np.linspace(self.elem_lines[i], self.elem_lines[i+1], num_pts)
			f = bernstein_horner(self.ext_rows[i], self.degree, s)
			plt.plot(x, f, color = color, **kwargs)
	