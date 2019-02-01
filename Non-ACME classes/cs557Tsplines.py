from __future__ import division
import sys
import numpy as np
from scipy import misc
from matplotlib import pyplot as plt

"""
Read in a T-spline surface and output a .obj file that contains a tessellation of the surface.
Pretty sure it's a cubic T-spline

TODO: Do we need normals?
"""

class TSpline(object):
    def __init__(self, m, n, rho, s_knots, t_knots, symbol_table, control_pts):
        self.m = m
        self.n = n
        self.rho = rho
        self.s_knots = s_knots
        self.t_knots = t_knots
        self.sym = symbol_table
        self.P = control_pts
        self.local_s_vecs = np.zeros((rho,5))
        self.local_t_vecs = np.zeros((rho,5))

def read_spline(filename):
    """
    File formatting:

    m, n, rho # len(knots in s-direction), len(knots in t-direction), number of control pts
    knot vector in s-direction
    knot vector in t-direction
        symbol table of T-mesh topology (m rows by n)
    list of control points (length rho)
        format: x, y, z, w, i, j
            (x, y, z, w) control point coordinates
            (i, j) index of where the control point lies in the control grid
    """
    with open(filename) as f:
        m, n, rho = map(int, f.readline().strip().lower().split())
        print m, " by ", n, "T-mesh with ", rho, "control points"

        s_knots = np.array(map(float, f.readline().strip().lower().split()))
        print "knots in s-direction: ", s_knots

        t_knots = np.array(map(float, f.readline().strip().lower().split()))
        print "knots in t-direction: ", t_knots
        
        symbol_table = []
        for i in xrange(n):
            symbol_table.append(f.readline().strip().lower())
        print "Symbol table:", symbol_table
        #print symbol_table[3][4]

        control_pts = []

        # DON'T FORGET TO CAST i AND j TO INTS AS NECESSARY
        for i in xrange(rho):
            control_pts.append(map(float,f.readline().strip().lower().split()))
        control_pts = np.array(control_pts)
        print "Control points:", control_pts, control_pts.shape

    return TSpline(m, n, rho, s_knots, t_knots, symbol_table, control_pts)

def get_local_knots(i,spline):
    # Find the local knot vectors for control point Pi using spline properties
    print spline.P.shape
    pt = spline.P[i,:]

    # Both local knot vectors will have five entries, since it's a cubic spline
    local_s, local_t = np.zeros(5), np.zeros(5)
    i0 = int(pt[4])
    j0 = int(pt[5])
    
    local_s[2] = spline.s_knots[i0]
    local_t[2] = spline.t_knots[j0]

    # s direction, left side
    knots_needed = 2
    k = 1
    while knots_needed > 0:
        #if spline.sym[i0-k][j0] == 'o' or spline.sym[i0-k][j0] == '-':
        #print spline.sym[j0][i0-k]

        if spline.sym[j0][i0-k] == 'o' or spline.sym[j0][i0-k]=='|':
            local_s[knots_needed-1] = spline.s_knots[i0-k]
            knots_needed -= 1
        k += 1

    # s direction, right side
    knots_needed = 2
    k = 1
    while knots_needed > 0:
        #print spline.sym[j0][i0+k]
        #if spline.sym[i0+k][j0] == 'o' or spline.sym[i0+k][j0] == '-':
        if spline.sym[j0][i0+k] == 'o' or spline.sym[j0][i0+k]=='|':
            local_s[5-knots_needed] = spline.s_knots[i0+k]
            knots_needed -= 1
        k += 1
    
    # t direction, left side
    knots_needed = 2
    k = 1
    while knots_needed > 0:
        #print spline.sym[j0-k][i0]
        if spline.sym[j0-k][i0] == 'o' or spline.sym[j0-k][i0] == '-':
            local_t[knots_needed-1] = spline.t_knots[j0-k]
            knots_needed -= 1
        k += 1

    # t direction, right side
    knots_needed = 2
    k = 1
    while knots_needed > 0:
        if spline.sym[j0+k][i0] == 'o' or spline.sym[j0+k][i0] == '-':
            local_t[5-knots_needed] = spline.t_knots[j0+k]
            knots_needed -= 1
        k += 1

    #print i0, j0, local_s, local_t
    spline.local_s_vecs[i,:] = local_s
    spline.local_t_vecs[i,:] = local_t

def one_basis_func(p,m,U,i,u):
    """Compute the basis function N_ip.

    p is the degree. U is the knot vector. m is the index of the furthest knot,
    so m = len(U) - 1, getting the ith basis function, at the point u.
    """

    # Special cases
    #if (i == 0 and u == U[0]) or (i == m-p-1 and u == U[m]): 
    #    return 1.0

    # Local property
    if u < U[i] or u >= U[-1]:
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

def blend_function(i,s,t,spline):
    # Evaluate the blending function for control pt Pi at (s,t), given local knot vectors
    
    # Get local knot vectors in s and t direction
    local_s = spline.local_s_vecs[i,:]
    local_t = spline.local_t_vecs[i,:]

    find_s = np.argwhere(local_s<s)
    find_t = np.argwhere(local_t<t)

    return one_basis_func(3,4,local_s,0,s)*one_basis_func(3,4,local_t,0,t)

def eval_surface(s,t,spline):
    # Add up blend_function(i,s,t,spline) for all the i's
    # NOT ALL OF THEM. ONLY THE ONES THAT CONTRIBUTE. WILL BE FASTER.
    # yeah but that's a trickier question for T splines

    # TODO: Weighted or unweighted? Divide out weights in blend_function or here? Try both
    tot = np.zeros(6)
    for i in xrange(spline.rho):
        tot += spline.P[i,:]*blend_function(i,s,t,spline)
    #print tot
    return tot[:3]

def main():
    k = 1
    # Ensure command line input is correct
    if len(sys.argv) < 3:
            print "Syntax:\n python cs557Tsplines.py testfile outputfile"
            return

    # Read in the data file and store the control points and symbol table
    spline = read_spline(sys.argv[1])
    for i in xrange(spline.n):
        print spline.sym[i]
    print spline.s_knots
    print spline.t_knots

    # For each control pt, collect and store local knot vectors
    for i in xrange(spline.rho):
        get_local_knots(i,spline)

    # Draw the T-spline surface by tessellating into a 40x40 grid of quadrilaterals
    # Also this is where you write to output file
    sdom = np.linspace(np.min(spline.s_knots),np.max(spline.s_knots),41)
    tdom = np.linspace(np.min(spline.t_knots),np.max(spline.t_knots),41)
    
    output = ["g"]
    face_output = []

    out_pts = np.zeros((41,41,3))
    for i in xrange(41):
        for j in xrange(41):
            out_pts[i,j,:] = eval_surface(sdom[i],tdom[j],spline)

    for i in xrange(40):
        for j in xrange(40):
            output.append("v {} {} {}".format(*out_pts[i,j,:3]))
            output.append("v {} {} {}".format(*out_pts[i+1,j,:3]))
            output.append("v {} {} {}".format(*out_pts[i+1,j+1,:3]))
            output.append("v {} {} {}".format(*out_pts[i,j+1,:3]))
            face_output.append("f {} {} {} {}".format(k,k+1,k+2,k+3))
            k += 4

    # Write to output file    
    with open(sys.argv[2], 'w') as f:
        f.write("\n".join(output)+"\n")
        f.write("\n".join(face_output)+"\n")
    print "all done"

main()