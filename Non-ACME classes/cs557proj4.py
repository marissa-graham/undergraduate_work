from __future__ import division
import sys
import numpy as np
from matplotlib import pyplot as plt

def deCasteljua(P,n,u):
    """Better because doesn't require calculating Bernstein schtuff."""
    Q = np.copy(P)
    #print "Q.shape, n: ", Q.shape, n

    # Need each column of Q to be a point
    # P is an array of points; each column is a point
    for k in xrange(1,n+1):
        for i in xrange(n-k+1):
            Q[:,i] = (1.0-u)*Q[:,i] + u*Q[:,i+1]
    return Q[:,0].T

def deCasteljua2(P,n,m,u0,v0):
    """n rows by m columns. MUST BE N ROWS BY M COLUMNS."""
    #print "P.shape: ", P.shape
    #print "n, m = ", n, m
    if n <= m:
        # 3 because 3d points
        Q = np.zeros((3,n+1))
        for j in xrange(n+1):
            # Go through the rows
            Q[:,j] = deCasteljua(P[j,:,:].T,m,u0)
        return deCasteljua(Q,n,v0)
    else:
        #print "Case where m > n:"
        Q = np.zeros((3,m+1))
        for i in xrange(m+1):
            # Go through the columns
            Q[:,i] = deCasteljua(P[:,i,:].T,n,v0)
        return deCasteljua(Q,m,u0)

def bezier_surface(A):
    """
    Plot the nonrational bezier surface defined by the (n+1)x(m+1)
    net of points defined by the (n+1)x(m+1)x3 matrix A.
    """
    n, m, z = A.shape
    n, m = n-1, m-1
    res = 10
    B = np.zeros((res,res,3))

    u = np.linspace(0,1,res)
    v = np.linspace(0,1,res)
    for i in xrange(res):
        for j in xrange(res):
            B[i,j,:] = deCasteljua2(A,n,m,u[i],v[j]) 

    return B

def get_normals(pts, deg):
    # Get the partial derivatives control points in the x direction
    # n*(P_i+1 - P_i)
    D_x = deg[1]*(pts[:,1:,:] - pts[:,:-1,:])
    print D_x.shape

    # Get the partial derivatives control points in the y direction
    D_y = deg[0]*(pts[1:,:,:] - pts[:-1,:,:])
    print D_y.shape

    # Get the partial derivative values at each sample point in the patch
    B_x = bezier_surface(D_x)
    B_y = bezier_surface(D_y)

    # Cross the partial derivatives to get the normal
    return np.cross(B_x,B_y)

class Patch(object):
    def __init__(self, degree, pts):
        # Get the i,jth point with pts[i,j,:]
        # Get the x coordinates with pts[:,:,0], etc.
        # pts[:,i,:] gives you the ith column of points in the patch
        # pts[j,:,:] gives you the jth row of points in the patch
        self.degree = degree
        self.points = pts

def get_patches(filename):
    patches = {}
    with open(filename) as f:
        # Get number of patches
        N = int(f.readline().strip().lower())

        for i in xrange(N):
            # Get the degree of the patch
            deg = map(int, f.readline().strip().lower().split())

            # Read in the points for the patch
            pts = []
            for j in xrange((deg[0]+1)*(deg[1]+1)):
                pts.append(np.array(map(float, f.readline().strip().lower().split())))

            # Each point is a row
            pts = np.array(pts).reshape((deg[0]+1,deg[1]+1,3))
            patches[i] = Patch(deg, pts)

    return N, patches

def main():
    k = 1
    # Ensure command line input is correct
    if len(sys.argv) < 3:
            print "Syntax:\n python cs557proj4.py testfile outputfile"
            return

    # Read in and store patches
    N, patches = get_patches(sys.argv[1])
    
    # Create the output file.
    output = ["g"]
    face_output = []

    # Go through the patches
    for n in xrange(N):
        print patches[n].points.shape
        B = bezier_surface(patches[n].points)
        Vn = get_normals(patches[n].points,patches[n].degree)

        for i in xrange(9):
            for j in xrange(9):
                # Add the face 
                output.append("v {} {} {}".format(*B[i,j,:]))
                output.append("v {} {} {}".format(*B[i+1,j,:]))
                output.append("v {} {} {}".format(*B[i+1,j+1,:]))
                output.append("v {} {} {}".format(*B[i,j+1,:]))
                output.append("vn {} {} {}".format(*Vn[i,j,:]))
                output.append("vn {} {} {}".format(*Vn[i+1,j,:]))
                output.append("vn {} {} {}".format(*Vn[i+1,j+1,:]))
                output.append("vn {} {} {}".format(*Vn[i,j+1,:]))
                face_output.append("f {}//{} {}//{} {}//{} {}//{}".format(k,k,k+1,k+1,k+2,k+2,k+3,k+3))
                k += 4

    # Write to output file    
    with open(sys.argv[2], 'w') as f:
        f.write("\n".join(output)+"\n")
        f.write("\n".join(face_output)+"\n")

main()