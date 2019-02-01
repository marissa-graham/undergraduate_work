from mpl_toolkits import mplot3d as a3
#from scipy import linalg as la
import pylab
import numpy as np
from matplotlib import pyplot as plt
import intersect2d
from intersect3d import *

def check_plane(P):
    n = np.cross(P[:,1] - P[:,0], P[:,2] - P[:,1])
    
    for i in xrange(3,P.shape[1]):
        nprime = np.cross(P[:,i-1] - P[:,i-2], P[:,i] - P[:,i-1])
        if np.allclose(nprime,n) or np.allclose(nprime, -1.0*n):
            pass
        else:
            return False
    
    return True

def display_polypoly(i, j):
    print '\n\t', poly_names[i], poly_names[j], '\n'
    P = polygons[i]
    Q = polygons[j]
    
    ax = a3.Axes3D(pylab.figure())
    ax.scatter(P[0,:], P[1,:], P[2,:], c='C0', s=50)
    poly = a3.art3d.Poly3DCollection([P.T]) 
    poly.set_edgecolor('k')
    poly.set_facecolor([0.5,0.5,1.0,0.4])
    ax.add_collection(poly)
    
    
    ax.scatter(Q[0,:], Q[1,:], Q[2,:], c='C0', s=50)
    poly = a3.art3d.Poly3DCollection([Q.T]) 
    poly.set_edgecolor('k')
    poly.set_facecolor([0.2,0.9,0.5,0.4])
    ax.add_collection(poly)
    
    pieces = polygon_intersection(P,Q)
    for piece in pieces:
        #print "\n3d piece:\n",piece, '\n'
        print piece.shape
        if piece.shape[1] == 1:
            ax.scatter(piece[0], piece[1], piece[2], c='purple', marker='*', s=150)
        else:
            ax.scatter(piece[0,:], piece[1,:], piece[2,:], c='purple', marker='*', s=150)
            ax.plot(piece[0,:], piece[1,:], piece[2,:], c='purple',linewidth=5)
            ax.plot( [piece[0,-1],piece[0,0]], [piece[1,-1],piece[1,0]], [piece[2,-1],piece[2,0]], c='purple',linewidth=5)
        
    ax.set_zlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xlim(0,1)
    
    plt.title(poly_names[i] + ' + ' + poly_names[j])
    pylab.show()

def display_polyplane(i, j):
    print '\n\t', poly_names[i], poly_names[j], '\n'
    P = polygons[i]
    Q = polygons[j]

    # Initialize the figure and get the intersection
    ax = a3.Axes3D(pylab.figure())
    pieces = plane_polygon_intersection(P,Q)
    n2 = np.cross(Q[:,1] - Q[:,0], Q[:,2] - Q[:,1])
    n2 = n2/la.norm(n2)
    
    # Display the polygon
    ax.scatter(P[0,:], P[1,:], P[2,:], c='C0', s=50)
    poly = a3.art3d.Poly3DCollection([P.T]) 
    poly.set_edgecolor('k')
    poly.set_facecolor([0.5,0.5,1.0,0.4])
    ax.add_collection(poly)
    
    # Display the plane

    # Get the line between P and Q
    p1, p2 = plane_intersection(P,Q)
    M, Minv = to_2d(P)
    
    # Get the intersections of P with the line between P and Q
    int_coords2d, edges = segment_intersect3d(P, p1, p2)
    thirdrow = np.dot(M,p1)[-1]*np.ones((1,int_coords2d.shape[1]))
    int_coords = np.dot(Minv, np.append(int_coords2d, thirdrow, axis=0))

    line_min3d = np.argmin(int_coords[0,:])
    line_max3d = np.argmax(int_coords[0,:])
    line_min = int_coords[:,line_min3d]
    line_max = int_coords[:,line_max3d]

    d = la.norm(line_max - line_min)
    x = (line_max - line_min)/d

    y = np.cross(x, n2)
    y = y/la.norm(y)

    x1 = line_min - 0.1*x - 0.5*d*y
    x2 = line_max + 0.1*x - 0.5*d*y
    x3 = line_max + 0.1*x + 0.5*d*y 
    x4 = line_min - 0.1*x + 0.5*d*y 

    Qprime = np.array([x1,x2,x3,x4]).T 
    if intersect2d.area2d(Qprime) < 0:
        Qprime = Qprime[:,::-1]

    # Display the normal vector
    mid = np.mean(Qprime,axis=1)
    normal = np.array([mid, mid + 0.25*n2]).T 
    ax.plot(normal[0,:], normal[1,:], normal[2,:])
    ax.scatter(normal[0,1], normal[1,1], normal[2,1],marker='o',s=50)
    
    poly = a3.art3d.Poly3DCollection([Qprime.T]) 
    poly.set_edgecolor('k')
    poly.set_facecolor([0.2,0.9,0.5,0.4])
    ax.add_collection(poly)


    # Display the intersection
    for piece in pieces:
        #print "\n3d piece:\n",piece, '\n'
        #print piece.shape
        if piece.shape[1] == 1:
            ax.scatter(piece[0], piece[1], piece[2], c='purple', marker='*', s=150)
        else:
            ax.scatter(piece[0,:], piece[1,:], piece[2,:], c='purple', marker='*', s=150)
            ax.plot(piece[0,:], piece[1,:], piece[2,:], c='purple',linewidth=5)
            ax.plot( [piece[0,-1],piece[0,0]], [piece[1,-1],piece[1,0]], [piece[2,-1],piece[2,0]], c='purple',linewidth=5)
        
    # Formatting and title
    ax.set_zlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xlim(0,1)
    plt.title(poly_names[i] + ' + ' + poly_names[j])
    pylab.show()
    
simplex = np.eye(3)

tri = np.array([[1,0.6,0],
                [0,0.6,0.5],
                [0,1,0]],dtype=np.float)

quad3 = np.array([[1,1,0,0],
                 [0,1,1,0],
                 [0,0,1,1]],dtype=np.float)

quad = np.array([[0.01,0.9,0.25,0.35],
                 [0.01,0.9,0.25,0.35],
                 [0.01,0.55,0.4,0.95]])

quad2 = np.array([[0.4,0.85,0.5,0.25],
                  [0.05,0.8,0.45,0.5],
                  [0.05,0.9,0.65,1.0]])
                 
trap = np.array([[0.9,0.1,0.1,0.9],
                 [0.9,0.1,0.1,0.9],
                 [0,0,0.9,0.1]])
                 
vert_tri = np.array([[1,1,0],
                     [0,1,1],
                     [0,0,1]],dtype=np.float)
                     
weirdie = np.array([[0.1,0.2,0.8,0.9,0.8,0.6,0.8,0.75,0.5,0.25,0.2,0.4,0.2],
                    [0.7,0.1,0.1,0.7,0.9,0.8,0.7,0.4,0.2,0.4,0.7,0.8,0.9],
                    [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]])

wprime = np.array([[0.2,0.4,0.2,0.25,0.5,0.75,0.8,0.6,0.8,0.9,0.8,0.2,0.1],
                   [0.1,0.2,0.3,0.6,0.8,0.6,0.3,0.2,0.1,0.3,0.9,0.9,0.3],
                   [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]])

spikey = np.array([[0.2,0.8,0.8,0.6,0.55,0.5,0.45,0.4,0.2],
                   [0.2,0.2,0.6,0.6,0.5,0.6,0.5,0.6,0.6],
                   [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]])

spikey2 = np.array([[0.2,0.8,0.8,0.6,0.55,0.5,0.45,0.4,0.2],
                   [0.2,0.2,0.6,0.5,0.5,0.6,0.5,0.6,0.6],
                   [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]])

spikey3 = np.array([[0.2,0.8,0.8,0.6,0.55,0.5,0.45,0.4,0.2],
                   [0.2,0.2,0.5,0.5,0.5,0.6,0.5,0.6,0.6],
                   [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]])
                   
vert_square = np.array([[0,1,1,0],
                        [0.4,0.4,0.8,0.8],
                        [0,0,1,1]])
                        
#              0      1    2      3      4     5       6         7       8       9          10        11        12
polygons = [simplex, tri, quad, quad2, quad3, trap, vert_tri, weirdie, wprime, spikey, vert_square, spikey2, spikey3]
poly_names = ['simplex', 'tri', 'quad', 'quad2', 'quad3', 'trap', 'vert_tri', 'weirdie', 'wprime', 'spikey', 'vert_square',
                'spikey2', 'spikey3']

display_polyplane(0,7)
for pair in [(2,4), (1,4), (7,10), (0,4), (4,10), (0,7), (9,10), (2,5)]:
    display_polyplane(pair[0],pair[1])
for i in xrange(len(polygons)):
    for j in xrange(len(polygons)):
        display_polyplane(i,j)
        #display_polypoly(i,j)
        print ''