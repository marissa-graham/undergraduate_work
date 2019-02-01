from __future__ import division
import numpy as np
from scipy import linalg as la
import intersect2d

"""
Throughout- 
Planes are defined by an array P that contains at least three points
in the plane, assumed to be in counterclockwise order.

All polygons are 3 x n.
"""

def leftright(x, n, Q):
    """
    Signed distance from the plane
    1: on the "right" (same direction as the normal vector)
    0: on the plane
    -1: on the "left" (opposite direction as the normal vector)
    """
    h = np.dot(n, Q[:,0] - x)
    if np.isclose(h, 0):
        return 0
    elif h > 0:
        return -1
    elif h < 0:
        return 1

def to_2d(P):
    """
    Return the matrix that transforms (x,y,z) coordinates into coordinates into
    the basis for the coordinate system for the plane.
    
    The origin is the closest point in P to the actual origin.
    
    The normal vector forms the z direction.
    
    On the plane itself: normal vector direction is zero, can ignore 
    """
    b1 = (P[:,1] - P[:,0])/la.norm(P[:,1] - P[:,0])
    n = np.cross(P[:,1] - P[:,0], P[:,2] - P[:,1])
    b3 = n/la.norm(n)
    b2 = np.cross(b1,b3)
    b2 = b2/la.norm(b2)
    
    A = np.vstack([b1,b2,b3]).T
    return la.inv(A), A

def set_diameter(all_pts):
    """
    Get the maximum distance between any two points in the
    set of points given by all_pts.
    """
    max_d = 0
    for i in xrange(all_pts.shape[1]):
        for j in xrange(i, all_pts.shape[1]):
            d = np.sum((all_pts[:,i] - all_pts[:,j])**2)
            if d > max_d:
                max_d = d
    return np.sqrt(max_d)

def plane_intersection(P, Q):
    """
    Get the line of intersection between the planes P and Q,
    returning two representative points. No guarantees about the points 
    except that they both lie on the line.
    
    Assume they're not parallel or coincident, we already checked for that.
    """
    
    # Take the cross product of the normal vectors to get a vector that's 
    # parallel to the line
    n1 = np.cross(P[:,1] - P[:,0], P[:,2] - P[:,1])
    n2 = np.cross(Q[:,1] - Q[:,0], Q[:,2] - Q[:,1])
    v = np.cross(n1, n2)
    v = v/la.norm(v) 
    
    # Get a point that's actually on the line 
    d1 = np.dot(n1, P[:,0])
    d2 = np.dot(n2, Q[:,0])
    
    # Intersection with the xy-plane
    denom = n1[0]*n2[1] - n1[1]*n2[0]
    if denom != 0:
        x = (n2[1]*d1 - n1[1]*d2) / denom
        y = (n1[0]*d2 - n2[0]*d1) / denom
        p = np.array([x,y,0])
        return p, p + v
        
    # If you can't get the xy-plane, get the yz-plane
    elif n1[1]*n2[2] - n1[2]*n2[1] != 0:
        denom = n1[1]*n2[2] - n1[2]*n2[1]
        y = (n2[2]*d1 - n1[2]*d2) / denom
        z = (n1[1]*d2 - n2[1]*d1) / denom
        p = np.array([0,y,z])
        return p, p + v
      
    # If you can't get either of those, get the xz-plane  
    elif n1[0]*n2[2] - n1[2]*n2[0] != 0:
        denom = n1[0]*n2[2] - n1[2]*n2[0]
        x = (n2[2]*d1 - n1[2]*d2) / denom
        z = (n1[0]*d2 - n2[0]*d1) / denom
        p = np.array([x,0,z])
        return p, p + v
     
    # I have no idea what this case is for   
    else:
        return P[:,0], P[:,1]

def segment_intersect3d(P, p1, p2):
    """
    Find all the intersection points between P and the LINE (not segment)
    from p1 to p2.
    """

    # We know the line from p1 to p2 is in the same plane as P
    # Call to_2d on P and then use intersect2d.segment_intersection
    M, Minv = to_2d(P)
    P2d = np.dot(M,P)
    p1_2d = np.dot(M,p1)
    p2_2d = np.dot(M,p2)
    
    int_coords, edges = intersect2d.segment_intersection(
                            P2d[:-1,:], p1_2d[:-1], p2_2d[:-1], segment=False)
    
    return int_coords, edges

    
def get_distances(num_ints, coords, p1, p2):
    distances = dict()
    d_p1p2 = np.sum((p1-p2)**2)
    
    for i in xrange(num_ints):
        d_p1 = np.sum((coords[:,i]-p1)**2)
        d_p2 = np.sum((coords[:,i]-p2)**2)
        if d_p1 <= d_p2 and d_p2 > d_p1p2:
            distances[i] = -d_p1
        else:
            distances[i] = d_p1
            
    return distances    

def get_segments(P, Q):

    # Get the line of intersection between P and Q (p1 and p2)
    p1, p2 = plane_intersection(P,Q)
    P_M, P_Minv = to_2d(P)
    Q_M, Q_Minv = to_2d(Q)
    
    # Call segment_intersect3d for both P and Q
    # Need to return in 3d since otherwise they'll be different coordinate systems
    
    P_int_coords, edges = segment_intersect3d(P, p1, p2)
    Q_int_coords, edges = segment_intersect3d(Q, p1, p2)
    
    P_num_ints = P_int_coords.shape[1]
    Q_num_ints = Q_int_coords.shape[1]
    
    P_3d = np.append(P_int_coords, np.dot(P_M,p1)[-1]*np.ones((1,P_num_ints)),axis=0)
    P_3d_coords = np.dot(P_Minv, P_3d)
    Q_3d = np.append(Q_int_coords, np.dot(Q_M,p1)[-1]*np.ones((1,Q_num_ints)),axis=0)
    Q_3d_coords = np.dot(Q_Minv, Q_3d)
    
    # Get the signed distance of the intersections from p1
    P_distances = get_distances(P_num_ints, P_3d_coords, p1, p2)
    Q_distances = get_distances(Q_num_ints, Q_3d_coords, p1, p2)
           
    line_verts = []
    from_P = []
    from_Q = []
    
    # Put the intersection points in order of signed distance from p1
    while max(len(P_distances), len(Q_distances)) > 0:

        # If there are no P intersections left, append the next Q intersection
        if len(P_distances) == 0:
            Qnext = min(Q_distances, key=Q_distances.get)
            line_verts.append(Q_3d_coords[:,Qnext])
            del Q_distances[Qnext]
            from_P.append(0)
            from_Q.append(1)

        # If there are no Q intersections left, append the next P intersection
        elif len(Q_distances) == 0:
            Pnext = min(P_distances, key=P_distances.get)
            line_verts.append(P_3d_coords[:,Pnext])
            del P_distances[Pnext]
            from_P.append(1)
            from_Q.append(0)

        else:
            # Look at the next P intersection and Q intersection
            Pnext = min(P_distances, key=P_distances.get)
            Qnext = min(Q_distances, key=Q_distances.get)

            # If they're the same distance, the intersection point occurs on 
            # both polygons
            if np.isclose(P_distances[Pnext], Q_distances[Qnext]):
                line_verts.append(P_3d_coords[:,Pnext])
                del P_distances[Pnext]
                del Q_distances[Qnext]
                from_P.append(1)
                from_Q.append(1)

            # If the P intersection is closer, append that one
            elif P_distances[Pnext] < Q_distances[Qnext]:
                line_verts.append(P_3d_coords[:,Pnext])
                del P_distances[Pnext]
                from_P.append(1)
                from_Q.append(0)

            # Otherwise, append the Q intersection
            else:
                line_verts.append(Q_3d_coords[:,Qnext])
                del Q_distances[Qnext]
                from_P.append(0)
                from_Q.append(1)
    
    # Go through and delete duplicates within tolerance
    i = 0
    while i < len(from_P):
        j = i + 1
        while j < len(from_P):
            if np.allclose(line_verts[i], line_verts[j]):
                del line_verts[i]
                del from_P[i]
                del from_Q[i]
            else:
                j += 1
        i += 1
    
    line_verts = np.array(line_verts).T
    from_P = np.array(from_P)
    from_Q = np.array(from_Q)
    
    return line_verts, from_P, from_Q, p1, p2, P_M, Q_M
    
def delete_duplicates(int_coords, edges, num_ints):
    """
    Delete the duplicate intersections that occur at vertex
    intersections, ensuring that the edge number stored is the edge
    which the vertex is a base of, not an endpoint of.
    """

    i = 0
    new_edges = []

    # Determine which non-unique edge should be stored
    def edge_max(vals):
        if min(vals) + 1 == max(vals):
            return max(vals)
        elif min(vals) == max(vals):
            return vals[0]
        else:
            return 0

    # Walk through the intersections
    while i < num_ints:

        # Keep track of all the edges a given point occurs on
        j = i + 1
        vals = [edges[i]]
        while j < num_ints:
            if np.allclose(int_coords[:,i], int_coords[:,j]):
                vals.append(edges[j])
                num_ints -= 1
                int_coords = np.delete(int_coords, i, axis=1)
                edges = np.delete(edges, j)
            else:
                j += 1
                
        # Store the correct edge    
        new_edges.append(edge_max(vals))
        
        i += 1
        
    return int_coords, np.array(new_edges), num_ints

def get_line_min(int_coords, int_coords2d):
    """
    Get the minimum and maximum points along the line shared
    by int_coords. Return the 2d versions of the points.

    Accounts for lines that change along only one or two axes.
    """

    # Get the minimum and maximum x values
    line_min3d = np.argmin(int_coords[0,:])
    line_max3d = np.argmax(int_coords[0,:])
    line_min = int_coords2d[:,line_min3d]
    line_max = int_coords2d[:,line_max3d]

    # If they are the same, get the minimum and maximum y values
    if np.allclose(line_min, line_max):
        line_min3d = np.argmin(int_coords[1,:])
        line_max3d = np.argmax(int_coords[1,:])
        line_min = int_coords2d[:,line_min3d]
        line_max = int_coords2d[:,line_max3d]

    # If they are still the same, get the min and max z values
    if np.allclose(line_min, line_max):
        line_min3d = np.argmin(int_coords[2,:])
        line_max3d = np.argmax(int_coords[2,:])
        line_min = int_coords2d[:,line_min3d]
        line_max = int_coords2d[:,line_max3d]

    return line_min, line_max
 

def polygon_intersection(P, Q):
    """ Get the intersection between P and Q. """
    
    # Get the normal vectors for each polygon
    n1 = np.cross(P[:,1] - P[:,0], P[:,2] - P[:,1])
    n2 = np.cross(Q[:,1] - Q[:,0], Q[:,2] - Q[:,1])
    n1 = n1/la.norm(n1)
    n2 = n2/la.norm(n2)
    
    pieces = []

    ############# COINCIDENT OR PARALLEL PLANES #####################

    if np.allclose(np.cross(n1,n2),0):
        
        # If they're coincident
        if np.dot(n1, P[:,0]) == np.dot(n2, Q[:,0]):

            # Get the matrix that converts to 2d
            M, Minv = to_2d(P)
            
            # Convert P and Q to 2d and ensure they are counterclockwise
            P2d = np.dot(M, P)
            Q2d = np.dot(M, Q)
            
            if intersect2d.area2d(P2d[:-1,:]) < 0:
                P2d = P2d[:,::-1]
            if intersect2d.area2d(Q2d[:-1,:]) < 0:
                Q2d = Q2d[:,::-1]
                
            # Intersect the 2d polygons
            polygons = intersect2d.polygon_intersection(P2d[:-1,:], Q2d[:-1,:])
            
            # Convert back to 3d
            for polygon in polygons:
                piece = np.zeros((3,polygon.shape[1]))
                piece[:-1,:] = polygon
                piece[-1,:] = P2d[-1,0]
                pieces.append(np.dot(Minv, piece))
            
    ################## SINGLE LINE OF INTERSECTION ######################
    else: 

        # Get the unique points of intersection between P and Q and which
        # polygon they came from
        line_verts, from_P, from_Q, p1, p2, P_M, Q_M = get_segments(P, Q)
        
        # Convert P and Q to 2d and ensure they are counterclockwise
        P2d = np.dot(P_M, P)
        Q2d = np.dot(Q_M, Q)
        
        if intersect2d.area2d(P2d[:-1,:]) < 0:
            P2d = P2d[:,::-1]
        if intersect2d.area2d(Q2d[:-1,:]) < 0:
            Q2d = Q2d[:,::-1]
        
        # Get the midpoints and test whether they are on each polygon    
        num_ints = len(from_P)
        mid_P = -1*np.ones(num_ints + 1)
        mid_Q = -1*np.ones(num_ints + 1)
        
        for i in xrange(1,num_ints):
            midpt = 0.5*(line_verts[:,i-1] + line_verts[:,i])
            mid_P[i] = intersect2d.genIO(np.dot(P_M, midpt)[:-1], P2d[:-1,:])
            mid_Q[i] = intersect2d.genIO(np.dot(Q_M, midpt)[:-1], Q2d[:-1,:])

        # Double check whether each intersection is on each polygon
        for i in xrange(num_ints):
            from_P[i] = intersect2d.genIO(np.dot(P_M, line_verts[:,i])[:-1], P2d[:-1,:])
            from_Q[i] = intersect2d.genIO(np.dot(Q_M, line_verts[:,i])[:-1], Q2d[:-1,:])
        
        # Walk through the intersections
        segment = []
        for i in xrange(num_ints):

            # If both mid_prev and both mid_next are all true, add the intersection
            # point to the segment you currently have
            if mid_P[i] >= 0 and mid_P[i+1] >= 0 and mid_Q[i] >= 0 and mid_Q[i+1] >= 0:
                segment.append(line_verts[:,i])
            
            # If both mid_prev but not both mid_next, finish the segment and append
            elif mid_P[i] >= 0 and mid_Q[i] >= 0:
                segment.append(line_verts[:,i])
                pieces.append(np.array(segment).T)
                segment = []
            
            # If both mid_next but not both mid_prev, start a segment
            elif mid_P[i+1] >= 0 and mid_Q[i+1] >= 0:
                segment = [line_verts[:,i]]
            
            # If you have both intersection points but neither mid_next or mid_prev
            # are both true, you need to add a solo intersection point
            else:
                if from_P[i] >= 0 and from_Q[i] >= 0:
                    pieces.append(np.array([line_verts[:,i]]).T)
                    segment = []
    
    # Return the list of polygons/line segments
    return pieces

def plane_polygon_intersection(P, Q):
    """
    Get the list of polygons that forms the intersection between P and 
    the plane defined by the first three points in Q.
    
    The normal vector for Q is defined to be 
        np.cross(Q[:,1] - Q[:,0], Q[:,2] - Q[:,1])
    and "left" is relative to that (left is opposite direction from the
    normal vector)
    """

    n = P.shape[1]
    n2 = np.cross(Q[:,1] - Q[:,0], Q[:,2] - Q[:,1])
    n2 = n2/la.norm(n2)
    
    # Tag all the vertices as left, right, or on the plane
    verts_LR = np.zeros(n)
    num_left = 0
    left_verts = []
    for i in xrange(n):
        verts_LR[i] = leftright(P[:,i], n2, Q)
        if verts_LR[i] < 0:
            num_left += 1
            left_verts.append(P[:,i])

    ################# TRIVIAL CASES ###############################
    
    # If all the vertices are on the right, the intersection is empty
    if np.all(verts_LR == 1):
        return []
        
    # If all the vertices are on the left or on the plane, return P
    elif np.all(verts_LR <= 0):
        return [P]
    
    # Get the line between P and Q
    p1, p2 = plane_intersection(P,Q)
    M, Minv = to_2d(P)

    ############# SETUP FOR NON-TRIVIAL CASES #####################
    
    # Get the intersections of P with the line between P and Q
    int_coords2d, edges = segment_intersect3d(P, p1, p2)
    num_ints = int_coords2d.shape[1]
    
    # Convert back to 3d
    thirdrow = np.dot(M,p1)[-1]*np.ones((1,num_ints))
    int_coords = np.dot(Minv, np.append(int_coords2d, thirdrow, axis=0))
    
    # Delete duplicate intersections
    int_coords, edges, num_ints = delete_duplicates(int_coords, edges, num_ints)
    
    # If edges are in reverse order, reverse
    if edges[0] > edges[-1]:
        edges = edges[::-1]
        int_coords = int_coords[:,::-1]
    int_coords2d = np.dot(M, int_coords)[:-1,:]

    ################# SEGMENT INTERSECTION CASE ####################
    
    # If all the vertices are either on the right or on the plane, then you'll
    # get the line segments like in the polygon_polygon case
    if np.all(verts_LR >= 0):
        xmin = np.argmin(int_coords[0,:])
        xmax = np.argmax(int_coords[0,:])
        Qprime = np.array([int_coords[:,xmin],int_coords[:,xmax],
                           int_coords[:,xmax]+n2,int_coords[:,xmin]+n2]).T
        return polygon_intersection(P, Qprime)

    ############### POLYGON INTERSECTION CASE #######################
    
    # Convert left-hand vertices to 2d
    left_verts = np.array(left_verts).T
    P2d = np.dot(M, P)
    left_verts2d = np.dot(M, left_verts)[:-1,:]

    # Get the minimuma and maximum points along the line
    line_min, line_max = get_line_min(int_coords, int_coords2d)
    
    # Get the distances to use for the bounding box
    max_d = set_diameter(np.append(left_verts2d, int_coords2d, axis=1))
    dn = 1/la.norm(line_min - line_max)

    # Get the points for the bounding box
    x1 = line_min - max_d*dn*(line_max - line_min)
    x2 = line_max + max_d*dn*(line_max - line_min)

    pt = np.array([line_max[1] - line_min[1], line_min[0] - line_max[0]])
    
    x3 = x2 - max_d*dn*pt
    x3p = P2d[-1,0]*np.ones(3)
    x3p[:-1] = x3
    x4 = x1 - max_d*dn*pt 

    # Make sure the bounding box is on the left
    if leftright(np.dot(Minv,x3p),n2,Q) > 0:
        x3 += 2*max_d*dn*pt
        x4 += 2*max_d*dn*pt

    # Make sure the converted 2d arrays are still clockwise
    Qprime = np.array([x1,x2,x3,x4]).T
    if intersect2d.area2d(Qprime) < 0:
        Qprime = Qprime[:,::-1]
    if intersect2d.area2d(P2d[:-1,:]) < 0:
        P2d = P2d[:,::-1]
    
    # Uncomment these lines if you want to see the bounding box
    #thirdrow = np.dot(M,p1)[-1]*np.ones((1,4))
    #Qprime3d = np.dot(Minv,np.append(Qprime,thirdrow,axis=0))
    #pieces = [Qprime3d]

    # Get the intersection and convert back to 3d
    pieces = []
    pieces2d = intersect2d.polygon_intersection(P2d[:-1,:],Qprime)
    
    for i in xrange(len(pieces2d)):
        piece = np.zeros((3,pieces2d[i].shape[1]))
        piece[-1,:] = P2d[-1,0]
        piece[:-1,:] = pieces2d[i]
        piece = np.dot(Minv, piece)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
        pieces.append(piece)
        
    return pieces

