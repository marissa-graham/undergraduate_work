from __future__ import division
import pylab
import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg as la

"""
All polygons are defined by an array P that contains at least three points, 
assumed to be in counterclockwise order (non-counterclockwise order breaks everything).
"""

def area2d(P):
	"""
	Find the area of the (not necessarily convex) polygon defined
	by the points in the 2 x n array P, using the Green's theorem
	formula.
	"""
	
	n = P.shape[1]
	
	tot = 0
	for i in xrange(n):
		# Inefficiency: Instead of doing a modulus, we could do an
		# if i + 1 == n, j = 0
		j = (i + 1) % n
		tot += P[0,i]*P[1,j] - P[0,j]*P[1,i]
		
	return 0.5*tot

def leftright(x, p1, p2):
	"""
	Determine whether the point x is on the left or right of the line
	going from p1 to p2 by determining if the triangle p1, p2, x goes
	clockwise or counterclockwise.
	
	Return 1 if on the left (counterclockwise)
	Return 0 if on the line
	Return -1 if on the right (clockwise)
	"""
	
	pts = [p1, p2, x]
	
	# Inefficiency: putting pts in list and then numpy array and transposing
	# could all be avoided
	tot = area2d(np.array(pts).T)
	
	if tot == 0:
		return 0
	elif tot > 0:
		return 1
	else:
		return -1
	
def genIO(x, P):
	"""
	Determine whether the 2-d point x is inside or outside of the
	polygon P, not necessarily convex.
	
	P is assumed to be 2 x n.
	Return 1 if inside
	Return 0 if on the polygon
	Return -1 if outside
	"""
	n = P.shape[1]
	num_intersections = 0
	tol = 1e-14
	
	# For each edge of the polygon, check whether the edge intersects
	# the vertical line passing through x and lies below x.
	for i in xrange(n):
		if np.allclose(x, P[:,i]):
			return 0
			
		# Name the points for convenience purposes
		x1 = P[0,i]
		x2 = P[0,(i+1)%n]
		y1 = P[1,i]
		y2 = P[1,(i+1)%n]
		
		# If the line is vertical and x is on it
		if abs(x1 - x[0]) < tol and abs(x2 - x[0]) < tol:
			if x[0] <= max(y1,y2) + tol and x[0] >= min(y1,y2) - tol:
				return 0
			else:
				return -1
		
		# Check if x is inside the projection of the line segment
		# onto the x axis
		if x[0] > min(x1,x2) and x[0] <= max(x1,x2):
			
			# Get the equation of the line
			on_line = y1 + (y2 - y1)/(x2 - x1)*(x[0] - x1)
			
			# If x is on the line
			if np.isclose(x[1], on_line):
				return 0
			
			# If x is above the line
			if x[1] > on_line:
				num_intersections += 1
		
	
	if num_intersections % 2 == 0:
		return -1
	else:
		return 1
		
def line_intersection(p1, p2, p3, p4, segment=True):
	"""
	Find the intersection point between the line segments
	from p1 to p2 and from p3 to p4.
	"""
	
	x1, y1 = p1
	x2, y2 = p2
	
	x3, y3 = p3
	x4, y4 = p4
	
	denom = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4)
	
	# If the line segments are literally the same
	if np.allclose(p1,p3) and np.allclose(p2,p4):
		return p1

	elif np.allclose(p1,p4) and np.allclose(p2,p3):
		return p1

	# Still should probably handle this case
	elif np.isclose(denom, 0):
		if np.isclose(x1,x2):
			# The line from p1 to p2 is vertical
			pass
		elif np.isclose(x3,x4):
			pass
		
	else:
		x = ((x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*y4 - y3*x4)) / denom
		y = ((x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*y4 - y3*x4)) / denom
		
		tol = 1e-14
		if x >= min(x1,x2) - tol and x <= max(x1, x2) + tol: 
			if y >= min(y1,y2) - tol and y <= max(y1,y2) + tol:
				if segment == False:
					return np.array([x,y])
				else:
					if x >= min(x3,x4) - tol and x <= max(x3, x4) + tol: 
						if y >= min(y3,y4) - tol and y <= max(y3,y4) + tol:
							return np.array([x,y])
	
	return []
			
def sort_intersections(int_coords, int_edges, p1, p2, segment=True, return_distances = False):
	"""
	Sort the points in int_coords by signed distance from p1.
	"""
	n = len(int_edges)
	distances = np.zeros(n)
	
	if segment:
		for i in xrange(n):
			distances[i] = np.sum((int_coords[i] - p1)**2)
	else:
		d_p1p2 = np.sum((p1-p2)**2)
		for i in xrange(n):
			d_p1 = np.sum((int_coords[i]-p1)**2)
			d_p2 = np.sum((int_coords[i]-p2)**2)
			if d_p1 <= d_p2 and d_p2 > d_p1p2:
				distances[i] = -d_p1
			else:
				distances[i] = d_p1
				
	sort_inds = np.argsort(distances)
	new_int_coords = np.zeros((2,n))
	new_int_edges = np.zeros(n)
	for i in xrange(n):
		new_int_coords[:,i] = int_coords[sort_inds[i]]
		new_int_edges[i] = int_edges[sort_inds[i]]
		
	if return_distances:
		return new_int_coords, new_int_edges, distances
		
	return new_int_coords, new_int_edges

def segment_intersection(P, p1, p2, segment=True):
	"""
	Get the intersections between the polygon P and the line segment
	from p1 to p2.
	"""
	
	int_coords = []
	int_edges = []
	n = P.shape[1]
	
	for i in xrange(n):
		j = (i+1) % n
		int_pt = line_intersection(P[:,i], P[:,j], p1, p2, segment=segment)
		if len(int_pt) > 0:
			int_coords.append(int_pt)
			int_edges.append(i)
			
	return sort_intersections(int_coords, int_edges, p1, p2, segment)
	
def delete_duplicates(int_coords, P_edges, Q_edges):
	"""
	Delete the duplicate intersection points that occur when an 
	intersection occurs at a vertex.
	"""
	
	# Iterate through the intersection points
	i = 0
	new_P_edges = []
	new_Q_edges = []

	# Rule for picking which edge to keep:
	# Pick the biggest, unless diff > 1 and there's a zero
	def edge_max(P_vals):
		if min(P_vals) + 1 == max(P_vals):
			return max(P_vals)
		elif min(P_vals) == max(P_vals):
			return P_vals[0]
		else:
			return 0
	
	while i < len(int_coords):
		j = i + 1
		P_vals = [P_edges[i]]
		Q_vals = [Q_edges[i]]
		del_inds = []
		
		# Get all the corresponding edge points while deleting duplicates
		# Can only have two values for each
		while j < len(int_coords):
			if np.allclose(int_coords[i], int_coords[j]):
				P_vals.append(P_edges[j])
				Q_vals.append(Q_edges[j])
				
				del int_coords[j]
				del Q_edges[j]
				del P_edges[j]
			else:
				j += 1
		
		new_P_edges.append(edge_max(P_vals))
		new_Q_edges.append(edge_max(Q_vals))
		
		i += 1
		
	return int_coords, np.array(new_P_edges), np.array(new_Q_edges)
	
def get_intersections(P, Q):
	"""
	Get all the points where the edges in P and the edges in Q intersect.
	"""
	n = P.shape[1]
	m = Q.shape[1]
	int_coords = []
	P_edges = []
	Q_edges = []
	
	for i in xrange(m):
		i_ints, Pi_edges = segment_intersection(P, Q[:,i], Q[:,(i+1)%m])
		int_coords.extend(i_ints.T)
		P_edges.extend(Pi_edges)
		for j in xrange(len(Pi_edges)):
			# Append the next edge if it's a vertex
			if np.allclose(np.array(i_ints)[:,j], Q[:,(i+1)%m]):
				Q_edges.append( 1.0*((i+1)%m) )
			# Otherwise append the edge
			else:
				Q_edges.append(i*1.0)
	
	int_coords, P_edges, Q_edges = delete_duplicates(int_coords, P_edges, Q_edges)
	
	return np.array(int_coords).T, np.array(P_edges), np.array(Q_edges)
	
def same_edge(p1, p2, p3, p4):
	""" 
	Return True if p1 == p3 and p2 == p4
	"""
	
	if np.allclose(p1, p3) and np.allclose(p2, p4):
		return True
	else:
		return False
	
def move_along(P, v, n, i, Pi, only_P, P_edges, int_coords):
	""" 
	Move along the polygon P, starting from the point v, which is either
	a vertex or an intersection point or both.
	
	Return:
	
	vnext - Next point in that direction
	i - New intersection index
	Pi - New edge index
	only_P - Return True if we are not at any intersection point.
	"""
	# Check if there are any intersection points on the upcoming edge
	edge_ints = np.where(P_edges == Pi)[0]
	
	# And if any of them are forward
	forward_ints = []
	distances = []
	tol = 1e-14
	d_p1v = np.sum((v - P[:,Pi])**2)
	for j in xrange(len(edge_ints)):
		d = np.sum((int_coords[:,edge_ints[j]] - P[:,Pi])**2)
		if d > d_p1v + tol:
			forward_ints.append(j)
			distances.append(d)
		
	distances = np.array(distances)
	
	# If so, move to the first forward intersection
	if len(forward_ints) > 0:
		# We just got vnext, and the process of finding it gives you i 
		j = np.argmin(distances)
		i = edge_ints[forward_ints[j]]
		vnext = int_coords[:,i]
		
		# If this intersection is also a vertex, update the edge index
		if np.allclose(P[:,(Pi+1)%n], vnext):
			Pi = (Pi+1)%n
		# only_P is false 
		return vnext, i, Pi, False
	
	# Otherwise, move to the vertex
	else:
		Pi = (Pi + 1) % n
		# Check if it's a vertex
		only_P = True
		for j in xrange(len(P_edges)):
			if np.allclose(P[:,Pi], int_coords[:,j]):
				only_P = False
				i = j
		return P[:,Pi], i, Pi, only_P
		
def get_vnext(P, Q, i, Pi, Qi, v, n, m, num_ints, int_coords, only_P, only_Q, P_edges, Q_edges, vprev, turnleft=True):
	"""
	Get the potential next vertex by determining whether to move along P or along Q and doing so.
	"""

	# If we're on only one polygon but not both
	if only_P:
		vnext, i, Pi, only_P = move_along(P, v, n, i, Pi, only_P, P_edges, int_coords)
	elif only_Q:
		vnext, i, Qi, only_Q = move_along(Q, v, m, i, Qi, only_P, Q_edges, int_coords)
	
	# If we're at an intersection point
	else:

		Pnext = P[:,(Pi+1)%n]
		Qnext = Q[:,(Qi+1)%m]
		
		Pleft = leftright(Pnext, vprev, v)
		Qleft = leftright(Qnext, vprev, v)
		
		left = 'Q'
		if Pleft == 1 and Qleft == 1:
			LR = leftright(Pnext, v, Qnext)
			if LR == 1:
				left = 'P'
			#else:
			#	left = 'Q'

		elif Pleft == 0 and Qleft == 0:
			left = 'P'
				
		#elif Qleft > Pleft:
		#	left = 'Q'

		elif Pleft > Qleft:
			left = 'P'

		else:
			LR = leftright(Pnext, v, Qnext)
			if LR == 1:
				left = 'P'
		#	else:
		#		left = 'Q'

		if left == 'P':
			if turnleft:
				vnext, i, Pi, only_P = move_along(P, v, n, i, Pi, only_P, P_edges, int_coords)
			else:
				vnext, i, Qi, only_Q = move_along(Q, v, m, i, Qi, only_Q, Q_edges, int_coords)

		elif left == 'Q':
			if turnleft:
				vnext, i, Qi, only_Q = move_along(Q, v, m, i, Qi, only_Q, Q_edges, int_coords)
			else:
				vnext, i, Pi, only_P = move_along(P, v, n, i, Pi, only_P, P_edges, int_coords)



	# If we are at an intersection point (only_P and only_Q both false)
	if only_P == False and only_Q == False:
		# Update to make sure Pi and Qi are both right
		Pi = P_edges[i]
		Qi = Q_edges[i]
	
	return vnext, i, Pi, Qi, only_P, only_Q
	
def get_inside_vertices(P, Q):
	"""
	Find the vertices of P that are inside Q and the vertices
	of Q that are inside P.
	"""
	P_inside = []
	Q_inside = []
	
	# Check all the vertices in P to see if they're in Q
	for i in xrange(P.shape[1]):
		if genIO(P[:,i], Q) == 1:
			P_inside.append(i)
	
	# And likewise for Q
	for i in xrange(Q.shape[1]):
		if genIO(Q[:,i], P) == 1:
			Q_inside.append(i)
			
	return np.array(P_inside), np.array(Q_inside)

def polygon_intersection(P, Q):
	"""
	Get the list of polygons that forms the intersection between P and Q.
	
	Return each item in the list as a 2 x n numpy array.
	"""
	
	n = P.shape[1]
	m = Q.shape[1]
	
	# Find the intersection points of the polygon
	int_coords, P_edges, Q_edges = get_intersections(P,Q)
	
	num_ints = len(P_edges)
	
	P_inside, Q_inside = get_inside_vertices(P,Q)
	P_inside_fixed, Q_inside_fixed = np.copy(P_inside), np.copy(Q_inside)

	# If the intersection is empty
	if num_ints == 0:
		# You know you're either disjoint or included
		
		# If the inside vertices of P or Q is all of them
		if len(P_inside) == n:
			return [P]
		elif len(Q_inside) == m:
			return [Q]
		else:		
			return []
			
	pieces = []
	visited = np.zeros(num_ints)
	
	num_iters = 0
	
	# Keep making new pieces until all vertices have been visited.
	while np.any(visited == 0.0) or np.any(P_inside != -1) or np.any(Q_inside != -1) and num_iters < 10:
	
		# Start a new polygon
		if np.any(visited == 0.0):
			i = np.where(visited == 0)[0][0]
			Pi = P_edges[i]
			Qi = Q_edges[i]
			only_P = False
			only_Q = False
			
			start = int_coords[:,i]
			visited[i] = 1
			# Initialize vprev
			vprev = P[:,Pi]
			if np.allclose(vprev, start):
				vprev = P[:,(Pi-1)%n]
				
		elif np.any(P_inside != -1):
			Pi = P_inside[np.where(P_inside != -1)[0][0]]
			P_inside[np.where(P_inside != -1)[0][0]] = -1
			vprev = P[:,(Pi-1)%n]
			start = P[:,Pi]
			only_P = True
			only_Q = False
			Qi = 0
			
		else:
			Qi = Q_inside[np.where(Q_inside != -1)[0][0]]
			Q_inside[np.where(Q_inside != -1)[0][0]] = -1
			vprev = Q[:,(Qi-1)%n]
			start = Q[:,Qi]
			only_Q = True
			only_P = False
			Pi = 0
			
		S = [start]
		v = start
		stop = False
		iters = 0
		while stop == False and iters < 25:
			
			# Turn left (get vnext)
			vnext, i, Pi, Qi, only_P, only_Q = get_vnext(
				P, Q, i, Pi, Qi, v, n, m, num_ints, int_coords, only_P, only_Q, 
				P_edges, Q_edges, vprev)
			
			if iters == 0:
				s1 = vnext
			
			for j in xrange(len(P_inside)):
				if np.allclose(P[:,P_inside[j]], vnext):
					P_inside[j] = -1
			for j in xrange(len(Q_inside)):
				if np.allclose(Q[:,Q_inside[j]], vnext):
					Q_inside[j] = -1

			
			allowed = True
			if only_P or only_Q:
				allowed = False 

				for j in xrange(len(P_inside)):
					if Pi == P_inside_fixed[j]:
						allowed = True
				for j in xrange(len(Q_inside)):
					if Qi == Q_inside_fixed[j]:
						allowed = True
				

			# Is this the same edge we started with?
			if iters > 0 and same_edge(start, s1, v, vnext):
				stop = True
				S = S[:-1]
			
			elif allowed == False:
				stop = True

			# Otherwise, add the new point to the polygon 
			else:
				S.append(vnext)
				vprev = v
				v = vnext
				# Mark v as visited
				visited[i] = 1
			
			iters += 1
		
		if np.allclose(S,S[-1]) and len(S) > 1:
			S = S[:-1]
		# Add the polygon to your pieces
		if len(S) > 0:
			pieces.append(np.array(S).T)

		num_iters += 1
	
	return pieces
	
def polygon_union(P, Q):	
	"""
	Get the list of polygons that forms the union of P and Q.
	
	Return each item in the list as a 2 x n numpy array.
	"""
	n = P.shape[1]
	m = Q.shape[1]
	
	int_coords, P_edges, Q_edges = get_intersections(P,Q)
	num_ints = len(P_edges)
	visited = np.zeros(num_ints)
	P_inside, Q_inside = get_inside_vertices(P,Q)
	
	# If there are no intersection points
	if num_ints == 0:
		# Either they're disjoint or one is contained in the other
		if len(P_inside) == n:
			return [Q], []
		elif len(Q_inside) == m:
			return [P], []
		else:		
			return [P,Q], []
	
	
	P_tovisit = np.setdiff1d(range(n), P_inside)
	Q_tovisit = np.setdiff1d(range(m), Q_inside)
	
	Pmin = np.argmin(P[0,:])
	Qmin = np.argmin(Q[0,:])
	i = 0
	if P[0,Pmin] <= Q[0,Qmin]:
		Pi = Pmin
		start = P[:,Pi]
		vprev = P[:,(Pi-1)%n]
		Qi = 0
		only_P = True
		only_Q = False
		for j in xrange(len(P_tovisit)):
			if Pi == P_tovisit[j]:
				P_tovisit[j] = -1
		for j in xrange(num_ints):
			if np.allclose(start, int_coords[:,j]):
				i = j
				Qi = Q_edges[j]
				only_P = False
	else:
		Qi = Qmin
		start = Q[:,Qi]
		vprev = Q[:,(Qi-1)%m]
		Pi = 0
		only_Q = True
		only_P = False
		for j in xrange(len(Q_tovisit)):
			if Qi == Q_tovisit[j]:
				Q_tovisit[j] = -1
		for j in xrange(num_ints):
			if np.allclose(start, int_coords[:,j]):
				i = j
				Pi = P_edges[j]
				only_Q = False
	
	S = [start]
	v = start
	stop = False
	iters = 0
	
	# Get the boundary
	while stop == False:
			
		# Turn right (get vnext)
		vnext, i, Pi, Qi, only_P, only_Q = get_vnext(
			P, Q, i, Pi, Qi, v, n, m, num_ints, int_coords, only_P, only_Q, 
			P_edges, Q_edges, vprev, turnleft=False)
		
		# Check your intersection points
		for j in xrange(num_ints):
			if np.allclose(vnext, int_coords[:,j]):
				visited[j] = 1
				Pi = P_edges[j]
				Qi = Q_edges[j]
		
		# Check your vertices to see if you visited any
		for j in xrange(len(P_tovisit)):
			if P_tovisit[j] != -1:
				if np.allclose(vnext, P[:,P_tovisit[j]]):
					P_tovisit[j] = -1
					
		for j in xrange(len(Q_tovisit)):
			if Q_tovisit[j] != -1:
				if np.allclose(vnext, Q[:,Q_tovisit[j]]):
					Q_tovisit[j] = -1
		
		if iters == 0:
			s1 = vnext
			
		if iters > 0 and same_edge(start, s1, v, vnext):
			stop = True
			S = S[:-1]
			boundary = [np.array(S).T]
			
		else:
			S.append(vnext)
			vprev = v
			v = vnext
			
		iters += 1
		
	# Get the holes
	holes = []
	allowed_pts = []
	for i in xrange(len(visited)):
		if visited[i] == 0.0:
			allowed_pts.append(int_coords[:,i])
	for i in xrange(len(P_tovisit)):
		if P_tovisit[i] != -1:
			allowed_pts.append(P[:,P_tovisit[i]])
	for i in xrange(len(Q_tovisit)):
		if Q_tovisit[i] != -1:
			allowed_pts.append(Q[:,Q_tovisit[i]])
	
	num_iters = 0
	while np.any(visited == 0.0) or np.any(P_tovisit != -1) or np.any(Q_tovisit != -1):
			
		only_P = False
		only_Q = False
		
		# Start a new piece
		if np.any(visited == 0.0):
			i = np.where(visited == 0)[0][0]
			Pi = P_edges[i]
			Qi = Q_edges[i]
			start = int_coords[:,i]
			visited[i] = 1
			vprev = P[:,Pi]
			if np.allclose(vprev, start):
				vprev = P[:,(Pi-1)%n]
				
		elif np.any(P_tovisit != -1):
			Pi = P_tovisit[np.where(P_tovisit != -1)[0][0]]
			P_tovisit[np.where(P_tovisit != -1)[0][0]] = -1
			vprev = P[:,(Pi-1)%n]
			start = P[:,Pi]
			only_P = True
			Qi = 0
		else:
			Qi = Q_tovisit[np.where(Q_tovisit != -1)[0][0]]
			Q_tovisit[np.where(Q_tovisit != -1)[0][0]] = -1
			vprev = Q[:,(Qi-1)%n]
			start = Q[:,Qi]
			only_Q = True
			Pi = 0
		
		S = [start]
		v = start
		stop = False
		iters = 0
		
		# Turn right until the place you're trying to go has already been visited.
		while stop == False and iters < 25:
				
			# Turn right (get vnext)
			vnext, i, Pi, Qi, only_P, only_Q = get_vnext(
				P, Q, i, Pi, Qi, v, n, m, num_ints, int_coords, only_P, only_Q, 
				P_edges, Q_edges, vprev, turnleft=False)
			
			allowed = False
			for j in xrange(len(allowed_pts)):
				if np.allclose(allowed_pts[j], vnext):
					allowed = True
			if allowed == False:
				j = 0
				while j < len(allowed_pts):
					if np.allclose(allowed_pts[j], v):
						del allowed_pts[j]
					else:
						j += 1
				break
			
			# Check your intersection points
			for j in xrange(num_ints):
				if np.allclose(vnext, int_coords[:,j]):
					visited[j] = 1
					Pi = P_edges[j]
					Qi = Q_edges[j]
			
			# Check your vertices to see if you visited any
			for j in xrange(len(P_tovisit)):
				if P_tovisit[j] != -1:
					if np.allclose(vnext, P[:,P_tovisit[j]]):
						P_tovisit[j] = -1
			for j in xrange(len(Q_tovisit)):
				if Q_tovisit[j] != -1:
					if np.allclose(vnext, Q[:,Q_tovisit[j]]):
						Q_tovisit[j] = -1
			
			if iters == 0:
				s1 = vnext
				
			if iters > 0 and same_edge(start, s1, v, vnext):
				stop = True
				S = S[:-1]
				holes.append(np.array(S).T)
			else:
				S.append(vnext)
				vprev = v
				v = vnext
			
			iters += 1
		num_iters += 1
	
	return boundary, holes	
