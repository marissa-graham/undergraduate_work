from __future__ import division
import sys
import numpy as np
from matplotlib import pyplot as plt

def specs():
	"""
	Read in a .obj file of a closed polyhedron and output a .obj file that contains
	a polyhedron produced by performing a single Catmull-Clark refinement on the original
	polyhedron. Vertices will be unique and the files will contain only vertices and faces.

	1. Produce four successive refinements of cube0.obj (cube1.obj, cube2.obj, cube3.obj, cube4.obj)
	2. Produce four successive refinements of tet0.obj (tet1.obj, tet2.obj, tet3.obj, tet4.obj)
	"""

class Vertex(object):
	def __init__(self, coords):
		"""
		x: Cartesian coordinates (x,y,z) of the vertex, as read in from the .obj file.
		valence: Number of faces/edges that meet at the vertex.
		faces: Vector of integers containing the face vector indices of all faces that share that vertex.
		edges: Vector of integers containing the edge vector indices of all edges that connect to the vertex.
		indices: Cartesian coordinates of the computed vertex point corresponding to that vertex.
		i0: Integer indicating where the vertex point occurs in the output file.
		"""
		self.x = np.array(coords,dtype=np.float)
		self.valence = None
		self.faces = None
		self.edges = None
		self.vert_pt = None
		self.i0 = None

class Edge(object):
	def __init__(self, i1, i2):
		"""
		i1, i2: Indices of the vertices of the endpoints
		f1, f2: Indices of the faces that meet at the edge
		vert_pt: Cartesian coordinates of the computed edge point.
		i0: Integer indicating where edge point is listed as a vertex in the output file.
		"""
		self.i1, self.i2 = i1, i2
		self.f1, self.f2 = None, None
		self.vert_pt = None
		self.i0 = None

class Face(object):
	def __init__(self, vertices, edges):
		"""
		n: An integer indicating the number of vertices/edges on the face.
		verts: A vector of integers of length n indicating which vertices lie on the face, as read
			in from the .obj file.
		edges: A vector of integers of length n indicating which edges outline the face, such that
			the neighboring edges in the vector are neighboring edges on the face.
		f_vert: Cartesian coordinates of the computed face point.
		i0: Integer indicating where the face point is listed as a vertex in the outputted .obj file.
		"""
		self.n = len(vertices)
		self.verts = np.array(vertices)
		self.edges = np.array(edges)
		self.f_vert = None
		self.i0 = None

def main():
    # Ensure command line input is correct
    if len(sys.argv) < 3:
            print "Syntax:\n python cs557proj5.py testfile outputfile"
            return

    # Read in input file and store vertices in V and faces in F (x for vertices, n, verts for faces)
    V = [0] # need 1-indexing, so dummy zero
    F = [] # TODO: need 1-indexing?
    E = []
    edge_index = 0
    line = "dummy"

    with open(sys.argv[1]) as f:
    	while line:
			line = f.readline().split()
			if line:
				if line[0] == 'v':
					vert = Vertex([float(line[i]) for i in xrange(1,len(line))])
					V.append(vert)

				elif line[0] == 'f':
					edges = []

					for i in xrange(1,len(line)):
						i1, i2 = int(line[i]),int(line[i%(len(line)-1)+1])
						new_edge = Edge(min(i1,i2),max(i1,i2))
						unique = True
						for j in xrange(len(E)):
							if E[j].i1 == new_edge.i1 and E[j].i2 == new_edge.i2:
								unique = False
								edges.append(j)
						#print "edge index is currently ", edge_index
						if unique:
							#print new_edge.i1, new_edge.i2
							E.append(new_edge)
							edges.append(edge_index)
							edge_index += 1

					newface = Face([int(line[i]) for i in xrange(1,len(line))],edges)
					F.append(newface)

	#print len(F)
	#print len(E)

	# Go through and figure out which faces meet at each edge
	for i in xrange(len(E)):
		faces = []
		for j in xrange(len(F)):
			# Check if the edge is associated with the face
			for k in xrange(len(F[j].edges)):
				#print "testing pair"
				#print E[i].i1, E[i].i2
				#print E[F[j].edges[k]].i1, E[F[j].edges[k]].i2
				if E[i].i1 == E[F[j].edges[k]].i1 and E[i].i2 == E[F[j].edges[k]].i2:
					#print "match found"
					#print E[i].i1, E[i].i2
					faces.append(j)

		#print faces
		E[i].f1, E[i].f2 = faces[0], faces[1]

	# Calculate valence, faces, and edges for each vertex
	for i in xrange(1,len(V)):
		faces = []
		edges = []
		for j in xrange(len(F)):
			# Check if the vertex is associated with the face
			if i in F[j].verts:
				faces.append(j)

		for j in xrange(len(E)):
			# Check if the vertex is associated with the edge
			if i == E[j].i1 or i == E[j].i2:
				edges.append(j)

		if len(faces) != len(edges):
			print "Error: Edge valence should match face valence."
			break
		valence = len(faces)

		V[i].edges = edges
		V[i].valence = valence
		V[i].faces = faces
					
    # Compute and store the face point Cartesian coordinates for each face
    for face in F:
    	tot = 0
    	for i in xrange(face.n):
    		tot += V[face.verts[i]].x
    	tot /= face.n
    	face.f_vert = tot

    # Compute and store the edge point Cartesian coordinates for each edge
    for edge in E:
    	edge.vert_pt = 0.25*(V[edge.i1].x + V[edge.i2].x + F[edge.f1].f_vert + F[edge.f2].f_vert)

    # Compute and store the vertex point Cartesian coordinates for each vertex
    for i in xrange(1,len(V)):
    	v = V[i]
    	n = v.valence
    	tot = (n-2)/n*v.x
    	for j in xrange(n):
    		tot += (F[v.faces[j]].f_vert + E[v.edges[j]].vert_pt)/n**2
    	V[i].vert_pt = tot
    for i in xrange(1,len(V)):
    	#print i
    	#v = V[i]
    	#n = v.valence
    	#tot = (n-2)/n*v.x
    	#for j in xrange(n):
    	#	tot += (F[v.faces[j]].f_vert + E[v.edges[j]].vert_pt)/n**2
    	#V[i].vert_pt = tot
    	#print len(V)
    	#print tot
    	#print V[i].vert_pt

	k = 1
	with open(sys.argv[2], 'w') as f:
		# Output to the new .obj file the Cartesian coordinates for each computed face point, edge point,
		# and vertex point. While doing so, assign the appropriate value to the variable io for each face
		# point, edge point, and vertex point.
		for face in F:
			#print k, "v {} {} {}".format(*face.f_vert)
			f.write("v {} {} {}\n".format(*face.f_vert))
			face.i0 = k
			k += 1
		for edge in E:
			#print k, "v {} {} {}".format(*edge.vert_pt)
			f.write("v {} {} {}\n".format(*edge.vert_pt))
			edge.i0 = k
			k += 1 
		for j in xrange(1,len(V)):
			vertex = V[j]
			#print vertex.vert_pt
			#print k, "v {} {} {}".format(*vertex.vert_pt)
			f.write("v {} {} {}\n".format(*vertex.vert_pt))
			V[j].i0 = k
			k += 1 

		# Output the faces to the new .obj file. Do this by visiting each face object in vector F. 
		# Each original face will be split into n four-sided faces in the output file. The four integers 
		# to be output for each face will be the io number for a vertex point, followed by the io number 
		# for an edge point, then the io number for the face point, then the io number for another 
		# edge point.

		# For each face
		for face in F:
			# For each vertex
			for i in xrange(face.n):
				# Get edges associated with both face and vertex
				c_edges = np.intersect1d(face.edges,V[face.verts[i]].edges)

				next_vert = face.verts[(i+1)%face.n]
				if next_vert == E[c_edges[0]].i1 or next_vert == E[c_edges[0]].i2:
					right_edge = c_edges[0]
					left_edge = c_edges[1]
				else:
					right_edge = c_edges[1]
					left_edge = c_edges[0]

				v1 = V[face.verts[i]].i0
				v2 = E[right_edge].i0
				v3 = face.i0
				v4 = E[left_edge].i0
				f.write("f {} {} {} {}\n".format(v1,v2,v3,v4))
				#print face.verts[i],',', face.verts[(i+1)%face.n], ':', E[common_edges[0]].i1, E[common_edges[0]].i2
				#print '  ', E[common_edges[1]].i1, E[common_edges[1]].i2

				# Figure out which is right-hand edge, aka the edge whose other vertex is the index
					# after the current vertex

				# Then output:
					# vertex pt
					# edge point associated with right-hand edge
					# face point
					# other edge point
			

			#for i in xrange(1,face.n+1):
			#	print "Edges associated with vertex: ", V[i].edges
			#	print len(E)
			#	for j in xrange(len(V[i].edges)):
			#		temp = V[i].edges
			#		#print temp[j]
			#		print "Vertices associated with edge", j, ":",E[temp[j]].i1, E[temp[j]].i2


main()