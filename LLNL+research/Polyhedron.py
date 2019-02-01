import mpl_toolkits.mplot3d as a3
from scipy import linalg as la
import pylab
import numpy as np
import intersect2d
import intersect3d

class Vertex(object):
	def __init__(self, coords):
		"""
		x: Cartesian coordinates of the vertex.
		faces: indices of the faces the vertex is part of
		faceinds: index of the vertex with respect to the corresponding face
			in faceinds
		"""
		self.x = np.array(coords)
		self.faces = []
		self.faceinds = []
		
class Edge(object):
	def __init__(self, p1, p2):
		"""
		p1, p2: Indices of the vertices that define the edge.
		f1, f2: Indices of the faces that meet at the edge.
		"""
		self.p1, self.p2 = p1, p2
		self.f1, self.f2 = None, None
		# Edges should know their edge number on each face

# ?????
def Intersection(object):
	def __init__(self):
		# coordinates
		# edge of occurrence
		# face(s) of occurence?
		pass
		
class Face(object):
	def __init__(self, vertices):
		"""
		n:     Number of vertices/edges on the face.
		verts: Vector of integers containing the indices of the
			   vertices that define the face (counterclockwise order).
		coords: Array with the coordinates of the vertices.
		edges: Vector of integers containing the indices of the 
			edges that define the face, such that neighboring
			indices in the vector correspond to neighboring edges
			on the face.
		"""
		self.n = len(vertices)
		self.verts = np.array(vertices)
		self.coords = []
		self.edges = None

	def get_normal(self, pts):
		self.coords = pts
		n = np.cross(pts[:,1] - pts[:,0], pts[:,2] - pts[:,1])
		self.normal = n/la.norm(n)
		self.M, self.Minv = intersect3d.to_2d(pts)
		
		P = np.dot(self.M, pts)
		if intersect2d.area2d(P[:-1,:]) < 0:
			P = P[:,::-1]
		self.d = P[-1,0]
		self.coords2d = P[:-1,:]
 
	def area(self):
		""" Compute the area of the face. """
		self.A = intersect2d.area2d(self.coords2d)
		#print self.A
		return self.A

	def on_face(self, x):
		""" 
		Determine whether the point x is on the face or not.
		"""
		p = np.dot(M, x)[:-1]
		return intersect2d.genIO(p, self.coords2d)
		
	def face_distance(self, x):
		""" 
		Compute the distance from the face to a given point pt.
		"""
		return np.dot(self.normal, self.coords[:,0] - x)
	
class Path(object):
	"""
	A path on the polyhedron.
	
	Feed it the start and end point and the vertices visited,
	in order.
	"""
	def __init__(self, P, p1, p2, vertices):
		self.P = P
		self.p1 = p1
		self.p2 = p2
		self.verts = vertices
		
	def on_path(self, x):
		"""
		Determine whether the point x is on the path.
		"""
		pass
	
	def __len__(self):
		"""
		Determine the length of the path.
		"""
		pass

class Polyhedron(object):
	"""
	It doesn't matter what order the vertices, edges, or faces 
	are in. If line information is not provided, the order of the
	edges is the order in which they are encountered while reading
	the faces.
	
	Edges are bidirectional.
	
	The order of the edges in each face DOES matter and is assumed
	to be counterclockwise.
	"""
	
	def __init__(self):
		"""
		V: List of Vertices
		F: List of Faces
		E: List of Edges

		Vcoords: numpy array of all the coordinates of the vertices
		(for bounding boxes, choosing planes to intersect with, etc.)
		"""
		self.nv = 0
		self.ne = 0
		self.nf = 0
		
		self.V = [] 
		self.Vcoords = []
		self.F = []
		self.E = []

	def read(self, filename):
		
		self.filename = filename
		
		surface_normals = False

		print "reading in file"
		
		for readLine in open(filename):
		
			line = readLine.split()
			if line:	
				if line[0] == 'v':
					vert = Vertex([float(line[i]) for i in xrange(1,len(line))])
					self.V.append(vert)
					self.Vcoords.append(vert.x)
					self.nv += 1
				
				elif line[0] == 'l':
					p1, p2 = int(line[1]) - 1, int(line[2]) - 1
					self.E.append(Edge(min(p1,p2),max(p1,p2)))
					self.ne += 1

				elif line[0] == 'vn':
					surface_normals = True
				
				elif line[0] == 'f':
					if surface_normals:
						n = len(line)
						l = line[i]
						face = Face([int(l.rpartition('/')[0][:-1]) - 1 for i in xrange(1,n)])
					else:
						face = Face([int(line[i]) - 1 for i in xrange(1,len(line))])
					
					vcoords = []
					for i in xrange(face.n):
						self.V[face.verts[i]].faces.append(self.nf)
						self.V[face.verts[i]].faceinds.append(i)
						vcoords.append(self.V[face.verts[i]].x)
					face.get_normal(np.array(vcoords).T)
					self.F.append(face)
					self.nf += 1

		self.Vcoords = np.array(self.Vcoords).T
		print "number of faces:", self.nf	
		print "number of vertices:", self.nv		
		
		edges = set()
		edge_dict = dict()
		for i in xrange(self.nf):

			# Go through the edges in the face and tag the appropriate edges
			for j in xrange(self.F[i].n):
				
				# Get the edge
				v1, v2 = self.F[i].verts[j], self.F[i].verts[(j+1)%self.F[i].n]
				minv = min(v1, v2)
				maxv = max(v1, v2)
				tup = (minv, maxv)

				# Tag the edge with the face number (f1 if f1 == None, else f2)
				unique = True
				if tup in edges:
					unique = False
					k = edge_dict[(minv, maxv)]
					if self.E[k].f1 == None:
						self.E[k].f1 = i
					elif self.E[k].f2 == None:
						self.E[k].f2 = i
					else:
						print "warning: too many faces per edge"
		
				# If the edge isn't in the list, we need to add it
				if unique:
					edges.add(tup)
					edge_dict[tup] = self.ne
					new_edge = Edge(minv, maxv)
					new_edge.f1 = i
					self.E.append(new_edge)
					self.ne += 1

		print "number of edges:", self.ne
		print "check for topological consistency"
		consistent = True
		for e in xrange(self.nf):
			if self.E[k].f2 == None:
				consistent = False
			elif self.E[k].f1 == self.E[k].f2:
				consistent = False
		if consistent == False:
			raise ValueError('Polyhedron is not topologically consistent')
		else:
			print "polyhedron is topologically consistent"

	def print_poly(self):
		"""
		Print the vertices, edges, and faces of the polyhedron.
		"""
		print "\n", self.filename, '\n'
		
		print "Vertices:"
		for i in xrange(self.nv):
			v = self.V[i]
			print i, v.x
			print v.faces, v.faceinds

		print "Vcoords array:\n", self.Vcoords
			
		print "\nEdges:"
		for i in xrange(self.ne):
			e = self.E[i]
			print 'Edge', i, ', defined by vertices', e.p1, 'and', e.p2, ', is attached to faces', e.f1, 'and', e.f2
			
		print "\nFaces:"
		for i in xrange(self.nf):
			face = self.F[i]
			print i, face.verts
			
	def write(self, fname):
		"""
		Write the polyhedron to an .obj file called fname.
		"""
		with open(fname, 'w') as f:
			for i in xrange(self.nv):
				f.write('v {} {} {}\n'.format(*self.V[i].x))
			for i in xrange(self.ne):
				f.write('l {} {}\n'.format(self.E[i].p1 + 1, self.E[i].p2 + 1))
			for i in xrange(self.nf):
				fstring = 'f '
				for j in xrange(self.F[i].n):
					fstring += str(self.F[i].verts[j] + 1) + ' '
				fstring += '\n'
				f.write(fstring)
		
	def display(self, labels=True):
		"""
		Display the polyhedron in 3d, with vertices, edges
		and faces labeled with their indices.
		
		Display face numbers in the middle of the face
		Display edge numbers in the middle of the edge
		"""
		ax = a3.Axes3D(pylab.figure())
			
		#for i in xrange(self.nv):
		#	v = self.V[i].x
		#	ax.scatter(v[0], v[1], v[2], c='C0',s=1)
		#	if labels:
		#		ax.text(v[0], v[1], v[2], str(i), fontsize=15)
		ax.scatter(self.Vcoords[0,:], self.Vcoords[1,:], self.Vcoords[2,:],c='C0',s=1)
			
		#for i in xrange(self.ne):
		#	x1, x2 = self.V[self.E[i].p1].x, self.V[self.E[i].p2].x
		#	p = 0.5*(x1 + x2)
		#	ax.plot([x1[0],x2[0]], [x1[1],x2[1]], [x1[2],x2[2]], c='C2')
		#	if labels:
		#		ax.text(p[0], p[1], p[2], str(i), fontsize=15, color='green')
			
		for i in xrange(self.nf):
			#print self.F[i].area()
			poly = a3.art3d.Poly3DCollection([self.F[i].coords.T])
			#p = np.mean(self.F[i].coords,axis=1)
			#p = self.F[i].get_centroid()
			#if labels:
				#ax.text(p[0], p[1], p[2], str(i), fontsize=15, color='red')
			poly.set_edgecolor('k')
			poly.set_facecolor([0.5, 0.5, 1.0, 0.4])
			ax.add_collection(poly)
		
		pylab.show()
	
	def volume(self):
		"""
		Calculate the volume of the polyhedron.
		"""
		p = 0
		for i in xrange(self.nv):
			p += self.V[i].x
		p /= self.nv
		
		vol = 0
		for i in xrange(self.nf):
			A = self.F[i].area()
			h = self.F[i].face_distance(p)
			
			vol += A*h/3

		self.volume = vol
		return vol
	




