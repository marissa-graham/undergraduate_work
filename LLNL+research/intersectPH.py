from Polyhedron import *

def tag_LR(P, Q):
	"""
	Determine whether the vertices of the polyhedron P are left, right, or on
	the plane defined by the first three points of Q. Then go through the faces 
	and determine whether they are left, right, or intersecting the plane.
	"""
	n = np.cross(Q[:,1] - Q[:,0], Q[:,2] - Q[:,1])
	n = n/la.norm(n)

	pre_LR = np.dot((Q[:,0].reshape((3,1)) - P.Vcoords).T, n.reshape((3,1)))
	LR = np.zeros_like(pre_LR)
	LR[pre_LR < 0] = 1
	LR[pre_LR > 0] = -1
	LR = LR.reshape(LR.shape[0],)

	faceLR = np.zeros(P.nf)
	for i in xrange(P.nf):
		# Does it still count if there's only one non-L/R vertex?
		fLR = LR[P.F[i].verts]
		if np.all(fLR == 1):
			faceLR[i] = 1
		elif np.all(fLR == -1):
			faceLR[i] = -1
		else:
			faceLR[i] = 0

	Pleft = Polyhedron()
	Pleft.filename = 'Pleft.obj'
	Pcenter = Polyhedron()
	Pcenter.filename = 'Pcenter.obj'
	vcoordsC = []
	vcoordsL = []
	for i in xrange(P.nv):
		if LR[i] < 0:
			Pleft.V.append(P.V[i])
			vcoordsL.append(P.V[i].x)
			Pleft.nv += 1
		elif LR[i] == 0:
			Pcenter.V.append(P.V[i])
			vcoordsC.append(P.V[i].x)
			Pcenter.nv += 1

	Pleft.Vcoords = np.array(vcoordsL).T
	Pcenter.Vcoords = np.array(vcoordsC).T

	for i in xrange(P.nf):
		if faceLR[i] < 0:
			Pleft.F.append(P.F[i])
			Pleft.nf += 1
		elif faceLR[i] == 0:
			Pcenter.F.append(P.F[i])
			Pcenter.nf += 1

	#Pleft.print_poly()
	Pleft.display()
	Pleft.write('Pleft.obj')

