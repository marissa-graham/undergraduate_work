from intersect2d import *

# Convenience functions for displaying the results of polygon_intersection and polygon_union
def display_intersection(P, Q, Pname, Qname):
	plt.scatter(P[0,:], P[1,:],c='C0')
	plt.plot(P[0,:], P[1,:],c='C0')
	plt.plot([P[0,-1],P[0,0]],[P[1,-1],P[1,0]],c='C0')
	
	plt.scatter(Q[0,:], Q[1,:],c='C1')
	plt.plot(Q[0,:], Q[1,:],c='C1')
	plt.plot([Q[0,-1],Q[0,0]],[Q[1,-1],Q[1,0]],c='C1')

	plt.xticks([])
	plt.yticks([])
	
	polygons = polygon_intersection(P, Q)
	for l in xrange(len(polygons)):
		polygon = polygons[l]
		plt.plot(polygon[0,:], polygon[1,:], c='green', linewidth=5)
		plt.plot([polygon[0,-1],polygon[0,0]],[polygon[1,-1],polygon[1,0]],c='green', linewidth=5)
		
	plt.title(Pname + ' + ' + Qname + ' intersection')
	plt.show()
	
def display_union(P, Q, Pname, Qname):
	plt.scatter(P[0,:], P[1,:],c='C0')
	plt.plot(P[0,:], P[1,:],c='C0')
	plt.plot([P[0,-1],P[0,0]],[P[1,-1],P[1,0]],c='C0')
	
	plt.scatter(Q[0,:], Q[1,:],c='C1')
	plt.plot(Q[0,:], Q[1,:],c='C1')
	plt.plot([Q[0,-1],Q[0,0]],[Q[1,-1],Q[1,0]],c='C1')

	plt.xticks([])
	plt.yticks([])
	
	boundary, holes = polygon_union(P, Q)
	for l in xrange(len(boundary)):
		polygon = boundary[l]
		plt.plot(polygon[0,:], polygon[1,:], c='green', linewidth=5)
		plt.plot([polygon[0,-1],polygon[0,0]],[polygon[1,-1],polygon[1,0]],c='green', linewidth=5)
		
	for l in xrange(len(holes)):
		polygon = holes[l]
		plt.plot(polygon[0,:], polygon[1,:], c='purple', linewidth=5)
		plt.plot([polygon[0,-1],polygon[0,0]],[polygon[1,-1],polygon[1,0]],c='purple', linewidth=5)
		
	plt.title(Pname + ' + ' + Qname + ' union')
	plt.show()
				
square = np.array([[0.2,0.8,0.8,0.2],
				   [0.2,0.2,0.8,0.8]])

bigsquare = np.array([[0.0,1,1,0],[0,0,1,1]])
				   
quad = np.array([[0.1,0.9,0.5,0.5],
				 [0.1,0.5,0.5,0.9]])

trek = np.array([[0.1,0.5,0.9,0.9,0.5,0.1],
				 [0.1,0.3,0.1,0.5,0.8,0.5]])

star = np.array([[0.1,0.3,0.2,0.5,0.8,0.7,0.9,0.65,0.5,0.35],
				 [0.7,0.4,0.1,0.3,0.1,0.4,0.7,0.65,0.9,0.65]])

triangle = np.array([[0.1,0.8,0.5],
					 [0.1,0.5,0.8]])
triangle2 = np.array([[0.1,0.9,0.5],
					  [0.1,0.5,0.9]])

hexagon = np.array([[0.2,0.1,0.4,0.8,0.8,0.6],
					[0.75,0.3,0.05,0.35,0.9,0.95]])

weirdie = np.array([[0.1,0.2,0.8,0.9,0.8,0.6,0.8,0.75,0.5,0.25,0.2,0.4,0.2],
					[0.7,0.1,0.1,0.7,0.9,0.8,0.7,0.4,0.2,0.4,0.7,0.8,0.9]])

shapes = [triangle, bigsquare, square, quad, trek, star, hexagon, weirdie]
shapenames = ['triangle', 'bigsquare', 'square', 'quad', 'trek', 'star', 'hexagon', 'weirdie']

wprime = np.array([[0.2,0.4,0.2,0.25,0.5,0.75,0.8,0.6,0.8,0.9,0.8,0.2,0.1],
				   [0.1,0.2,0.3,0.6,0.8,0.6,0.3,0.2,0.1,0.3,0.9,0.9,0.3]])

tri = np.array([[0,0.5,0.5,0.3,0.5,0.5],[0.5,0,0.3,0.5,0.7,1]])
backtri = np.array([[0.5,1,0.5,0.5,0.7,0.5],[0,0.5,1,0.7,0.5,0.3]])

display_union(hexagon, bigsquare, 'hexagon', 'bigsquare')

# Bunch of pictures for slides
if False:
	pass
	"""
	P = tri
	Q = np.copy(backtri)

	plt.scatter(P[0,:], P[1,:],c='C0')
	plt.plot(P[0,:], P[1,:],c='C0')
	plt.plot([P[0,-1],P[0,0]],[P[1,-1],P[1,0]],c='C0')

	plt.scatter(Q[0,:], Q[1,:],c='C1')
	plt.plot(Q[0,:], Q[1,:],c='C1')
	plt.plot([Q[0,-1],Q[0,0]],[Q[1,-1],Q[1,0]],c='C1')

	plt.xticks([])
	plt.yticks([])

	polygons = polygon_intersection(P, Q)
	for l in xrange(len(polygons)):
		polygon = polygons[l]
		print polygon.shape
		#if polygon.shape[1] == 1:
		#	plt.scatter(polygon[0], polygon[1], c='green', marker='*', s=150)
		plt.plot(polygon[0,:], polygon[1,:], c='green', linewidth=5)
		plt.plot([polygon[0,-1],polygon[0,0]],[polygon[1,-1],polygon[1,0]],c='green', linewidth=5)
		
	plt.title('Edge-Only Intersection')
	plt.show()

	P = np.array([[0,0.5,0],[0,0.5,1]])
	Q = np.array([[0.5,1,1],[0.5,0,1]])

	plt.scatter(P[0,:], P[1,:],c='C0')
	plt.plot(P[0,:], P[1,:],c='C0')
	plt.plot([P[0,-1],P[0,0]],[P[1,-1],P[1,0]],c='C0')

	plt.scatter(Q[0,:], Q[1,:],c='C1')
	plt.plot(Q[0,:], Q[1,:],c='C1')
	plt.plot([Q[0,-1],Q[0,0]],[Q[1,-1],Q[1,0]],c='C1')

	plt.xticks([])
	plt.yticks([])

	polygons = polygon_intersection(P, Q)
	for l in xrange(len(polygons)):
		polygon = polygons[l]
		print polygon.shape
		if polygon.shape[1] == 1:
			plt.scatter(polygon[0], polygon[1], c='green', marker='*', s=150)
		plt.plot(polygon[0,:], polygon[1,:], c='green', linewidth=5)
		plt.plot([polygon[0,-1],polygon[0,0]],[polygon[1,-1],polygon[1,0]],c='green', linewidth=5)
		
	plt.title('Single Point of Intersection')
	plt.show()

	P = np.copy(bigsquare)
	Q = np.copy(bigsquare)
	Q[0,:] = Q[0,:] + 1 

	display_intersection(P,Q,'p','q')

	plt.xticks([])
	plt.yticks([])

	int_coords, P_edges, Q_edges = get_intersections(P,Q)
	#plt.scatter(int_coords[0,:], int_coords[1,:], c='green', marker='*',s=125,zorder=10)
	#plt.scatter(P[0,:], P[1,:],c='C0')
	plt.plot(P[0,:], P[1,:],c='C0')
	plt.plot([P[0,-1],P[0,0]],[P[1,-1],P[1,0]],c='C0')

	#plt.scatter(Q[0,:], Q[1,:],c='C1')
	plt.plot(Q[0,:], Q[1,:],c='C1')
	plt.plot([Q[0,-1],Q[0,0]],[Q[1,-1],Q[1,0]],c='C1')
	plt.title('Disjoint Polygons')
	plt.show()


	P = star
	Q = bigsquare

	plt.scatter(P[0,:], P[1,:],c='C0')
	plt.plot(P[0,:], P[1,:],c='C0')
	plt.plot([P[0,-1],P[0,0]],[P[1,-1],P[1,0]],c='C0')

	plt.scatter(Q[0,:], Q[1,:],c='C1')
	plt.plot(Q[0,:], Q[1,:],c='C1')
	plt.plot([Q[0,-1],Q[0,0]],[Q[1,-1],Q[1,0]],c='C1')

	plt.xticks([])
	plt.yticks([])

	polygons = polygon_intersection(P, Q)
	for l in xrange(len(polygons)):
		polygon = polygons[l]
		plt.plot(polygon[0,:], polygon[1,:], c='green', linewidth=5)
		plt.plot([polygon[0,-1],polygon[0,0]],[polygon[1,-1],polygon[1,0]],c='green', linewidth=5)
		
	plt.title('Intersection with an Interior Polygon')
	plt.show()
		

	display_intersection(trek, triangle, 'trek', 'triangle')
	display_intersection(triangle, trek, 'triangle', 'trek')
	if False:
		for i in xrange(len(shapes)):
			display_union(weirdie, shapes[i], 'weirdie', shapenames[i])
	"""

if True:
	for i in xrange(len(shapes)):
		for j in xrange(i,len(shapes)):
			print '\n\n', shapenames[i], shapenames[j]
			#offset_shape = np.copy(shapes[j])
			#offset_shape[0,:] = offset_shape[0,:] + 1
			#display_intersection(shapes[i], offset_shape, shapenames[i], 'offset '+shapenames[j])
			display_intersection(shapes[i], shapes[j], shapenames[i], shapenames[j])
			display_union(shapes[i], shapes[j], shapenames[i], shapenames[j])