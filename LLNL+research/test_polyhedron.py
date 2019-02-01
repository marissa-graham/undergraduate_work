# coding: utf-8

from Polyhedron import *
from intersectPH import *
"""
Things you need to do:
 
- When you read in the file, be sure to tag the vertex and edge connectivity (minimal slowdown)
- The vertices need to know which faces they’re attached to. You CANNOT have a vertex in the middle 
  of an edge. That would imply more than two faces per edge, which violates your consistency condition.
- Make sure the polyhedron has a numpy array of the vertices. See if there’s a way you can pass the 
  ENTIRE L/R operation for all the vertices to numpy to speed it up, since it’s just a subtraction and a
  dot product.
- Face LR function:
	* Return 1 if all ones
	* Return -1 if all -1
	* Return 0 otherwise
- For now, make your slicing function return a polyhedron consisting of all of the faces which do not 
  have any vertices with a -1 (write it to an .obj file- easier than the nasty display function)
- Figure out a better way or better library for displaying a LOT of polygons at once (scatter the entire 
  array of vertices instead of one at a time, that might help, maybe adjust the bounding box instead of 
  scattering, I dunno)
- Write the exit_pt function
	* Get all the intersection points for the faces not entirely on the left or right. Tag the vertex 
	  intersections as well as the edge intersections. Make a hashable thing for the vertices to make 
	  this easier.
	* How to handle the intersection along a boundary edge case?
	* How to handle the coplanar faces?
	* How to handle the vertex intersection cases in general?
- How to get good planes for testing purposes?
- The walking function should be pretty similar to the 2d intersection case
"""
# Flat plane height 0
plane0 = np.array([[0.1,0.9,0.5,0.5],
				   [0.1,0.5,0.5,0.9],
				   [0, 0, 0, 0]])

plane1 = np.array([[0.1,0.9,0.5,0.5],
				   [0.1,0.5,0.5,0.9],
				   [1,1,1,1]])

# 					0 					1				2				3
smalltests = ['tetrahedron.obj', 'hexahedron.obj', 'cross.obj', 'icosahedron.obj',
			  'minispikey.obj', 'edgeshare.obj']
#					4					5

# 	  			 0				1 			  2			3             
bigtests = ['shuttle.obj', 'trumpet.obj', 'al.obj', 'bunny.obj', 'bigbunny.obj']
# faces:        393           11362	        3442 	   4968         69666
# vertices:     310           11908         3618       2503         34835
# edges:        701           23212         7018       7473         104499
def display(filename):
	print '\n\n\t', filename
	P = Polyhedron()
	P.read(filename)
	#P.print_poly()
	P.display(labels=False)

P = Polyhedron()
P.read('icosahedron.obj')
tag_LR(P, plane0) 

P = Polyhedron()
P.read('edgeshare.obj')
tag_LR(P, plane1)
#display('icosahedron.obj')
#display('shuttle.obj')
#display('bunny.obj')
#for f in smalltests: display(f)
#for f in bigtests: display(f)