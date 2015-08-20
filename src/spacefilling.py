#!/usr/bin/python
#
# usage: subdivider mesh.ply [--method stl] [--level 1] [--output file_out_1.py]
#
# method:
# stl - all faces are treated as unconnected
# ply - all the vertices are connected and indexed
#
# level: subdivision to use
#
# by yours truly,
# Dalai Felinto
#
# April 2011
# www.dalaifelinto.com

# ##########################
# Defines
# ##########################
class defines:
	class method:
		stl = "stl"
		ply = "ply"

# ##########################
# Imports
# ##########################
import os
import sys
import math
import time
import math
import numpy as np
import random
import argparse
from sys import argv
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

def is_good_file(filename):
	"""XXX error messages"""
	if not os.path.exists(filename):
		return False
	if not filename.endswith(".ply"):
		return False

	return True, 

def parse_args():
	import argparse  # to parse options for us and print a nice help message

	# get the args passed to blender after "--", all of which are ignored by
	# blender so scripts may receive their own arguments
	argv = sys.argv
	file_input = None

	if len(argv) < 2 or not is_good_file(argv[1]): 
		argv = []  # as if no args are passed
	else:
		file_input = argv[1]
		argv = argv[2:]  # get all args after "myfile.ply"

	# When --help or no args are given, print this help
	usage_text = \
	"usage:\n" \
	"subdivider mesh.ply [--method stl] [--level 1] [--output file_out_1.py]\n" \
	"\n" \
	"stl - all faces are treated as unconnected\n" \
	"ply - all the vertices are connected and indexed\n" \
	"\n" \
	"level: subdivision to use\n" \
	"\n"
	parser = argparse.ArgumentParser(description=usage_text)

	# Example utility, add some text and renders or saves it (with options)
	# Possible types are: string, int, long, choice, float and complex.
	parser.add_argument("-m", "--method", dest="method", type=str,
			help="Face structure method")

	parser.add_argument("-l", "--level", dest="level", type=int,
			help="Subdivision level")

	parser.add_argument("-o", "--output", dest="output", type=str,
			help="Output file")

	args = parser.parse_args(argv)  # In this example we wont use the args

	if not argv:
		parser.print_help()
		return 0

	if not file_input:
		print("Error: file argument not given, aborting.")
		parser.print_help()
		return 0

	# Run the example function
	return validate_function(file_input, args.output, args.method, args.level)

	print("Parsing Process Finished.")

def validate_function(input, output, method, level):
	"""For now simply bypass and return the passed values"""

	return input, output, method, level

# ##########################
# I/O
# ##########################

def import_ply(input):
	fileR = open(input, "r")

	inside_header=True
	element_vertex = 0
	element_face = 0
	
	verts = []
	faces = []

	line=fileR.readline()
	while line and not line.startswith("end_header"):

		if line.startswith("end_header"): break

		if line.startswith("element vertex"):
			element_vertex = int(line.split(' ')[2])

		if line.startswith("element face"):
			element_face = int(line.split(' ')[2])

		line=fileR.readline()

	for i in range(element_vertex):
		line = fileR.readline().split(' ')
		verts.append((float(line[0]), float(line[1]), float(line[2])))

	for i in range(element_face):
		line = fileR.readline().split(' ')
		faces.append((int(line[1]), int(line[2]), int(line[3])))

	print("Imported file: %s" % (input))
	print("Number of vertices: %d" % (element_vertex))
	print("Number of faces: %d" % (element_face))
	print("")

	return verts, faces

def export_ply(output, verts, faces):
	
	print("")
	print("Exporting file: %s" % (output))
	print("Number of vertices: %d" % (len(verts)))
	print("Number of faces: %d" % (len(faces)))
	print("")

	fileW = open(output, "w")
	
	fw=fileW.write

	for vert in verts:
		fw("%.6f %.6f %.6f" % vert)
		fw("\n")

	#for face in faces:
		#fw("%d %d %d" % face[:3])
		#var = int(face[0])
		#for vertice in range(len(verts)):
		#vertice = var






	fileW.close()
# ##########################
# Convert
# ##########################

class Vertex(tuple):
	def normalize(self):
		x,y,z = self[0], self[1], self[2]
		length = math.sqrt(x ** 2 + y ** 2 + z ** 2)

		x /= length
		y /= length
		z /= length

		return Vertex((x,y,z))

	def __add__(self, v1):
		return Vertex((self[0] + v1[0], self[1] + v1[1], self[2] + v1[2]))

def convert_stl(raw_verts, raw_faces):
	faces = []
	for face in raw_faces:
		faces.append((Vertex(raw_verts[face[0]]), \
					Vertex(raw_verts[face[1]]), \
					Vertex(raw_verts[face[2]])))

	return faces

def convert_ply(raw_verts, raw_faces):
	def get_edge(v1, v2, edges):
		"""creates initial edge indexing"""
		key = "%020d%020d" % (min(v1,v2), max(v1,v2))

		id = edges.get(key)
		if id == None:
			id = len(edges)
			edges[key] = id
		
		return id

	edges = dict()
	faces = []
	
	for face in raw_faces:
		v1,v2,v3 = face
		e1=get_edge(v2,v3, edges)
		e2=get_edge(v3,v1, edges)
		e3=get_edge(v1,v2, edges)
		
		faces.append((face + (e1, e2, e3)))

	verts =  list(Vertex(vert).normalize() for vert in raw_verts)

	return verts, faces

def convert(raw_verts, raw_faces, method):
	if method == defines.method.stl:
		return None, convert_stl(raw_verts, raw_faces)
	else:
		return convert_ply(raw_verts, raw_faces)

# ##########################
# Subdivide
# ##########################

def subdivide_stl(faces, level):
	for i in range(level):
		new_faces = []
		for v1, v2, v3 in faces:
			v4 = (v1+v2).normalize()
			v5 = (v2+v3).normalize()
			v6 = (v3+v1).normalize()

			new_faces.append(( v1, v4, v6))
			new_faces.append(( v2, v5, v4))
			new_faces.append(( v3, v6, v5))
			new_faces.append(( v4, v5, v6))

		faces = new_faces

	verts = []
	for v1,v2,v3 in faces:
		verts.append(v1)
		verts.append(v2)
		verts.append(v3)

	total=len(faces)
	faces= []
	for i in range(total):
		# to do: create ids for all faces
		off = i * 3
		faces.append((off,off+1,off+2))

	return verts, faces

def subdivide_ply(verts, faces, level):
	for i in range(level):
		time_initial = time.time()

		old_verts = len(verts)
		old_edges = int(len(faces) * 1.5)
		edge_inc = old_edges * 2

		verts.extend(list(None for i in range(old_edges)))
		new_faces = []

		for v1, v2, v3, ea, eb, ec in faces:
			# new vertex indices
			v4 = old_verts + ea
			v5 = old_verts + eb
			v6 = old_verts + ec

			# add new vertex only if non-existent
			if verts[v4] == None:
				verts[v4] = (verts[v2]+verts[v3]).normalize()

			if verts[v5] == None:
				verts[v5] = (verts[v3]+verts[v1]).normalize()

			if verts[v6] == None:
				verts[v6] = (verts[v1]+verts[v2]).normalize()

			# note, the following could be 'optmized' to run the
			# if only once . the gain is of 3%

			# new external edges
			e1 = ea + (old_edges if v2 > v3 else 0)
			e4 = ea + (old_edges if v2 < v3 else 0)

			e2 = eb + (old_edges if v3 > v1 else 0)
			e5 = eb + (old_edges if v3 < v1 else 0)

			e3 = ec + (old_edges if v1 > v2 else 0)
			e6 = ec + (old_edges if v1 < v2 else 0)

			# new internal edges
			e7 = edge_inc; edge_inc += 1
			e8 = edge_inc; edge_inc += 1
			e9 = edge_inc; edge_inc += 1

			new_faces.append(( v1, v6, v5, e7, e5, e3))
			new_faces.append(( v6, v2, v4, e1, e8, e6))
			new_faces.append(( v5, v4, v3, e4, e2, e9))
			new_faces.append(( v6, v4, v5, e9, e7, e8))

		faces=new_faces
		time_final= time.time() - time_initial
		#print("subdivide level %d : %4.4f" % (i+1, time_final))

	return verts, faces

def subdivide(verts, faces, level, method):
	if method == defines.method.stl:
		return subdivide_stl(faces, level)
	else:
		return subdivide_ply(verts, faces, level)

def density_variance(verts, faces):

	area = 0
	var = 0
	for m in range(len(faces)):
		tmp = faces[m]
		v1 = verts[tmp[0]]
		v2 = verts[tmp[1]]
		v3 = verts[tmp[2]]
		area_tri = math.sqrt(((v2[0]*v3[1])-(v3[0]*v2[1]))**2 + ((v3[0]*v1[1]) - (v1[0]*v3[1]))**2 + ((v1[0]*v2[1])-(v2[0]*v1[1]))**2)/2
		area += area_tri
	mean_area = area/float(len(faces))

	for n in range (len(faces)):
		var_tri = (area_tri - mean_area)**2
		var += var_tri

	mean_var = var/float(len(faces))
	#print mean_area
	#print mean_var
	print "Mean area: " + str(mean_area)
	print "Variability: " + str(mean_var)
	
	#print area

def scatter_icosahedron(filename):
    """ This function takes the points from the txt file created by subdiver.py to use for graphing """
    #script, filename = argv

    x, y, z = np.genfromtxt(filename, unpack=True)
    #print x
    return(x, y, z)

def scatter_plot3d(x, y, z, title=None):
    """Scatter plot a set of points in 3D"""
    pyplot.ioff()
    fig = pyplot.figure()
    axes = fig.gca(projection='3d')
    axes.scatter(xs=x, ys=y, zs=z)
    if title:
        pyplot.title(title)
        pyplot.show()


# ##########################
# Main
# ##########################

def main():
	# 1) parse args
	input, output, method,level = parse_args()

	# 2) import ply
	raw_verts, raw_faces = import_ply(input)
	
	# 3) convert to internal structure
	verts, faces = convert(raw_verts, raw_faces, method)
	
	# 4) subdivide
	verts, faces = subdivide(verts, faces, level, method)

	# 5) export ply
	export_ply(output, verts, faces)
	#print verts, faces

	#6 calculate mean area and variance of triangles
	density_variance(verts, faces)
	

	#7 scatter icosahedron graph
	filename = sys.argv[7]
	#print "Density mean = %g, variance = %g" %density_variance(x, y, z)
	(mx, my, mz) = scatter_icosahedron(filename)
	scatter_plot3d(mx, my, mz, "Icosahedron") 
	
if __name__ == 'main':
	main()

main()