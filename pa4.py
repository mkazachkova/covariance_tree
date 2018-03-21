"""
************************************************************************************
Programming 4 Assignment -- Main class
CIS
Michael Mudgett
Mariya Kazachkova
************************************************************************************
"""

import sys
import numpy as np
from basic_transforms_piv_calibration import *
from Covariance_Tree import *
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from timeit import default_timer


"""Input the body definition files."""

a3_file = open('PA234_Student_Data/Problem3-BodyA.txt','r')
b3_file = open('PA234_Student_Data/Problem3-BodyB.txt','r')

first_line_split_a = a3_file.readline().strip().split(' ')
num_markers_a = int(first_line_split_a[0])
first_line_split_b = b3_file.readline().strip().split(' ')
num_markers_b = int(first_line_split_b[0])
#get coordinated for A
markers_body_coords_A = read_coords(num_markers_a, a3_file)
A_tip = np.array(a3_file.readline().strip().split()).astype(float)
#get coordinates for B
markers_body_coords_B = read_coords(num_markers_b, b3_file)
B_tip = np.array(a3_file.readline().strip().split()).astype(float) 


"""For each sample frame k, input the values of ai,k and bi,k """
sample_file = open(sys.argv[1] + '-SampleReadingsTest.txt')
first_line_split_sample = sample_file.readline().strip().replace(',',' ').split()
n_total = int(first_line_split_sample[0]) #total number of markers
num_sample_frames = int(first_line_split_sample[1]) #total number of sample frames
n_dummy = n_total - num_markers_a - num_markers_b #number of dummmy markers

s_k_s = np.zeros((num_sample_frames,3))
d_k_s = np.zeros((num_sample_frames,3))

#initial guesses
R_F_reg = np.identity(3)
p_F_reg = np.array([0,0,0])

for i in range(num_sample_frames):
	a_coords = read_coords(num_markers_a, sample_file)
	b_coords = read_coords(num_markers_b, sample_file)
	dummy_coords = read_coords(n_dummy, sample_file)

	"""determine the poses FA,k and FB,k"""
	R_a_k, p_a_k = registration_3d_3d(markers_body_coords_A, a_coords)
	R_b_k, p_b_k = registration_3d_3d(markers_body_coords_B, b_coords)

	d_k = inv_transform(R_b_k, p_b_k, transform(R_a_k,p_a_k, A_tip)) 
	d_k_s[i] = d_k.reshape(1,3)
	"""Compute sample points s"""

mesh_file = open('PA234_Student_Data/Problem3MeshFile.sur', 'r')
first_line_split_mesh = mesh_file.readline().strip().replace(',',' ').split()
num_vertices = int(first_line_split_mesh[0]) #number of vertices
vertices_coords = read_coords(num_vertices,mesh_file) #number of coordinates

second_line_split_mesh = mesh_file.readline().strip().replace(',',' ').split()
num_triangles = int(second_line_split_mesh[0]) #the number of triangles
indices = read_coords(num_triangles,mesh_file)

"""******************************Brute force method******************************"""

for i in range(num_sample_frames):
	s_k = transform(R_F_reg, p_F_reg, d_k_s[i]) #should be same as dk because the transformation doesn't do anything currently
	s_k_s[i] = s_k.reshape((1,3)) #store our sample point values


#start_iterative = default_timer()

"""Now find the points C on the surface mesh that are closest to the sk""" 
c_closest = np.zeros(s_k_s.shape)
min_distances = np.zeros((num_sample_frames,1))
c_k_s = np.zeros((num_sample_frames,3))
min_triangle = 0
#s_k_s[len(s_k_s)-1] = np.array([1,1,1])
"""for s_k in range(len(s_k_s)):
	min_distance = 10000 #just start with a large number 
	for triangle in indices:
			#c_temp is the possible closest point
			c_temp = find_closest_point(s_k_s[s_k], vertices_coords[int(triangle[0])], vertices_coords[int(triangle[1])], vertices_coords[int(triangle[2])])
			temp_distance = distance(s_k_s[s_k],c_temp)
			#check if c_temp is the current closest point
			if temp_distance < min_distance:
				min_distance = temp_distance
				min_distances[s_k] = min_distance
				#set c_temp as the current closest point
				c_closest[s_k] = c_temp.reshape((1,3))
				min_triangle = int(triangle[0])	"""						

"""Print all results in the desired format (matching provided output file"""
#print_part3(s_k_s,c_closest,min_distances,num_sample_frames, sys.argv[1])
#print c_closest

#duration_iterative = default_timer() - start_iterative #stop timing
#print "One iteration of iterative method took: " + str(duration_iterative) + " seconds."





"""******************************Covariance Tree method******************************"""

#need to create our list of triangles
tris = []
for t in indices:
	tri = Triangle(vertices_coords[int(t[0])],vertices_coords[int(t[1])],vertices_coords[int(t[2])])
	tris.append(tri)

node = Cov_Tree_Node(tris,len(tris))

all_true = False
old_distances = np.full((len(s_k_s),1),10000) #store previous distances
new_distances = np.zeros((len(s_k_s),1)) #store new calculated distances

average_distances = [] #to make a graph

#initial guesses
R_F_reg = np.identity(3)
p_F_reg = np.array([0,0,0])

#get starting time
start_covaraince = default_timer()

while not all_true:
	for i in range(num_sample_frames):
		s_k = transform(R_F_reg, p_F_reg, d_k_s[i]) 
		s_k_s[i] = s_k.reshape((1,3)) #store our sample point values


	"""find closest point on triangle"""
	c_closest = np.zeros(s_k_s.shape)

	vals = []
	for s_k in range(len(s_k_s)):
		d = distance(np.array(s_k_s[s_k]),np.array([10000,10000,10000]))
		#actually search the tree
		c_closest[s_k] = node.search_tree(np.array(s_k_s[s_k]),np.array([10000,10000,10000]),d).reshape((1,3))

	#find new F reg values
	R_F_reg,p_F_reg  = registration_3d_3d(d_k_s,c_closest)

	all_true = True
	for x in range(len(s_k_s)):
		new_distances[x] = distance(s_k_s[x], c_closest[x])

	average_distances.append(np.mean(np.array(new_distances)))

	"""check to see if need to keep iterating"""
	for x in range(len(s_k_s)):
		s_k_s[x] = transform(R_F_reg, p_F_reg, d_k_s[x]).reshape((1,3))
		if abs(1 - (old_distances[x]/new_distances[x])) > .005: 
			old_distances = new_distances.copy()
			all_true = False
			break

duration_covariance = default_timer() - start_covaraince #stop timing
#print "Covariance method until convergence took: " + str(duration_covariance) + " seconds."

	
print_part3(s_k_s,c_closest,new_distances,num_sample_frames, sys.argv[1])
#pyplot.plot(average_distances)
#pyplot.ylabel('Average distance')
#pyplot.xlabel('Iterations')
#pyplot.grid(True)
#pyplot.title('Average distance Between Sample Points and Closest Points versus Iterations using Covariance Tree')
#pyplot.show()
