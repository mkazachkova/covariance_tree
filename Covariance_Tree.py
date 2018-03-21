
"""
************************************************************************************
Programming 4 Assignment  -- Covariance Tree class
CIS
Michael Mudgett
Mariya Kazachkova
************************************************************************************
"""
import numpy as np
import scipy as sp
from basic_transforms_piv_calibration import *
from scipy import linalg
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from numpy import median

"""
************************************************************************************
Represents a single node in our covariance tree.
************************************************************************************
"""
class Cov_Tree_Node:

	minCount = 25 #the minimum number of triangles each tree MUST have in self.things

	"""
	************************************************************************************
	Our constructor for a single node. 
	Takes in an array of traingles and the length of that array.
	************************************************************************************
	"""
	def __init__(self,ts,nT):
		self.things = ts #array of triangles where each triangle holds 3 vectors
		self.num_things = nT
		self.points = None
		self.R_f = None
		self.p_f = None
		self.U = None
		self.left = None
		self.right = None
		self.haveSubTrees = False
		self.cc = None
		self.pos_cent = None
		self.neg_cent = None
		self.median = 0
	 	self.compute_cov_frame(self.things,self.num_things) #get our transformation
		self.construct_subtrees() #make left and right subtrees or finish 

	
	"""
	************************************************************************************
	This function computes a transformation into the current node's coordinate system.
	It takes an array of triangles and the length of that array.
	************************************************************************************
	"""	
	def compute_cov_frame(self,ts,nT):
		points,num_p = self.extract_points(ts,nT)
		self.points = np.array(list(points)) #change to np array from set
		cent = self.centroid(self.points)
		self.U = self.points - cent
	

		A = np.cov(self.U.transpose())
		lamda,Q = np.linalg.eig(A) #get eigenvalues and vectors
		index = lamda.tolist().index(max(lamda))
		new_Q=Q
		new_lamda = np.array([[lamda[0], 0, 0], [0, lamda[1], 0], [0, 0, lamda[2]]])

		#need to ennsure that we have either a rotation or reflection matrix
		if index==0:
			new_Q = Q
			
		if index==1:
			new_Q = np.array([Q[:,1], Q[:,0], Q[:,2]])	
			
		if index==2:
			new_Q = np.array([Q[:,2], Q[:,1], Q[:,0]])

		R = new_Q #our new matrix
  		self.R_f = R
  		self.p_f = cent 


  	"""
  	************************************************************************************
  	Transforms the array of triangles into a set of points and then returns this set.
  	Uses an array of triangles and the length of that array.
  	************************************************************************************
  	"""
	def extract_points(self,ts,nT):
		extracted = []
		for triangle in ts:
			extracted.append(triangle.p)
			extracted.append(triangle.r)
			extracted.append(triangle.q)
		extracted = np.array(extracted)
		extracted = np.unique(extracted, axis=0)
		return extracted, len(extracted)
		
		
	"""
	************************************************************************************
	Calculates the center of all of the points.
	Uses the extracted list of points from the node's array of triangles.
	************************************************************************************ 
	"""
	def centroid(self,points):
		summed = []
		for triangle in self.things:
			summed.append((triangle.p + triangle.r + triangle.q)/3.0)
		return np.mean(np.array(summed),axis=0) #axis?

	"""
	************************************************************************************
	Call Cov_Tree_Node constructor in order to create the left and right subtrees of 
	the current node.
	If the number of triangles is less than the minimum number allowed, the function is
	 returned and our node gets no subtrees.
	************************************************************************************
	"""
	def construct_subtrees(self):
		if self.num_things < self.minCount:
			self.haveSubTrees = False #this node has no subtrees
			return

		centers = []
		left_list = []
		right_list = []
		pos_centers = []
		neg_centers = []
		new_points = []
		for i in range(len(self.things)):
			centers.append(np.dot(np.linalg.inv(self.R_f),((self.things[i].p + self.things[i].r + self.things[i].q)/3.0) - self.p_f)) #maybe change to np.add
		for i in range(len(self.points)):
			new_points.append(np.dot(np.linalg.inv(self.R_f),((self.points[i]) - self.p_f)))
		self.median = 0
		
		for i in range(len(centers)):
			#get our triangle vertices in terms of the current node's coordinate system
			transformed_p = np.dot(np.linalg.inv(self.R_f),(self.things[i].p - self.p_f))
			transformed_q = np.dot(np.linalg.inv(self.R_f),(self.things[i].q - self.p_f))
			transformed_r = np.dot(np.linalg.inv(self.R_f),(self.things[i].r - self.p_f))

			#if our matrix is a rotation matrix
			if np.linalg.det(self.R_f) > 0:
				if transformed_p[0] <= self.median or transformed_q[0] <= self.median or transformed_r[0] <= self.median:
					left_list.append(self.things[i]) #append triangle
					neg_centers.append(centers[i])
				if transformed_p[0] > self.median or transformed_q[0] > self.median or transformed_r[0] > self.median:
					right_list.append(self.things[i]) #append triangle
					pos_centers.append(centers[i])
			#if our matrix is a reflection matrix
			else:
				if transformed_p[0] > self.median  or transformed_q[0] > self.median or transformed_r[0] > self.median: 
					left_list.append(self.things[i]) #append triangle
					neg_centers.append(centers[i])
				if transformed_p[0] <= self.median or transformed_q[0] <= self.median or transformed_r[0] <= self.median:
					right_list.append(self.things[i]) #append triangle
					pos_centers.append(centers[i])


		self.haveSubTrees = True #we will construct subtrees
		self.pos_cent = pos_centers
		self.neg_cent = neg_centers

		"""if Cov_Tree_Node.count > 12:
			self.three_d_plot(pos_centers, neg_centers)
			exit()"""

		self.left = Cov_Tree_Node(left_list,len(left_list)) #build left tree
		self.right = Cov_Tree_Node(right_list, len(right_list)) #build right tree


	"""
	************************************************************************************
	Plots each subtree and the current sample point in three dimensional space.
	Allows us to visually track which subtree is being entered.
	Takes two lists of points 3d points and a single 3d point.
	************************************************************************************
	"""
	def three_d_plot(self,l,r,v):
		x = []
		y = []
		z = []
		for val in l:
			x.append(val[0])
			y.append(val[1])
			z.append(val[2])
		"""for val in l:
			x.append(val.p[0])
			x.append(val.r[0])
			x.append(val.q[0])
			y.append(val.p[1])
			y.append(val.r[1])
			y.append(val.q[1])
			z.append(val.p[2])
			z.append(val.r[2])
			z.append(val.q[2])"""

		fig = pyplot.figure()
		ax = Axes3D(fig)
		ax.scatter(x,y,z)

		x = []
		y = []
		z = []
		for val in r:
			x.append(val[0])
			y.append(val[1])
			z.append(val[2])
		"""for val in r:
			x.append(val.p[0])
			x.append(val.r[0])
			x.append(val.q[0])
			y.append(val.p[1])
			y.append(val.r[1])
			y.append(val.q[1])
			z.append(val.p[2])
			z.append(val.r[2])
			z.append(val.q[2])"""

		ax.scatter(x,y,z)

		ax.scatter(v[0],v[1],v[2],c='Black',s=150)

		for angle in range(0, 360,2):
			ax.view_init(30, angle)
			pyplot.draw()
			pyplot.pause(.1)


	"""
	************************************************************************************
	Recursive function to go through the tree and return the point closest 
	to the sample point.
	Takes in a sample point, the last closest distance, and the distance between the 
	two aforementioned points.
	************************************************************************************
	"""
	def search_tree(self, a, c,dist): #a = sample point, c = previous closest point
		#put point in terms of current node's coordinate system
		b = inv_transform(self.R_f,self.p_f,a)		
		if self.haveSubTrees:
			#self.three_d_plot(self.pos_cent,self.neg_cent,b) #for viewing in 3d space
			#check if matrix is a rotation matrix
			if np.linalg.det(self.R_f) > 0:
				if b[0] > 0:
					return self.right.search_tree(a,c,dist)
				elif b[0] <= 0:
					return self.left.search_tree(a,c,dist)
			#if matrix is a refelection matrix
			else:

				if b[0] <= 0:
					return self.right.search_tree(a,c,dist)
				elif b[0] > 0:
					return self.left.search_tree(a,c,dist)
		else:
			min_i = 0
			for i in range(len(self.things)):
				cp = find_closest_point(a, self.things[i].p,self.things[i].q,self.things[i].r)
				dst = distance(cp,a)
				if dst < dist:
					dist = dst
					c = cp
					min_i = i
			return c #found the closest point

"""
************************************************************************************
A triangle consists of 3 3d points.
The 'things' in our self.things array.
************************************************************************************
"""
class Triangle:
	def __init__(self,p,q,r):
		self.p = p
		self.q = q
		self.r = r





