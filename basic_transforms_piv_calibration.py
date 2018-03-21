import numpy as np
from scipy import linalg
import math


def registration_3d_3d(A_1,B_1): #passed in as N x 3 since easier to work with
	A = A_1.copy()
	B = B_1.copy()

	A_sum = np.zeros((1,3))
	for a in A:
		A_sum = np.add(A_sum,a)
	a_mean = np.multiply(A_sum,(1.0/len(A))) #sum/N

	B_sum = np.zeros((1,3))
	for b in B:
		B_sum = np.add(B_sum,b)
	b_mean = np.multiply(B_sum,(1.0/len(B))) #sum/N

	#subtract average from each A
	for i in range(len(A)):
		A[i]= np.subtract(A[i],a_mean)
		
	#subtract average from each B
	for i in range(len(B)):
		B[i]= np.subtract(B[i],b_mean)
		
  	"""Start Quaternion method to solve for R"""
  	H = np.zeros((3,3))
  	for i in range(len(A)):
  		a = A[i] #3x3
  		b = B[i] #3x3
  		temp = np.array([  [  np.dot(a[0], b[0]) , np.dot(a[0],b[1])  , np.dot(a[0],b[2]) ], [ np.dot(a[1],b[0])  , np.dot(a[1],b[1])  , np.dot(a[1],b[2]) ],  [  np.dot(a[2],b[0]) , np.dot(a[2],b[1])  , np.dot(a[2],b[2])] ])
  		H = np.add(H,temp)

  	#step 2 -- getting trace and delta to create our G matrix
  	trace_h = np.trace(H)
  	delta_t = np.array([H[1,2] - H[2,1], H[2,0] - H[0,2], H[0,1] - H[1,0]])
  	lower_right = np.subtract(np.add(H,H.transpose()), np.multiply(trace_h, np.identity(3)))

  	G_top = np.append([trace_h],delta_t)
  	delta_t = delta_t.reshape(1,3)
  	G_bottom = np.concatenate((delta_t.transpose(), lower_right),axis=1)
  	#ending with a 4x4 matrix
  	G = np.concatenate((G_top.reshape(4,1).transpose(),G_bottom),axis=0)


  	#step 3 - eigenvalue decomposition
  	lamda,Q = np.linalg.eig(G)

  	#step 4 - find index of largest lambda and use that eigenvector
  	index = lamda.tolist().index(max(lamda))
  	q = Q[:,index] #unit quaternion

  	R = np.array([  [q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2,  2*(q[1]*q[2] - q[0]*q[3]), 2*(q[1]*q[3] + q[0]*q[2]) ]   , [2*(q[1]*q[2] + q[0]*q[3]), q[0]**2 - q[1]**2 + q[2]**2-q[3]**2, 2*(q[2]*q[3] - q[0]*q[1])], [ 2*(q[1]*q[3] - q[0]*q[2]), 2*(q[2]*q[3] + q[0]*q[1]) , q[0]**2-q[1]**2-q[2]**2+q[3]**2]  ])
  	p = np.subtract(b_mean.transpose(), np.dot(R, a_mean.transpose()))
  	return R,p



def transform(R,p,c): #where R is the rotation, p is the translation, c is the vector
	p = p.reshape(3,1)
	c = c.reshape(3,1)

	return np.add(np.dot(R,c),p)

def inv_transform(R,p,c): #where R is the rotation, p is the translation, c is the vector
	p = p.reshape(3,1)
	c = c.reshape(3,1)
	return np.subtract(np.dot(np.linalg.inv(R),c), np.dot(np.linalg.inv(R),p))


def read_coords(num,f): #num coords to read in, pointer to file 
	"""this has been changed for the spaces"""
	coords = np.zeros((num,3))
	for i in range(num):
		line = f.readline().replace(',','').strip()
		if len(line) == 0:
			continue
		line_list = np.array(line.split()).astype(np.float)
		coords[i] = line_list[0:3]
	return coords


def find_closest_point(a,p,q,r): #takes a 3d point and 3 vectors representing corners in a triangle 
	#make sure that they are 3 x 1
	a = a.reshape((3,1))
	p = p.reshape((3,1))
	q = q.reshape((3,1))
	r = r.reshape((3,1))
	b = np.subtract(a,p) #create b for least squares
	A = np.concatenate((np.subtract(q,p),np.subtract(r,p)),axis=1) #create A for least sqaures
	X = np.dot(np.linalg.pinv(A,rcond=1e-20),b) #X is a 2 by 1 matrix (first one lambda and second mu)
	c = np.add(np.add(p,X[0]*(np.subtract(q,p))),X[1]*(np.subtract(r,p))) #calculate what c will be if all lower conditions met
	if X[0] < 0:
		c = project_on_segment(c,r,p) #project on r p
	elif X[1] < 0:
		c = project_on_segment(c,p,q) #project on p q
	elif X[0] + X[1] > 1:
		c = project_on_segment(c,q,r) #project on q r
	return c

def project_on_segment(c,p,q):
	numer = np.dot(np.subtract(c,p).transpose(),np.subtract(q,p)) #(c - p) dot (q - p)
	denom = np.dot(np.subtract(q,p).transpose(),np.subtract(q,p)) #(q - p) dot (q - p)
	lamb = np.divide(numer,denom) #calculate our lambda
	lamb_star = max(0,min(lamb,1.0))
	c_star = np.add(p,lamb_star*np.subtract(q,p)) #find point on the border of the triangle
	return c_star

def distance(a,b): #takes in two points (3D) and computes distance between them (euclidian distance)
	return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)**(1.0/2)

def print_part3(s,c,mag, num, param):
	print num,
	print " ",
	print param + '-Output.txt'
	for i in range(len(c)):
		a = s[i]
		b = c[i]
		print "  ",
		print str(format(float(a[0]),'.2f')).rjust(6), 
		print "    ",
		print str(format(float(a[1]),'.2f')).rjust(6),
		print "    ",
		print str(format(float(a[2]),'.2f')).rjust(6),

		print "       ",
		print str(format(float(b[0]),'.2f')).rjust(6), 
		print "    ",
		print str(format(float(b[1]),'.2f')).rjust(6),
		print "    ",
		print str(format(float(b[2]),'.2f')).rjust(6),

		print "     ",
		print str(format(float(mag[i]),'.3f')).rjust(6)  





