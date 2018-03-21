import numpy as np
from basic_transforms_piv_calibration import *
from Covariance_Tree import *


def test_compute_cov_frame():
	tris = []
	tris.append(Triangle(np.array([1,3,5])  , np.array([3,4,8])  , np.array([-1,6,2])  ))
	tris.append(Triangle(np.array([2,4,9])  , np.array([-3,4,5])  , np.array([-2,4,8])  ))
	tris.append(Triangle(np.array([-1,6,-5])  , np.array([2,4,7])  , np.array([-3,5,2])  ))
	node = Cov_Tree_Node(tris,3)

	node.compute_cov_frame(tris,3)

	result = np.array([  [ 0.27602148 ,-0.960946 , 0.0198729 ]  , [  -0.17158428 ,-0.02892068 , 0.98474485 ]  ,  [  0.94571189 , 0.27522061 , 0.17286596  ] ])
	result_p = np.array([-0.22222222 , 4.44444444 , 4.55555556])
	for val in range(len(node.R_f)):
		for val2 in range(len(node.R_f[0])):
			if abs(result[val][val2] - node.R_f[val][val2]) > .01:
				return False

	for val in range(len(result_p)):
		if abs(result_p[val] - node.p_f[val]) > .01:
			return False
	
	return True

def test_extract_points():
	tris = []
	tris.append(Triangle(np.array([1,3,5])  , np.array([3,4,8])  , np.array([-2,4,8])  )) #duplicate to make sure it works
	tris.append(Triangle(np.array([2,4,9])  , np.array([-3,4,5])  , np.array([-2,4,8])  ))
	tris.append(Triangle(np.array([-1,6,-5])  , np.array([2,4,7])  , np.array([-3,5,2])  ))
	node = Cov_Tree_Node(tris,3)

	points = node.extract_points(tris,3)
	if len(np.array(list(points))[0]) != 8:
		return False #points not unique; had duplicate
	return True

def test_centroid():
	tris = []
	tris.append(Triangle(np.array([1,3,5])  , np.array([3,4,8])  , np.array([-2,4,8])  )) #duplicate to make sure it works
	tris.append(Triangle(np.array([2,4,9])  , np.array([-3,4,5])  , np.array([-2,4,8])  ))
	tris.append(Triangle(np.array([-1,6,-5])  , np.array([2,4,7])  , np.array([-3,5,2])  ))
	node = Cov_Tree_Node(tris,3)

	points = node.extract_points(tris,3)

	result = [-0.3333 , 4.2222 , 5.222]
	cent = node.centroid(points)

	for i in range(len(result)):
		if abs(result[i] - cent[i]) > .01:
			return False

	return True


if not test_compute_cov_frame():
	print "Compute Cov Frame Test failed."
if not test_extract_points():
	print "Extract Points test failed."
if not test_centroid():
	print "Centroid Test failed."

