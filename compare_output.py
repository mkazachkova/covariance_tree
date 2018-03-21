import sys
import itertools
import numpy as np
#calculate difference between our output file and the provided output file
our_file = open(sys.argv[1]) #our produced file
given_file = open(sys.argv[2]) #given file

first_line = given_file.readline().strip().split() 
trash = our_file.readline()

error1 = np.zeros((int(first_line[0]),3))
error2 = np.zeros((int(first_line[0]),3))
error3 = np.zeros(((int(first_line[0]),1)))
index = 0

for ours, given in itertools.izip(our_file, given_file):
	ours = np.array(ours.split()).astype(np.float)
	given = np.array(given.split()).astype(np.float)
	error1[index] = abs(np.subtract(ours[0:3],given[0:3])).reshape((1,3))
	error2[index] = abs(np.subtract(ours[3:6],given[3:6])).reshape((1,3))
	error3[index] = abs(np.subtract(ours[6:7],given[6:7])).reshape((1,1))
	index+=1

error1_mean = np.mean(error1, axis=0)
error2_mean = np.mean(error2, axis=0)
error3_mean = np.mean(error3, axis=0)
print str(format(float(error1_mean[0]),'.2f')) + ',',
print str(format(float(error1_mean[1]),'.2f')) + ',',
print str(format(float(error1_mean[2]),'.2f'))

print str(format(float(error2_mean[0]),'.2f')) + ',',
print str(format(float(error2_mean[1]),'.2f')) + ',',
print str(format(float(error2_mean[2]),'.2f'))

print str(format(float(error3_mean[0]),'.2f'))