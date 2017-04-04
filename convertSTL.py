# Convert ASCII STL file into list of points with normals

import sys
import numpy as np
import math


def trianglearea(x1, x2, x3):
	k = np.linalg.norm(crossproduct(x3-x1, x3-x2), 2)
	return 0.5*k*k

def crossproduct(u, v):
	return np.array([u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]])


filename = sys.argv[1]			# First argument is filename (without .stl extension)
Ntotal = int(sys.argv[2])		# Second argument is approximate number of desired points
STLfile = open(filename + '.stl')
PNfile = open(filename + '.PN', 'w')
lines = STLfile.readlines()


# Compute the total surface area on the STL file

area = 0.0
counter = 0

for i,line in enumerate(lines):

	words = line.split()
	if len(words) > 0:
		if words[0] == "facet":

			words2 = lines[i+2].split()
			words3 = lines[i+3].split()
			words4 = lines[i+4].split()
		
			x1 = np.array([float(words2[1]), float(words2[2]), float(words2[3])])
			x2 = np.array([float(words3[1]), float(words3[2]), float(words3[3])])
			x3 = np.array([float(words4[1]), float(words4[2]), float(words4[3])])
		
			area += trianglearea(x1,x2,x3)
			counter += 1


# Return to beginning of STL file and seed N particles per area

print("This STL has " + str(counter) + " faces, giving a total surface area of " + str(area) + ".")

particledensity = Ntotal/area

print("The target number of particles is " + str(Ntotal) + ", giving target particle density of " + str(particledensity) + ".")

counter = 0

for i,line in enumerate(lines):
	words = line.split()
	
	if len(words) > 0:
		if words[0] == "facet":
			
			normal = np.array([float(words[2]), float(words[3]), float(words[4])])
			
			words2 = lines[i+2].split()
			words3 = lines[i+3].split()
			words4 = lines[i+4].split()
			
			x1 = np.array([float(words2[1]), float(words2[2]), float(words2[3])])
			x2 = np.array([float(words3[1]), float(words3[2]), float(words3[3])])
			x3 = np.array([float(words4[1]), float(words4[2]), float(words4[3])])
			
			thistriarea = trianglearea(x1,x2,x3)
			
			
			# Now populate triangle with N markers - depends on area of this triangle
					
			N = math.ceil(thistriarea*particledensity)

			# First marker goes on the centroid

			newx = 0.333333*x1 + 0.333333*x2 + 0.333333*x3
			newline = str(normal[0]) + " " + str(normal[1]) + " " + str(normal[2]) + " " + str(newx[0]) + " " + str(newx[1]) + " " + str(newx[2])
			PNfile.write(newline)
			PNfile.write('\n')
			counter += 1
			thisN = 1
			
			while (thisN < N):

				# Randomly seed 'N' markers in this triangle with uniform distribution

				a = np.random.uniform()
				b = np.random.uniform()
				
				if (a + b <= 1.0):

					newx = x1 + a*(x2 - x1) + b*(x3 - x1)
					newline = str(normal[0]) + " " + str(normal[1]) + " " + str(normal[2]) + " " + str(newx[0]) + " " + str(newx[1]) + " " + str(newx[2])
					PNfile.write(newline)
					PNfile.write('\n')
					counter += 1
					thisN += 1
				

print("Created point cloud of " + str(counter) + " points.")

