#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
#	Title: Error Max Vs. Mesh Size
#	Description: 
#		Plots the maximum error against
#		mesh size.
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys

nel = 4

# Load true solution
array = np.loadtxt('array256.txt', delimiter=',')

ubar = np.zeros(9)
for i in xrange(256):
	for j in xrange(256):
		ubar+=array[i,j]

ubar = ubar/(256*256)

for x in xrange(1,8):
	# read in file...
	app = np.loadtxt('array'+ str(nel) +'.txt', delimiter=',')

# 	# Determine file length
	n = len(array)
	mean = np.zeros(n)
	u = 0.
	for i in xrange(nel):
		for j in xrange(nel):
			u+=app[i,j]

	u = u/(nel*nel)
	# Find max error...
	em = abs(ubar - u)

	dx = 1./nel
	print(dx,em,ubar,u)
	plt.plot(dx,em,'o',label=nel)
	nel = nel*2

# Add labels etc...
plt.title('Maximum Error')
plt.ylabel('E max')
plt.xlabel('delta x')
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(1, 0), loc=4, borderaxespad=0., title="Number of Elements")
plt.grid(True)
plt.savefig('error.svg')