#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
#	Title: Solution Plotter
#	Description: 
#		Plots the means of the Gauss-Seidel approx.
#		across X axis...
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys

nel = 4

for x in xrange(1,8):
	# read in file...
	array = np.loadtxt('array'+ str(nel) +'.txt', delimiter=',')
	
	# Determine file length
	n = len(array)
	ax = np.linspace(0,1,n)
	# print(ax[nel/2])
	mean = np.zeros(n)

	# Sum solutions
	for i in xrange(n):
		mean[i] = array[nel/2][i]
	
	plt.plot(ax,mean,label=nel,lw=2)
	nel = nel*2


# Add labels etc...
plt.title('Temperature at y = 0.5')
plt.ylabel('T (K)')
plt.xlabel('x axis')
plt.legend(bbox_to_anchor=(0, 0), loc=3, borderaxespad=0., title="Number of Elements")
# plt.show()
plt.savefig('plot.svg')
