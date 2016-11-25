#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
#	Title: residual with iterations
#	Description: 
#		Plots residual for every iteration
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys


# Load true solution
res = np.loadtxt('../data/residual32.txt', delimiter='\n')
val = np.loadtxt('../data/value32.txt', delimiter='\n')
n = len(res)

plt.plot(res)
plt.plot(val)

plt.figure()

# sp1
plt.subplot(121)
plt.plot(res)
plt.title('Residual change with iterations')
plt.ylabel('Residual')
plt.xlabel('Iterations')
plt.yscale('log')
# plt.legend(bbox_to_anchor=(1, 0), loc=4, borderaxespad=0., title="Number of Elements")
plt.grid(True)

# sp2
plt.subplot(122)
plt.plot(val)
plt.grid(True)

plt.savefig('residual_iteration.svg')