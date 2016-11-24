#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
#	Title: Gauss-Seidel Solver
#	Description: 
#		Solves the equation using Gauss-Seidel.
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import mesh_gen as mg
import bc_gen as bc
import solvers as run
import store_file as sve
import sys
import time

# Define your variables...
n = 16
nel = n+1 # num of elements
tnel = nel*nel # tot elements
inel = nel-2
tinel = inel*inel
ax = np.linspace(0,1,nel)
itar = 100000 # iterations

A = np.zeros([tnel,tnel])
# Initialise parameters
V, R = mg.csr_mesh(nel)
b = bc.bc_init(nel, ax)
x = np.zeros(tinel)
x[0] = 1. # set first element to 1 according to gauss-seidel

# Now perform Gauss-Seidel...
t0 = time.time()

# x, emax, t = run.gauss_seidel(V,R,x,b,itar,tinel)
x, emax, t = run.sor(V,R,x,b,itar,tinel)
t1 = time.time()
t_total = t1 - t0
titer = t_total/t

print('Time taken to solution(s): ' + str(t_total).rjust(10) + ' - Time per iteration (s): ' + str(titer).rjust(10) + ' - Error Max: ' + str(emax).rjust(15))
# S = np.linalg.solve(A,b)

# Put back into matrix form...
x_mat = x.reshape(inel,inel)
S = mg.embed(x_mat)
sve.save(S, nel)

# Create Subplot
ax = plt.subplot(111)
ax.set_title('Temperature (K)')
# colormap will be used for the contour lines
im = plt.imshow(x_mat,interpolation='bilinear', origin='lower',extent=(0,1,0,1))
levels = np.arange(-.1,2.1, 0.1)
CS = plt.contour(x_mat, levels,origin='lower',linewidths=.5,extent=(0,1,0,1),colors='gray')
plt.clabel(CS, inline=1, fontsize=10)

# We can still add a colorbar for the image, too.
CBI = plt.colorbar(im,orientation='horizontal', shrink=0.8)
plt.xlabel('x axis')
plt.ylabel('y axis')
# Save figure...
plt.savefig('./images/temp' + str(n) + '.svg')