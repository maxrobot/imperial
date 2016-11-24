#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
#	Title: Gauss-Seidel Solver
#	Description: 
#		Solves the equation using Gauss-Seidel.
#----------------------------------------------------------#
# import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import mesh_gen as mg
import bc_gen as bc
import solvers as run
import store_file as sve
import optparse
import time

def arg():
	# Define number of elements in x and y axes and also the numbe rof iterations in x1000
	usage = "usage: %prog [option1] arg1 [option2] arg2 [option3] arg3"
	parser = optparse.OptionParser(usage=usage)
	# parser = optparse.OptionParser()
	parser.add_option('-n', action="store", type="int",help="Number of elements - x,y")
	parser.add_option('-i', action="store", type="int",help="Maximum iterations (x1000)")
	parser.add_option('-s', action="store", type="int",help="Gauss-Seidel = 1: SOR = 2")
	options, args = parser.parse_args()
	return options.n, options.i, options.s

def main():
	# Define your variables...
	n, it, scheme = arg() # parsed from command line
	nel = n+1 # num of elements
	tnel = nel*nel # tot elements
	inel = nel-2
	tinel = inel*inel
	ax = np.linspace(0,1,nel)
	itar = it*1000 # iterations
	A = np.zeros([tnel,tnel])

	# Initialise parameters
	V, R = mg.csr_mesh(nel)
	b = bc.bc_init(nel, ax)
	x = np.zeros(tinel)
	x[0] = 1. # set first element to 1 according to gauss-seidel

	# Now perform Gauss-Seidel...
	t0 = time.time()

	if scheme == 1:
		x, emax, t = run.gauss_seidel(V,R,x,b,itar,tinel)
	elif scheme == 2:
		x, emax, t = run.sor(V,R,x,b,itar,tinel)
	else:
		print('Error: No scheme selected, 1 = Gauss-Seidel, 2 = Successive Overrelaxation')
		quit()
		
	t1 = time.time()
	t_total = t1 - t0
	titer = t_total/t

	print('Time taken to solution(s): ' + str(t_total).rjust(10) + ' - Time per iteration (s): ' + str(titer).rjust(10) + ' - Error Max: ' + str(emax).rjust(15))
	# t = np.linalg.solve(At,b)

	# Put back into matrix form...
	x_mat = x.reshape(inel,inel)
	S = mg.embed(x_mat)

	# Store data and plot the contour...
	sve.save(S, nel)
	sve.plot_contour(x_mat,n)

if __name__ == "__main__":
   main()