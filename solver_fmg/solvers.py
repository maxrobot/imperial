#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
#	Title: Solvers
#	Description: 
#		Solves the equation using Gauss-Seidel or SOR.
#		From the Compressed Sparse Rows, of course!
#----------------------------------------------------------#
import numpy as np
import math
# V- values, R - rows, x - solution, it -iterations, iels - int. eles.
def gauss_seidel(V, R, x, b, it, iels):
	e1 = 0.
	e2 = 0.	
	# Now perform Gauss-Seidel...
	for t in xrange(it):
		em = 0.
		for i in xrange(iels):
			s1 = 0.
			e2 = x[i]
			for j in xrange(R[i],R[i+1]):
				pnt = V[1][j]
				if i==pnt:
					diag = j
					dval = 1/V[0][diag]
				else:
					s1 += V[0][j]*x[pnt]
					x[i] = dval*(b[i]-s1)
			e1 = x[i]
			ep = abs(e1-e2)
			if ep>em:
				em = ep
		Res = em
		if t==1:
			emax = Res
		if t%250==0:
			print('Residual: ' + str(Res).rjust(25) + ' Tstep: ' + str(t).rjust(8))
		if Res<1E-6:
			print('Residual: ' + str(Res).rjust(25) + ' Tstep: ' + str(t).rjust(8))
			break
	return x, emax, t

# V- values, R - rows, x - solution, it -iterations, iels - int. eles.
def sor(V, R, x, b, it, iels):
	e1 = 0.
	e2 = 0.
	# Now perform Successive over relaxation
	# First ten steps calculate omega
	om = 1.
	for t in xrange(it):
		em = 0.
		for i in xrange(iels):
			s1 = 0.
			e2 = x[i]
			for j in xrange(R[i],R[i+1]):
				pnt = V[1][j]
				if i==pnt:
					diag = j
					dval = om/V[0][diag]
				else:
					s1 += V[0][j]*x[pnt]
					x[i] = dval*(b[i]-s1) + (1 - om)*x[i]
			e1 = x[i]
			ep = abs(e1-e2)
			if ep>em:
				em = ep
		Res = em
		if t==1:
			emax = Res
		if t%1==0:
			print('Residual: ' + str(em).rjust(25) + ' Tstep: ' + str(t).rjust(8))
		# if em<1E-8:
		# 	print('Residual: ' + str(em).rjust(25) + ' Tstep: ' + str(t).rjust(8))
		# 	break
	return x, emax, t