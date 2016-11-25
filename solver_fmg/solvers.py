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
def gauss_seidel(V, R, x, b, it, iels,n):
	e1 = 0.
	e2 = 0.	
	tres = []
	val = []
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
					dval = 1./V[0][diag]
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
		if t>1:
			tres.append(Res/emax)
			val.append(x[iels/2])
		if t%250==0:
			print('Residual: ' + str(Res).rjust(25) + ' Tstep: ' + str(t).rjust(8))
		# if Res<1E-6:
		# 	print('Residual: ' + str(Res).rjust(25) + ' Tstep: ' + str(t).rjust(8))
		# 	break
	np.savetxt('./data/residual' + str(n) + '.txt', tres, delimiter=',')
	np.savetxt('./data/value' + str(n) + '.txt', val, delimiter=',')
	return x, emax, t

# V- values, R - rows, x - solution, it -iterations, iels - int. eles.
def sor(V, R, x, b, it, iels,n):
	e1 = 0.
	e2 = 0.
	# Now perform Successive over relaxation
	# First ten steps calculate omega
	om = 1.
	dxk = 0.
	dxk2 = 0.
	tres = []

	# First loop calculates initial delta x
	for t in xrange(20):
		em = 0.
		xk1 = 0.
		xk2 = 0.
		for i in xrange(iels):
			s1 = 0.
			xk2 += x[i]
			e2 =	 x[i]
			for j in xrange(R[i],R[i+1]):
				pnt = V[1][j]
				if i==pnt:
					diag = j
					dval = om/V[0][diag]
				else:
					s1 += V[0][j]*x[pnt]
					x[i] = dval*(b[i]-s1) + (1. - om)*e2
			xk1 += x[i]
			e1 = x[i]
			ep = abs(e1-e2)
			if ep>em:
				em = ep
		Res = em
		# Calculate the delta x at various points...
		if t==1:
			emax = Res
		if t>1:
			tres.append(Res/emax)
		if t==8:
			ep1 = abs(xk2-xk1)
		if t==9:
			ep2 = abs(xk2-xk1)
			dxk = abs(ep2-ep1)
		if t==18:
			ep1 = abs(xk2-xk1)
		if t==19:
			ep2 = abs(xk2-xk1)
			dxk2 = abs(ep2-ep1)
	
	# Calculate optimal omega
	val = 1. - (dxk2/dxk)**(.1)
	om = 2/(1. + np.sqrt(val))

	# Now iterate rest of loops
	for t in xrange(20,it):
		em = 0.
		for i in xrange(iels):
			s1 = 0.
			e2 = x[i] # e2 is also useful for the SOR
			for j in xrange(R[i],R[i+1]):
				pnt = V[1][j]
				if i==pnt:
					diag = j
					dval = om/V[0][diag]
				else:
					s1 += V[0][j]*x[pnt]
					x[i] = dval*(b[i]-s1) + (1. - om)*e2
			e1 = x[i]
			ep = abs(e1-e2)
			if ep>em:
				em = ep
		Res = em
		tres.append(Res/emax)
		if t%250==0:
			print('Residual: ' + str(em).rjust(25) + ' Tstep: ' + str(t).rjust(8))
		if Res<1E-6:
			print('Residual: ' + str(Res).rjust(25) + ' Tstep: ' + str(t).rjust(8))
			break
	np.savetxt('./data/residual' + str(n) + '.txt', tres, delimiter=',')
	return x, emax, t