#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
#	Title: CSR Mesh Generator
#	Description: 
#		Creates CSR mesh for use in Gauss-Seidel 
#		method, using central difference equation.
#		Call elements - return CSR
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import bc_gen as bc
import numpy as np
import scipy as sp
import scipy.linalg
import math

# nel - number of elements in x and y axes...
def csr_mesh(nel):
# Isolate interior elements and split into blocks inel*inel
	inel = nel-2
	tinel = inel*inel

# Calculate the total number of values needed for a dense CSR matrix...
	tvals = tinel*5-(4*(inel-2))-8

# Init with nan so we can tell later if there has been a value unitiated.
# Also entry [0,:] will hold the value and [1,:] will hold the column pointer
	vals =  np.empty([2,tvals])*np.nan
	row = np.zeros(tinel+1)
	cnt = 0 
	pnt = 1
	icnt = 0

	row[0] = 0
# Now place the finite-difference approximation into the CSR form...
	for i in xrange(inel):
		for j in xrange(inel):
			if i==0 and j==0: # Define first corner - connectivity = 3
				vals[0][0] = -4
				vals[0][1] = 1
				vals[0][2] = 1
				vals[1][0] = cnt
				vals[1][1] = cnt+1
				vals[1][2] = cnt+inel
				icnt += 3
			if i==0 and 0<j and j<inel-1: # Define bottom row - connectivity = 4
				row[pnt] = icnt
				vals[0][icnt] = 1
				vals[0][icnt+1] = -4
				vals[0][icnt+2] = 1
				vals[0][icnt+3] = 1
				vals[1][icnt] = cnt-1
				vals[1][icnt+1] = cnt
				vals[1][icnt+2] = cnt+1
				vals[1][icnt+3] = cnt+inel
				icnt += 4
				pnt += 1
			if i==0 and j==inel-1: # Define second corner - connectivity = 3
				row[pnt] = icnt
				vals[0][icnt] = 1
				vals[0][icnt+1] = -4
				vals[0][icnt+2] = 1
				vals[1][icnt] = cnt-1
				vals[1][icnt+1] = cnt
				vals[1][icnt+2] = cnt+inel
				icnt += 3
				pnt += 1
			if j==0 and i>0 and i<inel-1: # Define first column - connectivity = 4
				row[pnt] = icnt
				vals[0][icnt] = -4
				vals[0][icnt+1] = 1
				vals[0][icnt+2] = 1
				vals[0][icnt+3] = 1
				vals[1][icnt] = cnt
				vals[1][icnt+1] = cnt+1
				vals[1][icnt+2] = cnt+inel
				vals[1][icnt+3] = cnt-inel
				icnt += 4
				pnt += 1
			if 0<i and i<inel-1 and 0<j and j<inel-1: # Define interior - connectivity = 5
				row[pnt] = icnt
				vals[0][icnt] = -4
				vals[0][icnt+1] = 1
				vals[0][icnt+2] = 1
				vals[0][icnt+3] = 1
				vals[0][icnt+4] = 1
				vals[1][icnt] = cnt
				vals[1][icnt+1] = cnt-1
				vals[1][icnt+2] = cnt+1
				vals[1][icnt+3] = cnt-inel
				vals[1][icnt+4] = cnt+inel
				icnt += 5
				pnt += 1
			if i==inel-1 and j==0: # Define third corner - connectivity = 3
				row[pnt] = icnt
				vals[0][icnt] = -4
				vals[0][icnt+1] = 1
				vals[0][icnt+2] = 1
				vals[1][icnt] = cnt
				vals[1][icnt+1] = cnt+1
				vals[1][icnt+2] = cnt-inel
				icnt += 3
				pnt += 1
			if i==inel-1 and j==inel-1: # Define final/fourth corner - connectivity = 3
				row[pnt] = icnt
				vals[0][icnt] = -4
				vals[0][icnt+1] = 1
				vals[0][icnt+2] = 1
				vals[1][icnt] = cnt
				vals[1][icnt+1] = cnt-1
				vals[1][icnt+2] = cnt-inel
				icnt += 3
				pnt += 1
			if i==inel-1 and 0<j and j<inel-1: # Define top row - connectivity = 4
				row[pnt] = icnt
				vals[0][icnt] = 1
				vals[0][icnt+1] = -4
				vals[0][icnt+2] = 1
				vals[0][icnt+3] = 1
				vals[1][icnt] = cnt-1
				vals[1][icnt+1] = cnt
				vals[1][icnt+2] = cnt+1
				vals[1][icnt+3] = cnt-inel
				icnt += 4
				pnt += 1
			if j==inel-1 and 0<i and i<inel-1: # Define second column - connectivity = 4
				row[pnt] = icnt
				vals[0][icnt] = 1
				vals[0][icnt+1] = -4
				vals[0][icnt+2] = 1
				vals[0][icnt+3] = 1
				vals[1][icnt] = cnt-1
				vals[1][icnt+1] = cnt
				vals[1][icnt+2] = cnt+inel
				vals[1][icnt+3] = cnt-inel
				icnt += 4
				pnt += 1
			cnt += 1
		row[pnt] = icnt
	# Convert matrix to integers...
	rows = np.int_(row)
	# val = np.int_(vals)
	return(vals, rows)

def test_mesh(nel):
	# Isolate interior elements and split into blocks inel*inel
	inel = nel-2
	tinel = inel*inel

	# Create matrix blocks diagonal, upper, and lower the sum
	Bd = -4*np.eye(inel)
	Bu = np.diag([1]*(inel-1),1)
	Bl = np.diag([1]*(inel-1),-1)
	B = Bd+Bu+Bl

	# Now we need to create entire matrix with inner blocks B
	# so list B into as many blocks as there are inel into matrix A
	Bl = [B]*inel
	A = sp.linalg.block_diag(*Bl)


	# Now fill in the upper stencil which has a offset of inel
	# then repeat for lower diagonals with offset minus inel.
	DiagU = np.diag(np.ones(inel*(inel-1)),inel)
	DiagL = np.diag(np.ones(inel*(inel-1)),-inel)

	# Sum this with main diagonal matirx A to create the linear algebra eq.
	A += DiagU + DiagL
	return(A)

def embed(T):
	N = T.shape[0] + 2
	Tfull = np.zeros([N,N])
	Tfull = bc.dom_init(Tfull, N)
	Tfull[1:-1, 1:-1] = T
	return Tfull