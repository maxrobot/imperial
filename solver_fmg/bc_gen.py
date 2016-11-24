#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
#	Title: BC Generator
#	Description: 
#		Creates boundary conditions for mesh...
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import numpy as np
import math

def bc_init(nel, dx):
	# create matrix/vector
	inel = nel-2 # interior elements
	tinel = inel*inel # total interior elements
	# initialiase the matrix b
	B = np.zeros(tinel)
	
	# Fill in boundary conditions
	cnt = 0
	for i in range(0,inel):
		for j in range(0,inel):
			if i==0:
				B[cnt]+=-(1+dx[j])
			if i==(inel-1):
				B[cnt]+=-1
			if j==0:
				B[cnt]+=-1
			if j==(inel-1):
				B[cnt]+=-(math.cos(1.5*math.pi*dx[i])+1)
			cnt+=1
	return B

def dom_init(fdom, nel):
	# Initialise the domain for visualisationsss
	ax = np.linspace(0,1,nel)
	for i in xrange(nel):
		for j in xrange(nel):
			if i==0:
				fdom[i][j]=(1+ax[j])
			if i==(nel-1):
				fdom[i][j]=1
			if j==0:
				fdom[i][j]=1
			if j==(nel-1):
				fdom[i][j]=(math.cos(1.5*math.pi*ax[i])+1)
	return fdom