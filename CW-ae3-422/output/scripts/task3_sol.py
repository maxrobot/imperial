#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
# 	Plots Analytical Solutions for AE3-422
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# Concetrated Load
for x in xrange(1,10):
	pnt = x*1000
	plt.title('Computed Deflection Task 1')
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.ylabel('Deflection (m)')
	plt.xlabel('Axial Position (m)')
	plt.grid(b=True, which='major', color='k', linestyle='--')
	crs = open('./output/data/output' + str(pnt) +'.txt','r')
	D = np.loadtxt(crs)
	x2 = np.linspace(0,10,len(D))
	plt.plot(x2,D[:,1],'k')
	plt.savefig('./output/images/task3',bbox_inches='tight')
	# plt.clf()
