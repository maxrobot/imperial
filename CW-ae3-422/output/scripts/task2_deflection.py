#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
# 	Plots Analytical Solutions for AE3-422
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


# Concetrated Load
for x in xrange(1,6):
	val = 2000*x
	print(val)
	plt.title('Mid-point displacement task 2')
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.ylabel('Deflection (mm)')
	plt.xlabel('Time (s)')
	plt.grid(b=True, which='major', color='k', linestyle='--')
	crs = open('./output/data/task2_deflection' + str(val) + '.txt','r')
	D = np.loadtxt(crs)
	x2 = np.linspace(0,2,len(D)) 
	plt.plot(x2,D,label=val)
	plt.legend(bbox_to_anchor=(1, 0), loc=4, borderaxespad=0., title="Solution")

plt.savefig('./output/images/task2_deflection',bbox_inches='tight')