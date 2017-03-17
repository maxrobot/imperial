#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
# 	Plots Analytical Solutions for AE3-422
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# Concentrated Load
plt.title('Computed Deflection Task 1')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Deflection (mm)')
plt.xlabel('Axial Position (m)')
plt.grid(b=True, which='major', color='k', linestyle='--')
crs = open('./output/data/task_1.txt','r')
D = np.loadtxt(crs)
x2 = np.linspace(0,10,len(D))
plt.plot(x2,D[:,1],'k')
plt.savefig('./output/images/task1',bbox_inches='tight')
