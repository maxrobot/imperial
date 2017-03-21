#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
# 	Plots Analytical Solutions for AE3-422
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# Define constant variables
E = 210000
b = 100
h = 120
I = (b*h**3)/12
Fy = 1000

# Distributed Load
x = np.linspace(0,10000,100)
u = np.zeros(100)
u = ((1*x*x*(10000-x)**2))/(24*E*I)


# Concentrated Load
x2 = np.linspace(0,5000,50)
u2 = np.zeros(50)
u2 = ((Fy*x2*x2*(30000-4*x2)))/(48*E*I)

# Sum of both dist and conc forces
u3 = np.zeros(100)
for x in xrange(1,50):
	u3[x] = u[x] + u2[x]
	u3[x+49] = u[x+49] + u2[50-x]
	
# Concentrated Load
plt.title('Computed Deflection Task 1')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Deflection (mm)')
plt.xlabel('Axial Position (m)')
plt.grid(b=True, which='major', color='k', linestyle='--')
crs = open('./output/data/task_1.txt','r')
D = np.loadtxt(crs)
x2 = np.linspace(0,10,len(D))
x3 = np.linspace(0,10,100)
plt.plot(x2,D[:,1],'k',lw=2,label='Computed')
plt.plot(x3,u3,'r--',lw=1.25,label='Analytical')
plt.legend(bbox_to_anchor=(.5, 0), loc=8, borderaxespad=0., title="Solution")
plt.savefig('./output/images/task1',bbox_inches='tight')
