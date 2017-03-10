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

plt.title('Deflection due to dist. load')
plt.ylabel('Deflection (m)')
plt.xlabel('Axial Position (m)')
plt.grid(b=True, which='major', color='k', linestyle='--')
plt.plot(x,u,'k')
plt.savefig('./output/images/dist_load')
plt.clf()

# Concentrated Load
x2 = np.linspace(0,5000,50)
u2 = np.zeros(500)

u2 = ((1000*x2*x2*(30000-4*x2)))/(48*E*I)

plt.title('Deflection due to concentrated load')
plt.ylabel('Deflection (m)')
plt.xlabel('Axial Position (m)')
plt.grid(b=True, which='major', color='k', linestyle='--')
plt.plot(x2,u2,'k')
plt.savefig('./output/images/conc_load')
plt.clf()

# Concetrated Load
plt.title('Computed Deflection Task 1')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Deflection (m)')
plt.xlabel('Axial Position (m)')
plt.grid(b=True, which='major', color='k', linestyle='--')
crs = open('./output/data/task1_.txt','r')
D = np.loadtxt(crs)
x2 = np.linspace(0,10,len(D))
plt.plot(x2,D[:,1],'k')
plt.savefig('./output/images/computed_load',bbox_inches='tight')
