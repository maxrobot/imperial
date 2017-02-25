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
plt.ylabel('Deflection (mm)')
plt.xlabel('Axial Position (mm)')
plt.grid(b=True, which='major', color='k', linestyle='--')
plt.plot(x,u,'k')
plt.savefig('dist_load')
plt.clf()

# Concetrated Load
x2 = np.linspace(0,5000,50)
u2 = np.zeros(500)

u2 = ((1000*x2*x2*(30000-4*x2)))/(48*E*I)

plt.plot(x2,u2,'k')
plt.savefig('conc_load')