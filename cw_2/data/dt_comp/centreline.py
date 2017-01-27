#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# Create Subplot
ax = plt.subplot(111)
crs = open('./output0003', 'r')
D = np.loadtxt(crs)
crs2 = open('./output0001','r')
E = np.loadtxt(crs2)
plt.plot(E[:][65],'r-',lw='2',label='Adams-Bashforth, CFL = 0.25')
plt.plot(D[:][65],'k--',lw='2.5',label='Runge-Kutta, CFL = 0.75')

# Add labels etc...
plt.axis([1,129,-.1,.05])
plt.title('Centreline Velocity')
plt.ylabel('Streamwise Velocity')
plt.xlabel('X axis')
plt.legend(bbox_to_anchor=(1, 0), loc=4, borderaxespad=0.)
plt.grid(True)
plt.savefig('./velocity')