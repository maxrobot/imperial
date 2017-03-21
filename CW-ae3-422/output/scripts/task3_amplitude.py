#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
# 	Plots Analytical Solutions for AE3-422
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

amp = np.zeros(11)
x2 = np.zeros(11)

# Concetrated Load
for x in xrange(1,11):
	val = 1000*x
	plt.title('Amplitude of oscillation against loading time - task 3')
	plt.ylabel('Amplitude (mm)')
	plt.xlabel('Load time (s)')
	plt.grid(b=True, which='major', color='k', linestyle='--')
	plt.grid(b=True, which='minor', color='gray', linestyle='--')
	crs = open('./output/data/task3_deflection' + str(val) + '.txt','r')
	D = np.loadtxt(crs)
	D2 = np.zeros(10000)
	# Look just at fully loaded area...
	for i in xrange(10000,len(D)):
		D2[i-10000] = D[i]

	amp[x] = max(D2) - min(D2)
	x2[x] = val/1000

plt.plot(x2,amp,'-ro',label=val)
plt.yscale('log')
plt.xscale('log')
plt.savefig('./output/images/task3_amplitude',bbox_inches='tight')