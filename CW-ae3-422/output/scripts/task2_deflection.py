#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
# 	Plots Analytical Solutions for AE3-422
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

crs = open('./output/data/task2_deflection.txt','r')
D = np.loadtxt(crs)

size = len(D)
# Plot Amplitudes
amp = np.zeros(size)
x2 = np.zeros(size)

plt.title('Amplitude of oscillation against loading time - task 2')
plt.ylabel('Amplitude (mm)')
plt.xlabel('Load time (s)')
plt.grid(b=True, which='major', color='k', linestyle='--')
plt.grid(b=True, which='minor', color='gray', linestyle='--')

for x in xrange(1,size):
	val = .001*x
	D2 = np.zeros(10000)
	# Look just at fully loaded area...
	for i in xrange(10000,20000):
		D2[i-10000] = D[x,i]

	amp[x] = (max(D2) - min(D2))/2
	mean = (max(D2) + min(D2))
	x2[x] = val
	plt.plot(x2,amp,'-',label=val)

plt.yscale('log')
plt.xscale('log')
plt.savefig('./output/images/task2_amplitude',bbox_inches='tight')

plt.clf()

# # Concetrated Load
for x in [2,4,6,8,10]:
	val = x*100
	plt.title('Mid-point displacement task 2')
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.ylabel('Deflection (mm)')
	plt.xlabel('Time (s)')
	plt.grid(b=True, which='major', color='k', linestyle='--')
	D2 = np.zeros(15000)
	for i in xrange(1,15000):
		D2[i] = D[val-1,i]
	x2 = np.linspace(0,1.5,15000) 
	plt.plot(x2,D2,label=val)
	plt.legend(bbox_to_anchor=(1, 0), loc=4, borderaxespad=0., title="Solution")

plt.axhline(y=10.75204443, xmin=0., xmax = 2., linewidth=1.5, color='k', ls='-')
plt.xlim([0,1.5])
plt.savefig('./output/images/task2_deflection',bbox_inches='tight')
plt.clf()