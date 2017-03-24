#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
# 	Compares time to solution for task 4 - AE3-422
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

cnt = [24,48,72,96,144,192,288,480,984]
values = np.zeros((2,9))

for i in range(0,2):
	pnt = 0
	for j in cnt:
		crs = open('./output/data/timing_task4/np' + str(i+1) +'el' + str(j),'r')
		D =  np.loadtxt(crs)
		avg = np.mean(D)
		values[i,pnt] = avg
		pnt += 1
	plt.plot(cnt[:],values[i,:],label='Procs: '+str(i+1), lw=2.5)

plt.title('Time to solution for Np =  1 and 2')
plt.xlabel('Number of elements')
plt.ylabel('Time (s)')
plt.grid(True)
plt.legend(bbox_to_anchor=(1, 0), loc=4, borderaxespad=0., title="Num. of proc.")
plt.savefig('./output/images/task4_timing',bbox_inches='tight')