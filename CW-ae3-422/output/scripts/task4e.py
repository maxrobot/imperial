#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
# 	Compares time to solution for task 4 - AE3-422
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

cnt = [24,48,96,144,192,288]
values = np.zeros((2,6))

for i in range(0,2):
	pnt = 0
	for j in cnt:
		crs = open('./output/data/timing_task4/np' + str(i+1) +'el' + str(j),'r')
		D =  np.loadtxt(crs)
		avg = np.mean(D, axis=0)
		values[i,pnt] = avg[0]
		pnt += 1
	plt.plot(cnt[:],values[i,:],label='Procs: '+str(i+1))

plt.title('Time to solution for Np =  1 and 2')
plt.xlabel('Number of elements')
plt.ylabel('Time (s)')
plt.legend(bbox_to_anchor=(1, 0), loc=4, borderaxespad=0., title="Num. of proc.")
plt.savefig('./output/images/task4_timing',bbox_inches='tight')



# # Concetrated Load
# for x in xrange(1,6):
# 	val = 2000*x
# 	print(val)
# 	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# 	plt.ylabel('Deflection (mm)')
# 	plt.xlabel('Time (s)')
# 	plt.grid(b=True, which='major', color='k', linestyle='--')
# 	crs = open('./output/data/task3_deflection' + str(val) + '.txt','r')
# 	D = np.loadtxt(crs)
# 	x2 = np.linspace(0,2,len(D)) 
# 	plt.plot(x2,D,label=val)
