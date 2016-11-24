#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
#	Title: File Saver
#	Description: 
#		Saves data from domain array.
#----------------------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt

def save(A ,nel):
	np.savetxt('./data/array' + str(nel-1) + '.txt', A, delimiter=',')

def plot_contour(x_mat,n):
	# Create Subplot
	ax = plt.subplot(111)
	ax.set_title('Temperature (K)')
	# colormap will be used for the contour lines
	im = plt.imshow(x_mat,interpolation='bilinear', origin='lower',extent=(0,1,0,1))
	levels = np.arange(-.1,2.1, 0.1)
	CS = plt.contour(x_mat, levels,origin='lower',linewidths=.5,extent=(0,1,0,1),colors='gray')
	plt.clabel(CS, inline=1, fontsize=10)
	# We can still add a colorbar for the image, too.
	CBI = plt.colorbar(im,orientation='horizontal', shrink=0.8)
	plt.xlabel('x axis')
	plt.ylabel('y axis')
	# Save figure...
	plt.savefig('./images/temp' + str(n) + '.svg')