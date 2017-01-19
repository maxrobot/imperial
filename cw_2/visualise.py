#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
#----------------------------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


for x in xrange(1,5):
	# Open the first file
	crs = open('./data/output000' + str(x),'r')
	D = np.loadtxt(crs)

	# Create Subplot
	ax = plt.subplot(111)
	ax.set_title('Vorticity')
	# quadmesh = ax.pcolormesh(X,Y,data)
	# colormap will be used for the contour lines
	im = plt.imshow(D,interpolation='bilinear', origin='lower',extent=(0,1,0,1))
	im.set_clim(vmin=-.25, vmax=.25)

	plt.xlabel('x axis')
	plt.ylabel('y axis')
	if x==1:
		CBI = plt.colorbar(im,orientation='horizontal', shrink=0.8)

	# We can still add a colorbar for the image, too.
	# Save figure...
	plt.savefig('./images/temp000' + str(x))