#----------------------------------------------------------#
#	Friedrich M. Grabner - 01220997
#	Title: File Saver
#	Description: 
#		Saves data from domain array.
#----------------------------------------------------------#
import numpy as np

def save(A ,nel):
	np.savetxt('./data/array' + str(nel-1) + '.txt', A, delimiter=',')
