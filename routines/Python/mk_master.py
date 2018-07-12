#this program will combine images to make a master frame

#if you use this code, please cite Oelkers & Stassun 2018

#import the relevant libraries for basic tools
import pyfits
import numpy as np
import scipy
from scipy import stats
from os import path
import math
import time

#import relevant libraries for a list
import glob, os
from os import listdir
from os.path import isfile, join, exists

###UPDATE HERE#####
#what field are you looking at?
camera = '2'
ccd = '2'
blknum = 50 #how many images go into each holder?

#useful directories
cdedir = '../code/master/' #code directory
clndir = '../clean/' #directory where the cleaned images reside

###END UPDATE###

#get the image list and the number of files which need reduction
files = np.loadtxt(cdedir+'images.dat', dtype = 'string', unpack = 1, usecols = (0))
nfiles = len(files)

#iterate through the files
cnt = 0 #counter for the number of images used
kk = 0 #the jumper for file placement

for ii in range(0,len(files)):

	if (ii == 0):  # get size on first iteration only
		nx = pyfits.getval(clndir+files[0], 'NAXIS2')
	        ny = pyfits.getval(clndir+files[0], 'NAXIS1')
	        all_data = np.ndarray(shape=(blknum,nx,ny))
		expt = np.zeros(blknum)

	#read in the image
	img_data = pyfits.getdata(clndir+files[ii])
	expt[cnt] = pyfits.getval(clndir+files[ii],'EXPOSURE')

	#add the image to the vector
	all_data[cnt] = img_data 
	cnt = cnt+1

	if (ii % 10 == 0) and (ii > 0):
		print 'Finished with 10 images at '+str(time.strftime("%a %d %b %Y %H:%M:%S"))+'.'

	if (ii == len(files)-1) or ((ii+1) % blknum == 0):
		
		#median combine the data
		combined_data = np.median(all_data,axis=0)

		# Write data to new file    
		new_image = pyfits.PrimaryHDU(combined_data)
		new_image.header.set('NUMCOMB', cnt)
		new_image.header.set('EXPTIME', np.median(expt))

		#print the file with the appropriate counter
		if (kk < 10):
			new_image.writeto(cdedir+'frames/'+camera+'_'+ccd+'_master_0'+str(kk)+'.fits',clobber=True)
		else:
			if (kk >= 10) and (kk < 100):
				new_image.writeto(cdedir+'frames/'+camera+'_'+ccd+'_master_'+str(kk)+'.fits',clobber=True)

		print "The master frame hold was created using a median of "+str(cnt)+" images."
		kk = kk+1
		cnt = 0

		#clear the data file
	        all_data = np.ndarray(shape=(blknum,nx,ny))
		expt = np.zeros(blknum)

del all_data, img_data # clear up some memory
