#this program will detrend stars based on similar magnitude

#if you use this code, please cite Oelkers & Stassun 2018 submitted

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

###UPDATE INFORMATION HERE
camera = '2'
ccd = '2'

#useful directories
lcdir = '../lc/' #the directory where the raw light curves 'live'
otdir = '../lc/detrend/' #the directory where the detrended light curves will be output to
cdedir = '../code/phot/' #directory where this code 'lives'
mstdir = '../code/master/fin/' #directory with the master frame

###END UPDATE INFORMATION

#read in the star information
nme, mg = np.loadtxt(mstdir+'2_2_master.ap', usecols = (0,4), unpack = 1)
nstars = len(nme) # get the total number of stars in the frame

#read in all light curves to decrease the detrending time
#make the light curve holder
big = np.zeros((nstars, 1348))
s = np.argsort(mg)

for ii in range(0, nstars):
	#read in the light curves based on magnitude
	m = np.loadtxt(lcdir+str(long(nme[s[ii]]))+'.lc', usecols = (1))
	big[ii,:] = m
	if (ii% 1000 == 0) and (ii > 0):
		print '1000 light curves read in at '+str(time.strftime("%a %d %b %Y %H:%M:%S"))+'.'


#now loop through the stars and build the trend
for ii in range(0, nstars):
	#read in the star to detrend
	j, m, e = np.loadtxt(lcdir+str(long(nme[s[ii]]))+'.lc', unpack = 1)
	print nme[s[ii]]
	#find the +/-500 stars to use in magnitude space
	udx = 1000-ii
	if (udx < 500):
		udx = ii+500
	if (udx > nstars-1):
		 udx = nstars-1
	ldx = ii-500
	if (ldx < 0):
		ldx = 0

	#initial check of rms values
	chk1 = np.std(m)

	#holder for the trend stars
	trd = np.zeros((1000,1348))

	#loop through and check for >10% improvement in rms
	for jj in range(ldx, udx): 
		if (jj != ii):
			#second check of rms values
			chk2 = np.std(m-(big[jj,:]-np.median(big[jj,:])))
			if (chk1/chk2 >= 1.1): #if it is better than 10%, include the star
				trd[jj-ldx,:] = big[jj,:]-np.median(big[jj,:])

	
	trend = np.zeros(1348) #make the holder for the final trend
	offs = np.zeros(1348) #make the holder to find the light curve offsets

	#now loop through and median combine the trends, if they exist
	for jj in range(0, 1348):
		vvv = np.where(trd[:,jj] != 0)
		vv = vvv[0]
		if (len(vv) != 0):
			trend[jj] = np.median(trd[vv,jj])
			if (jj > 0):
				offs[jj] = np.abs(trend[jj]-trend[jj-1])

	
	#check to see if a trend was found to remove
	chk3 = np.where(trend != 0)
	otrend = trend #keep the old trend just in case

	#if a trend was found, begin attempting to scale
	if len(chk3[0] != 0):
		#check for large deviations between subsequent data points
		cmg, clow, chigh = scipy.stats.sigmaclip(offs, low=3, high = 3) #do a 3  sigma clipping
		mm = np.median(cmg) #determine the mean magnitude
		ss = np.std(cmg) #determine the sigma(magnitude)

		brk = mm+5*ss#break occur at large offsets

		#get the indicies of the 'breaks'
		tmes = np.where(offs > brk)
		tme = np.concatenate(([0],tmes[0]), axis = 0)

		#if large devaitions exist, move through them and scale
		for jj in range(0, len(tme)):
	
			# get the lower and upper limits of the 'broken' light curve
			lw = tme[jj]
			if (jj < len(tme)-1):
				upp = tme[jj+1]-1
			if (jj == len(tme)-1):
				upp = 1347

			#determine the median offset of the light curve
			mm = np.median(m[lw:upp])-np.median(m)
			mt = np.median(trend[lw:upp])

			#make the offset
			off = mm-mt

			#check the std with and without the offset
			chk5 = np.std(m[lw:upp])
			chk4 = np.std(m[lw:upp] - (trend[lw:upp]+off))

			#updat the trend
			trend[lw:upp] = trend[lw:upp]+off

		#now we check to make sure the offset doesn't make the light curve more messy
		df1 = np.zeros(1348)
		df2 = np.zeros(1348)
		
		for jj in range(1,1348):
			df1[jj] = m[jj]-m[jj-1]
			df2[jj] = (m[jj]-trend[jj])-(m[jj-1]-trend[jj-1])

		#without the trend removed
		cdf1, clow, chigh = scipy.stats.sigmaclip(df1, low=2.5, high = 2.5) #do a 3 sigma clipping
		mm1 = np.median(cdf1) #determine the sky value
		ss1 = np.std(cdf1) #determine the sigma(sky)

		#with the trend removed
		cdf2, clow, chigh = scipy.stats.sigmaclip(df2, low=2.5, high = 2.5) #do a 3 sigma clipping
		mm2 = np.median(cdf2) #determine the sky value
		ss2 = np.std(cdf2) #determine the sigma(sky)

		#number of data points in large sigma regime
		num1_idx = np.where(np.abs(df1) > mm1+5*ss1)
		num1 = len(num1_idx[0])
		num2_idx = np.where(np.abs(df2) > mm2+5*ss2)
		num2 = len(num2_idx[0])

		#if more breaks are found in the light curve, keep the old trend
		if (num2 >= num1):
			trend = otrend

		#output the cleaned light curve
		output = oppen(otdir+nme[s[ii]]+'.lc', 'w')
		for jj in range(0, len(j)): 
			output.write('%f %f %f\n' % (np.round(j[jj], decimals = 6), np.round(m[jj]-trend[jj], decimals = 6), np.round(e[jj], decimals = 6)))
		output.close()
	
	#if no trend stars exist, the output the raw light curve
	if len(chk3[0] == 0): 
		output = oppen(otdir+nme[s[ii]]+'.lc', 'w')
		for jj in range(0, len(j)): 
			output.write('%f %f %f\n' % (np.round(j[jj], decimals = 6), np.round(m[jj], decimals = 6), np.round(e[jj], decimals = 6)))
		output.close()

	if (ii % 1000 == 0) and (ii > 0): 
		print 'Working on the next 1000 curves at '+str(time.strftime("%a %d %b %Y %H:%M:%S"))+'.'

print 'All done at '+str(time.strftime("%a %d %b %Y %H:%M:%S"))+'. See ya later alligator!'
