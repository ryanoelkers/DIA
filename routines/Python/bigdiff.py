#this program will run the difference image analysis from Oelkers et al. 2015

#if you use this code, please cite Oelkers et al. 2015, AJ, 149, 50, Alard et al. 2000 and Miller et al. 2008

#import the relevant libraries for basic tools
import numpy
import scipy
from scipy import stats
import scipy.ndimage as ndimage
import astropy
from astropy.stats import sigma_clipped_stats
import math
import time
from time import strftime
from photutils import aperture_photometry
from photutils import CircularAperture
from photutils import CircularAnnulus
import random
import pyfits

#for reading in fits files
from astropy.io import fits

#import relevant libraries for a list
import glob, os
from os import listdir
from os.path import isfile, join, exists

###UPDATE HERE#####
#compile the C differencing program
compdiff = os.system('gcc oisdifference.c -L/usr/local/lib -I/usr/local/include -lcfitsio -lm')

#useful directories
cdedir = '../Python/' #code directory
caldir = '../calib/' # directory for the location of the master frame
clndir = '../clean/'# directory of where the images are located
difdir = '../diff/' # directory to put the differenced images

#the optimal aperture to use from get_opt_rad.py
rad = 3.0

#size of the kernel, stamp and if you want an order
krnl = 3
stmp = 5
ordr = 0
nrstars = 200

#read in the master frame
mlist = fits.open(caldir+'master.fits')
mheader = mlist[0].header #get the header info
mast = mlist[0].data #get the image info
mean, median, std = sigma_clipped_stats(mast, sigma = 3.0, iters = 5)
expm_time = pyfits.getval(caldir+'master.fits','EXPTIME')

#subtract the background from the master frame
nmast = mast-median

#write the new master file
mhd = fits.PrimaryHDU(nmast, header=mheader)
shush = os.system('rm '+cdedir+'ref.fits')
mhd.writeto(cdedir+'ref.fits')
mlist.close()

#get the image list and the number of files which need reduction
os.chdir(clndir) #changes to the raw image direcotory
files = [f for f in glob.glob("*.fits") if isfile(join(clndir, f))] #gets the relevant files with the proper extension
nfiles = len(files)
os.chdir(cdedir) #changes back to the code directory

#read in the star list
ids, xx, yy = numpy.loadtxt(caldir+'starlist.txt', unpack = 1, delimiter = ',')
ids1, xm, ym, mflx, mflx_er = numpy.loadtxt(caldir+'master.flux', unpack = 1, delimiter = ',')

#begin with the algorithm to difference the images
for ii in range(0, nfiles):
	hld = files[ii].split('_')
	finnme = hld[0]+'_d'+hld[1]
	#check to see if the differenced file already exists
	if (os.path.isfile(difdir+finnme) == 0):
		#read in the image
		imglist = fits.open(clndir+files[ii])
		iheader = imglist[0].header #get the header info
		img = imglist[0].data #get the image info
		mean, median, std = sigma_clipped_stats(img, sigma = 3.0, iters = 5)
	
		#write the new image file
		nimg = img-median
		ihd = fits.PrimaryHDU(nimg, header=iheader)
		shush = os.system('rm '+cdedir+'img.fits')
		ihd.writeto(cdedir+'img.fits')
		exp_time = pyfits.getval(cdedir+'img.fits','EXPTIME')
		jd = pyfits.getval(cdedir+'img.fits','JD')-2454000.0
		print 'Getting magnitudes from the star list at '+strftime("%a, %d %b %Y %H:%M:%S")+'.'

		#determine the magnitudes, errors and distances to the objects
		#prepare the apertures
		positions = [xx,yy]
		apertures = CircularAperture(positions, r = rad)
		annulus_apertures = CircularAnnulus(positions, r_in=rad+2, r_out=rad+4)
		
		#get the photometry for the stars
		rawflux = aperture_photometry(img, apertures)

		#get the background
		bkgflux = aperture_photometry(img, annulus_apertures)		
		bkg_mean = bkgflux['aperture_sum']/annulus_apertures.area()
		bkg_sum = bkg_mean*apertures.area()

		#get the star flux and error & mag and error
		flx = rawflux['aperture_sum']-bkg_sum
		flx_er = numpy.sqrt(rawflux['aperture_sum'])
		mag = 25.-2.5*numpy.log10(flx)+2.5*numpy.log10(exp_time)
		mag_er = (2.5/numpy.log(10.))*(flx_er/flx)

		print 'Getting reference stars for the subtraction at '+strftime("%a, %d %b %Y %H:%M:%S")+'.'
		#star selecting stars that are not near any other bright stars	
		output = open(cdedir+'refstars.txt', 'w')
		cnt = 0
		itr = 0
		x = xx
		y = yy

		while (cnt < nrstars) and (itr < len(xx)):
			#select a random object
			jj = random.randint(0,len(x)-1)
			if (mag_er[jj] > 0) and (mag_er[jj] < 0.02):
				#get the nearest neightbors in 20 pix
				d = numpy.sqrt((x[jj]-x)**2+(y[jj]-y)**2)
				idx = numpy.where(d < 10)
				idxs = idx[0]
				#assuming the star is alone and it is not near an edge
				if (len(idxs) == 0) and (x[jj] < 3596) and (x[jj] > 500) and (y[jj] < 3596) and (y[jj] > 500):
					output.write(str(x[jj])+' '+str(y[jj]))
					cnt = cnt+1
				else:
					if (len(idxs) > 0) and (x[jj] < 3596) and (x[jj] > 500) and (y[jj] < 3596) and (y[jj] > 500):
						#check the magnitudes in case the neighbors are just faint stars
						dmag = mag[jj]-mag[idxs]
						cdmag = dmag[numpy.where(dmag != 0)]
						chk = numpy.where(cdmag > -1.5)
						if (len(chk[0]) == 0):
							output.write("%4d %4d\n" % (x[jj],y[jj]))
							cnt = cnt+1
			x = numpy.delete(x,jj)
			y = numpy.delete(y,jj)
			itr = itr+1
		nrstars = cnt #update the number of stars to use, just in case the maximum wasn't found but we ran out of stars
		output.close()
		imglist.close()

		#write the parameter file now that we have the stars
		output = open(cdedir+'parms.txt', 'w')
		output.write("%1d %1d %1d %4d\n" % (stmp, krnl, ordr, nrstars))
		output.close()

		output = open(cdedir+'ref.txt', 'w')
		output.write("ref.fits\n")
		output.close()

		output = open(cdedir+'img.txt', 'w')
		output.write("img.fits\n")
		output.close()


		#do the differencing!
		print 'Now starting the subtraction at '+strftime("%a, %d %b %Y %H:%M:%S")+'.'
		dodiff = os.system('./a.out')
		mvdiff = os.system('mv dimg.fits '+difdir+finnme)
		
		#get the photometry from the differenced image
		print 'Now starting the photometry at '+strftime("%a, %d %b %Y %H:%M:%S")+'.'

		#read in the image
		diflist = fits.open(difdir+finnme)
		iheader = diflist[0].header #get the header info
		dif = diflist[0].data #get the image info
	
		print 'Getting fluxes from the differenced file at '+strftime("%a, %d %b %Y %H:%M:%S")+'.'
		#determine the magnitudes, errors and distances to the objects
		#prepare the apertures
		positions = [xx,yy]
		apertures = CircularAperture(positions, r = rad)
		annulus_apertures = CircularAnnulus(positions, r_in=rad+2, r_out=rad+4)
		
		#get the photometry for the stars
		rawflux = aperture_photometry(dif, apertures)

		#get the background
		bkgflux = aperture_photometry(dif, annulus_apertures)		
		bkg_mean = bkgflux['aperture_sum']/annulus_apertures.area()
		bkg_sum = bkg_mean*apertures.area()

		#get the star flux and error & mag and error
		flx = rawflux['aperture_sum']-bkg_sum
		flx_er = numpy.sqrt(numpy.abs(rawflux['aperture_sum']))
		mag = 25.-2.5*numpy.log10(flx/exp_time+mflx/expm_time)
		mag_er = (2.5/numpy.log(10.))*(numpy.sqrt((flx_er/exp_time)**2+(mflx_er/expm_time)**2)/(flx/exp_time+mflx/expm_time))
		diflist.close()

		#print the flux information to the data file
		nme = finnme.split('.')
		output = open(difdir+nme[0]+'.flux', 'w')
		for jj in range(0, len(xx)):
			output.write(str(long(ids[jj]))+','+str(xm[jj])+','+str(ym[jj])+','+str(jd)+','+str(mag[jj])+','+str(mag_er[jj])+'\n')
		output.close()
		print 'Moving to the next file at '+strftime("%a, %d %b %Y %H:%M:%S")+'.'
print 'Done with the photometry at '+strftime("%a, %d %b %Y %H:%M:%S")+'.'
