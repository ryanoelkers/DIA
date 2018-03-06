#this program will do the photometry on the reference frame to get stars 
#for the subtraction and to get stars for later

#if you use this code, please cite Oelkers & Stassun 2018

#import the relevant libraries for basic tools
import numpy
import scipy
from scipy import stats
import scipy.ndimage as ndimage
import astropy
from astropy.stats import sigma_clipped_stats
import math
import time
from photutils import DAOStarFinder
from photutils import aperture_photometry
from photutils import CircularAperture
from photutils import CircularAnnulus

#for reading in fits files
from astropy.io import fits
from astropy.wcs import WCS

#import relevant libraries for a list
import glob, os
from os import listdir
from os.path import isfile, join, exists

###BEGIN UPDATE#####
camera = '2'
ccd = '2'

#useful directories
cdedir = '../code/master/' #code directory
caldir = '../code/master/fin/' # directory to output the flux files and location of master frame
clndir = '../clean/'#directory with the cleaned images
##END UPDATE####

#read in the master frame
mast, mheader = fits.getdata(caldir+camera+'_'+ccd+'_master.fits', header = True)

###POSSIBLE TO UNCOMMENT HERE IF YOU HAVE NO STAR LIST###

#get the positions in the master, if no star list is provided 
#mean, median, std = sigma_clipped_stats(mast, sigma = 3.0, iters = 5)
#daofind = DAOStarFinder(fwhm = 1.5, threshold = 5.*std)
#sources = daofind(mast-median)
#set up for aperture photometry on the master frame
#mean, median, std = sigma_clipped_stats(mast, sigma = 3.0, iters = 5)
#positions = (sources['xcentroid'], sources['ycentroid'])
#x = sources['xcentroid']
#y = sources['ycentroid']

###END POSSIBLE UNCOMMENT###

###POSSIBLE TO UNCOMMENT HERE IF YOU HAVE A STAR LIST###
#get the positions from the star list if one is provided ###COMMENT HERE###
ticid, tmag, ra, dec = numpy.loadtxt('ffi_test.csv', unpack = 1, delimiter = ',', skiprows =1)

#convert the ra/dec into positions
os.chdir(rawdir) #changes to the raw image direcotory
files = [f for f in glob.glob("*.fits") if isfile(join(rawdir, f))] #gets the relevant files with the proper extension
files.sort()
os.chdir(cdedir) #changes back to the code directory

#get wcs from first calibrated file
w = WCS(clndir+files[0])
x, y = w.all_world2pix(ra, dec, 0)
positions = (x, y)
###END POSSIBLE UNCOMMENT###

print 'Getting the photometry on the master frame.'
#apertures to test for optimal apeture size
rads = numpy.arange(2,5,.25) 

#do the aperture photometry and find the optimal aperture
apertures = [CircularAperture(positions, r=r) for r in rads]
phot_table = aperture_photometry(mast, apertures, method = 'exact')
idx = 0

offset = numpy.zeros((len(rads),len(x)))
for ii in range(0, len(x)):
	if (x[ii] > 0) and (x[ii] < 2048) and (y[ii] > 0) and (y[ii] < 2048):
		dist = numpy.sqrt((x[ii]-x)**2+(y[ii]-y)**2)
		chk = numpy.where(dist < 6.)
		if (len(chk[0]) == 1):
			for jj in range(1, len(rads)):
				mg1 = 25.-2.5*numpy.log10(phot_table[ii][jj+3])
				mg0 = 25.-2.5*numpy.log10(phot_table[ii][jj+2])
				offset[jj,ii] = mg1-mg0
prv = 1.
opt_rad = 10.
for ii in range(0, len(rads)):
	chk = numpy.median(offset[ii,:])	
	if (numpy.abs(chk-prv) <= 0.001) and (rads[ii] < opt_rad):	
		opt_rad = rads[ii]
		print 'The optimal aperture size is '+str(opt_rad)+'.'
	if (numpy.abs(chk-prv) > 0.001):
		prv = chk

#do the aperture photometry
apertures = CircularAperture(positions, r = opt_rad)
phot_table = aperture_photometry(mast, apertures, method = 'exact')

#get the background of the image
cimg, clow, chigh = scipy.stats.sigmaclip(mast, low=2.5, high = 2.5) #do a 2.5 sigma clipping
bkg_mean = numpy.median(cimg) #determine the sky value
sig = numpy.std(cimg) #determine the sigma(sky)

#convert to magnitudes
flx = phot_table['aperture_sum']-(bkg_mean*(numpy.pi*opt_rad**2))
flx_er = numpy.sqrt(phot_table['aperture_sum'])
x_pix = x
y_pix = y

#create the magnitudes from the flux
mag = 25.0-2.5*numpy.log10(flx)
err = (2.5/numpy.log(10.))*(flx_er/flx)

#write the magnitudes to a file
output = open(caldir+camera+'_'+ccd+'_master.ap', 'w')
for ii in range(0, len(phot_table['id'])):
	if (x_pix[ii] > 0) and (x_pix[ii] < 2048) and (y_pix[ii] > 0) and (y_pix[ii] < 2048) and (numpy.isnan(mag[ii]) != 1):
		output.write(str(long(ticid[ii]))+','+str(x_pix[ii])+','+str(y_pix[ii])+','+str(tmag[ii])+','+str(mag[ii])+','+str(err[ii])+'\n')
output.close()

#write the fluxes to a file
output = open(caldir+camera+'_'+ccd+'_master.flux', 'w')
for ii in range(0, len(phot_table['id'])):
	if (x_pix[ii] > 0) and (x_pix[ii] < 2048) and (y_pix[ii] > 0) and (y_pix[ii] < 2048) and (numpy.isnan(mag[ii]) != 1):
		output.write(str(long(ticid[ii]))+','+str(x_pix[ii])+','+str(y_pix[ii])+','+str(flx[ii])+','+str(flx_er[ii])+'\n')
output.close()

#write the star list to a file
output = open(caldir+camera+'_'+ccd+'_starlist.txt', 'w')
for ii in range(0, len(phot_table['id'])):
	if (x_pix[ii] > 0) and (x_pix[ii] < 2048) and (y_pix[ii] > 0) and (y_pix[ii] < 2048) and (numpy.isnan(mag[ii]) != 1):
		output.write(str(long(ticid[ii]))+','+str(long(x_pix[ii]))+','+str(long(y_pix[ii]))+'\n')
output.close()
