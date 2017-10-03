PRO refphot

;this program will do the photometry on the reference frame to get stars for the subtraction
;and to get stars for later

;if you use this code, please cite Oelkers et al. 2015, AJ, 149, 50

;;;;UPDATE INFORMATION HERE;;;;;

;useful directories
cdedir = '../IDL/' ;code directory
caldir = '../calib/' ; directory to output the flux files and location of master frame

;update the aperture size based on the images
apr = 3
tor = [apr+2,apr+4]
;;;END UPDATE INFORMATION HERE;;;;

ref = readfits(caldir+'master.fits', /silent)

print, 'Getting flux information from master frame at '+systime()+'.'
;find the sky information
sky, ref, skym, skys, /silent

;find the stars in the frame
find, ref, x,y, flux, sharp, round, 3*skys, 1.59, [-1.,1.], [.2,1.], /silent

;do the aperture photometry on those stars in both flux and magntiude space
aper, ref, x, y, mags, errap, sky, skyer, 1, apr, tor, [-500, 45000], /silent
aper, ref, x, y, flux, fluxer, fsky, fskyer, 1, apr, tor, [-500, 45000], /silent, /flux

;get rid of bad photomerty
good = where(mags lt 30. and mags gt 0. and errap gt 0 and errap lt 0.25)

;get rid of possible bad stars
x = x[good] & y = y[good] & flux = flux[good] & fluxer = fluxer[good] & sky = sky[good]
skyer = skyer[good]& fsky = fsky[good] & fskyer = fskyer[good] & mags = mags[good] & errap = errap[good]

;write magnitude information
openw, 1, caldir+'master.ap'
for i = 0L, n_elements(x) -1 do printf, 1, i, x[i], y[i], mags[i], errap[i],sky[i], skyer[i], format = '(i0, 1x, d0,1x, d0,1x,d0,1x,d0,1x,d0, 1x, d0)'
close, 1
;write the flux information
openw, 1, caldir+'master.flux'
for i = 0L, n_elements(x) -1 do printf, 1, i, x[i], y[i], flux[i], fluxer[i], fsky[i],fskyer[i], format = '(i0, 1x, d0,1x, d0,1x,d0,1x,d0,1x,d0, 1x, d0)'
close, 1
;write the star list of differencing 
openw, 1, caldir+'starlist.txt'
for i = 0L, n_elements(x) -1 do printf, 1, i, x[i], y[i], format = '(i0,1x,d0,1x,d0)'
close, 1

print, 'Finished at '+systime()+'. See ya later alligator!'

END
