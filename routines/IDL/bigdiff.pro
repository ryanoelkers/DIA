PRO bigdiff

;this program will run the difference image analysis

;if you use this code, please cite: Alard & Lupton 1998, Alard 2000, Miller+2018, Oelkers+2015 & Oelkers & Stassun 2018

;;;;UPDATE INFORMATION HERE;;;;;

;useful directories
cdedir = '../code/diff/' ;code directory
caldir = '../code/master/fin/' ; directory for the location of the master frame
clndir = '../clean/'; directory of where the images are located
difdir = '../dif/' ; directory to put the differenced images

;this step compiles the C code and may need to be tweaked for directory information
spawn, 'gcc oisdifference.c -L/usr/local/lib -I/usr/local/include -lcfitsio -lm -lcurl'

;change this informaiton to update the differencing code kernel information
stmp = 3 ; what is the size of the stamp you want? 2xstmp+1
krn = 2 ; what is the size of the kernel you want? 2xkrn+1
spc = 0 ; what order kernel do you want?

;set up the aperture photometry information
apr = 2.5
tor = [apr+2,apr+4]

camera = '2'
ccd = '2'
;;;END UPDATE INFORMATION;;;;

;get the master frame and photometry
ref = readfits(caldir+camera+'_'+ccd+'_master_idl.fits', rhead,/silent)
readcol, caldir+camera+'_'+ccd+'_master_idl.flux', id, xm, ym, mflux, mfluxe, format = '(l, f, f, f, f)', /silent
mexpt = sxpar(rhead, 'EXPTIME')*3600.*24.

;move the reference stars
spawn, 'cp '+caldir+camera+'_'+ccd+'_starlist_idl.txt ./ref_stars.txt'

;read in the reference stars
readcol, 'ref_stars.txt', idx, xr, yr, format = '(i,f,f)', /silent, count = nrstars

;get the reference sky and subtract
sky, ref, skym, skys, /silent
nref = ref-skym
sze = size(ref)

;prepare files for the C code
openw, 1, 'ref.txt'
printf, 1, 'ref.fits'
close, 1
openw, 1, 'img.txt'
printf, 1, 'img.fits'
close, 1

;write out the reference file
writefits, 'ref.fits', nref, rhead

;get the image list
cd, clndir
spawn, 'ls *.fits ', files
nfiles = n_elements(files)
cd, cdedir

for ii = 0L, nfiles-1 do begin
    	;get the name of the file
    	nm = strsplit(files[ii], '.', /extract)
	nmfin = nm[0]+'d.fits'
	nm2 = strsplit(nmfin, '.', /extract)

	;only difference files that aren't yet differenced
    	if file_test(difdir+nmfin) ne 1 then begin

        	;read in the image file
        	img = readfits(clndir+files[ii], header, /silent)
        	sky, img, bk, bks, /silent
        	nimg = img-bk

        	;get the centroid of the reference stars
        	aper, nimg, xr, yr, mg, meg, sk, sks, 1., 5, [7,9], [-500,45000], /silent, setskyval = skym

        	;recover the reference stars
        	xx = xr
        	yy = yr
		bd1 = fltarr(n_elements(xx))

        	;remove bad reference stars
		for jj = 0l, n_elements(xx)-1 do begin
			dif = sqrt((xx[jj]-xx)^2+(yy[jj]-yy)^2)
			vv = where(dif lt 3. and (mg[jj]-mg) gt 0)
			if vv[0] ne -1 then bd1[jj] = 1
 		endfor
		;remove stars with high error and are too crowded
        	bd2 = where(mg lt 0 or mg gt 99. or meg lt 0 or meg gt .05 or bd1 eq 1)
        	if bd2[0] ne -1 and n_elements(bd2) ne n_elements(xx) then remove, bd2, xx, yy, mg, meg

		sss = sort(mg)
		;write the refernce stars for the C code, limit to 500 brightest
		nst = n_elements(mg)
		if (nst gt 500) then nst = 500
        	openw, 1, './refstars.txt'
        	for kkk = 0L, nst-1 do begin
            		printf, 1, xx[sss[kkk]],yy[sss[kkk]], format = '(i4,1x,i4)'
        	endfor
        	close, 1

		; assuming there are more stars to use than bad stars, begin the difference
        	if nst gt 1 then begin
            		;write the image out
            		writefits, 'img.fits', nimg, header

            		;print the parameters file, fwhm, w, d and number of stars
            		openw, 1, 'parms.txt'
            		printf, 1, stmp, krn, spc, nst, format = '(i1,1x,i1,1x,i1,1x,i3)'
            		close, 1
            
            		print, 'Now beginning difference imaging at ', systime(),'.'

            		;difference the images
            		spawn, './a.out'
						
            		;read in the differenced image
            		sxaddpar, header, "DIFF", 'yes'
            		dimg = readfits('./dimg.fits', head,/silent)
            		writefits, difdir+nmfin, float(dimg), header

            		print, 'File '+files[ii]+' now differenced at '+systime()+'.'

			;do the photometry on the differenced file
			expt = sxpar(header, 'EXPOSURE')*3600.*24.
			jd1 = sxpar(header, 'TSTART')
			jd2 = sxpar(header, 'TSTOP')
			jd = mean([jd1, jd2])

			sky, dimg, bkd, bkds, /silent
			aper, dimg, xm, ym, flux, fluxer, fsky, fskyer, 1, apr, tor, [-45000, 45000], /silent, /flux, setskyval = bkd
			mg = 25.-2.5*alog10(flux/expt+mflux/mexpt)
			mge = ((2.5)/(alog(10.))*(sqrt((mfluxe/mexpt)^2+(fluxer/expt)^2.)/(flux/expt+mflux/mexpt)))

			;output the flux file
			openw, 1, difdir+nm2[0]+'.flux'
			for kk = 0l, n_elements(xm)-1 do printf, 1, id[kk], xm[kk], ym[kk], jd, mg[kk], mge[kk], format = '(i0, 1x, f8.3, 1x, f8.3, 1x, f11.6, 1x,f9.6, 1x, f8.6)'
			close, 1
        	endif
    		if n_elements(bd) eq n_elements(xx) then print, 'File '+file[ii]+' is a bad file...'
    	endif
endfor

print, 'All done at '+systime()+'. See ya later alligator.'

END
