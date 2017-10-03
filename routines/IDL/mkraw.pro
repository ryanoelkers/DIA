PRO mkraw

;this program will make the raw magnitude data files

;if you use this code, please cite Oelkers et al. 2015, AJ, 149, 50

;;;;UPDATE INFORMATION HERE;;;;;
;;;be sure to update the number of light curves on line 40, now for simplicity, it is set to 100
;useful directories
cdedir = '../IDL/' ;code directory
difdir = '../diff/' ; directory to put the differenced images
caldir = '../calib/' ; directory with master frame information
lcdir  = '../lightcurve/' ; directory to put the light curves

;;;END UPDATE INFORMATION;;;

;change into light curve directory and get the flux files
cd, difdir
spawn, 'ls *.flux', files
nfiles = n_elements(files)
cd, cdedir

;read in the master frame information to get the star count
readcol, caldir+'master.flux', id, xm, ym, mflux, mfluxe, format = '(l, f, f, f, f)', /silent, count = nstars
jd = fltarr(nfiles)
mg = fltarr(nfiles, nstars)-99.00
er = fltarr(nfiles, nstars)-99.00

print, 'Splitting photometry files at '+systime()+'.'
;loop through each file to combine the files
for ii = 0L, nfiles-1 do begin
	readcol, difdir+files[ii], idd, mx, my, jdh, mgh, mgeh, format = '(l, f, f, f, f, f)', /silent
	match, idd, id, aa, bb
	jd[ii] = jdh[0]
	mg[ii, bb] = mgh[aa]
	er[ii, bb] = mgeh[aa]
endfor

;loop through the stars make the light curves right now it is set to 100 to simiplify the code
for ii =0l, 100 do begin
	if ii lt 10 then nm = '0000'+strcompress(ii, /remove_all)+'.lc'
	if ii lt 100 and ii ge 10 then nm = '000'+strcompress(ii, /remove_all)+'.lc'
	if ii lt 1000 and ii ge 100 then nm = '00'+strcompress(ii, /remove_all)+'.lc'
	if ii lt 10000 and ii ge 1000 then nm = '00'+strcompress(ii, /remove_all)+'.lc'
	if ii lt 100000 and ii ge 10000 then nm = '0'+strcompress(ii, /remove_all)+'.lc'
	if ii ge 100000 then nm = strcompress(ii, /remove_all)+'.lc'

	openw, 1, lcdir+nm
	for jj =0l, nfiles -1 do printf, 1, jd[jj], mg[jj,ii], er[jj,ii], format = '(f11.6, 1x, f6.3, 1x, f5.3)'
	close, 1
endfor
print, 'All done at '+systime()+'. See ya later alligator.'
END
