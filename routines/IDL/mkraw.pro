PRO mkraw

;this program will make the raw magnitude data files

;if you use this code, please cite Oelkers & Stassun 2018

;;;;UPDATE INFORMATION HERE;;;;;
;;;be sure to update the number of light curves on line 40, now for simplicity, it is set to 100
;useful directories
cdedir = '../code/phot/' ;code directory
difdir = '../dif/' ; directory to put the differenced images
caldir = '../code/master/fin/' ; directory with master frame information
lcdir  = '../lc/' ; directory to put the light curves

;CCD information
camera = '2'
ccd = '2'
;;;END UPDATE INFORMATION;;;

;change into light curve directory and get the flux files
cd, difdir
spawn, 'ls *.flux', files
nfiles = n_elements(files)
cd, cdedir

;read in the master frame information to get the star count
readcol, caldir+camera+'_'+ccd+'_master_idl.ap', id, xm, ym, ttmg, mstmg, mstmge, format = '(l, f, f, f, f, f)', /silent, count = nstars
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
openw, 2, lcdir+'starinfo.txt'
for ii =0l, nstars-1 do begin
	nm=strcompress(id[ii],/remove_all)+'.lc'
	openw, 1, lcdir+nm
	for jj =0l, nfiles -1 do printf, 1, jd[jj], mg[jj,ii], er[jj,ii], format = '(f11.6, 1x, f9.6, 1x, f8.6)'
	close, 1
	printf, 2, nm, xm[ii], ym[ii], mstmg[ii], mstmge[ii], format ='(a15, 1x, f8.3, 1x, f8.3, 1x, f9.6, 1x, f8.6)'
endfor
close, 2
print, 'All done at '+systime()+'. See ya later alligator.'
END
