;reads SSM/I EASE-grid data:

function get_ease, date, hemi, ch89=ch89, basepath=basepath, $
		xgrid=xgrid, ygrid=ygrid, lon=lon, lat=lat

  if n_elements(basepath) ne 1 then basepath="/rinax/pmills/smosdata/ssmiease/"

  monthname=["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]

  t=convert_date(date)
  month=monthname[t.month-1]
  yr=t.year mod 100

  fname=string(t.year, format="(i4.4)")+"/"+month+"/"+$
	  	string(yr, t.month, t.day, format="(i2.2, i2.2, i2.2)")

  if hemi lt 0 then begin
    if keyword_set(ch89) then begin
      nx=632
      ny=664
      nchan=2
      dsg="s3a"
    endif else begin
      nx=316
      ny=332
      nchan=5
      dsg="s3b"
    endelse
    x0=-3950.
    xf=3950.
    y0=-3950.
    yf=4350.
  endif else begin
    if keyword_set(ch89) then begin
      nx=608
      ny=896
      nchan=2
      dsg="n3a"
    endif else begin
      nx=304
      ny=448
      nchan=5
      dsg="n3b"
    endelse
    x0=-3850.
    xf=3750.
    y0=-5350.
    yf=5850.
  endelse

  raw=uintarr(nchan, nx, ny)

  openr, lun, basepath+dsg+"/"+fname+"."+dsg, /get_lun
  readu, lun, raw

  free_lun, lun

  data=0.1*raw+50.

  xgrid=x0+findgen(nx)*(xf-x0)/(nx-1)
  ygrid=y0+findgen(ny)*(yf-y0)/(ny-1)

  if arg_present(lon) or arg_present(lat) then begin
    xarr=rebin(xgrid, nx, ny)
    yarr=rebin(transpose(ygrid), nx, ny)
    mapxy, xarr, yarr, lat, lon, hemi
  endif

  return, data
end

