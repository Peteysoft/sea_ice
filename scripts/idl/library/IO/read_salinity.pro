@time_lib

pro read_salinity, filename, z1, z2, s, error=err, $
		date=date, lon=lon, lat=lat


  on_ioerror, finish
  openr, lun, filename, /get_lun
  
  ;read the headers:
  header=""
  readf, lun, header

  ;just lop off the date field so it doesn't get in our way
  ;add a dummy field so we don't have to updated the rest of the code...
  field=[strmid(header, 0, 16), strsplit(strmid(header, 16), /extract)]

  ;field=strsplit(header, /extract)
  datestr=strtrim(field[0], 2)

  date={time_str, strmid(datestr, 6, 4), strmid(datestr, 3, 2), strmid(datestr, 0, 2), $
		strmid(datestr, 10, 2), strmid(datestr, 13, 2), 0}

  lon=float(field[1])
  lonlon=long(lon)
  lon=lonlon+(lon-lonlon)*100./60.

  lat=float(field[1])
  lonlat=long(lat)
  lat=lonlat+(lat-lonlat)*100./60.

  ;print, time_string(date), " (", lon, ", ", lat, ")"

  ;read the comments:
  ncom=0L
  readf, lun, ncom
  for i=0, ncom-1 do begin
    readf, lun, header
    ;print, header
  endfor

  ltype=0L
  dtype=0L
  n=0L
  while dtype ne 2 do begin
    readf, lun, ltype, dtype, n
    ;print, ltype, dtype, n
    if ltype eq 2 then begin
      data=fltarr(3, n)
      readf, lun, data
    endif else begin
      data=fltarr(2, n)
      readf, lun, data
    endelse
  endwhile

  z1=data[0, *]
  z2=data[1, *]
  s=data[2, *]

  err=0
  finish: if n_elements(dtype) ne 0 then begin
    if dtype ne 2 then err=1
  endif else begin
    err=2
  endelse

  free_lun, lun

end
