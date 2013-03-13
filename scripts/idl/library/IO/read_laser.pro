@time_lib

;read the laser altimetry data:
pro read_laser, fname, z, lon, lat, t, s

  ;extract the date from the file name:
  files=strsplit(fname, "/", /extract)
  date0str=(strsplit(files[n_elements(files)-1], ".", /extract))[0]

  date0={time_str}
  date0.year=fix(strmid(date0str, 0, 4))
  date0.month=fix(strmid(date0str, 4, 2))
  date0.day=fix(strmid(date0str, 6, 2))
  date0.hour=fix(strmid(date0str, 8, 2))
  date0.minute=fix(strmid(date0str, 10, 2))

  openr, lun, fname, /get_lun
  junk=""
  readf, lun, junk
  nmax=1000000
  data=fltarr(6, nmax)

  on_ioerror, finish
  readf, lun, data

  finish:
  ind=where(data[0, *] ne 0, cnt)
  if cnt eq 0 then begin
    message, "read_laser: error, no data found in "+fname
    stop
  endif

  lon=reform(data[2, ind])
  lat=reform(data[1, ind])
  z=reform(data[5, ind])
  s=reform(data[4, ind])

  t=replicate({time_str}, cnt)

  tval=reform(data[3, ind])

  t[*].year=date0.year
  t[*].month=date0.month
  t[*].day=date0.day
  t[*].hour=tval
  t[*].second=(tval*3600) mod 60
  t[*].minute=(tval*60) mod 60

  free_lun, lun

  ;stop

end
  

