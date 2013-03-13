@time_lib

pro read_thickness, tfile, thickness, lon, lat, t, flag, $
		remove_bad=remove_bad

  maxn=100000
  data1=fltarr(10, maxn)
  on_ioerror, finish1
  openr, lun, tfile, /get_lun
  readf, lun, data1

  finish1: ind=where(data1[0, *] ne 0, n1)
  data1=data1[*, ind]
  free_lun, lun

  lon=reform(data1[6, *])
  lat=reform(data1[7, *])
  flag=reform(fix(data1[9, *]))

  t=replicate({time_str}, n1)

  t.year=reform(data1[0, *])
  t.month=reform(data1[1, *])
  t.day=reform(data1[2, *])
  t.hour=reform(data1[3, *])
  t.minute=reform(data1[4, *])
  t.second=reform(data1[5, *])

  thickness=reform(data1[8, *])

  if keyword_set(remove_bad) then begin
    ind=where(flag eq 0)
    lon=lon[ind]
    lat=lat[ind]
    t=t[ind]
    thickness=thickness[ind]
    flag=flag[ind]
  endif

  ;stop

end

