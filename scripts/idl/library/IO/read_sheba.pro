@time_lib

pro read_sheba, filename, temp, z, t
  maxn=20000

  openr, lun, filename, /get_lun
  head=""
  readf, lun, head
  field=strsplit(head, ",", /extract)
  nz=n_elements(field)-1
  z=float(field[1:nz])

  data=dblarr(nz+1, maxn)

  on_ioerror, finish
  readf, lun, data

  finish:
    ind=where(data[0, *] ne 0, n)
    free_lun, lun

  temp=float(data[1:nz, ind])
  t=replicate({time_str}, n)

  traw=reform(data[0, ind])
  t.second=float(traw mod 100)
  traw=traw/100
  t.minute=fix(traw mod 100)
  traw=traw/100
  t.hour=fix(traw mod 100)
  traw=traw/100
  t.day=fix(traw mod 100)
  traw=traw/100
  t.month=fix(traw mod 100)
  t.year=traw/100

end
