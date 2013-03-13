@../my_model/ice_model
@sft_class
@~/../arash/pro/textoidl/textoidl

;if we are serious about doing this report, we ought to start with
;the complex permittivity:
;use our dark nilas profile without snow:
pro write_dn_profile, t, s, dbub=dbub, dgr=dgr, bratio=bratio, bangle=bangle, $
		filename=filename, thick=thick
  if n_elements(dbub) ne 1 then dbub=1.2
  if n_elements(dgr) ne 1 then dgr=10.
  if n_elements(bratio) ne 1 then bratio=1.
  if n_elements(bangle) ne 1 then bangle=0.
  if n_elements(thick) ne 1 then thick=0.01

  if n_elements(filename) ne 1 then filename="dum.dat"

  openw, lun, filename, /get_lun
  printf, lun, 7
  printf, lun, 88.065, 0.17132
  printf, lun, 80.248, 0.36076
  printf, lun, 67.624, 0.46790
  printf, lun, 51.735, 0.46790
  printf, lun, 33.840, 0.36076
  printf, lun, 14.933, 0.17132
  printf, lun, 0.0, 0.0

  printf, lun, '3,10,50,F,3,T'
  printf, lun, -9., 271.4, 1., 0., 1., 0., 0., 32., 0., 0., 1
  printf, lun, -thick, t, 0.8, 0., 0., dgr, dbub, s, bangle, bratio, 2
  printf, lun, 0., t, 0.8, 0., 0., dgr, dbub, s, bangle, bratio, 2

  free_lun, lun

end


