;creates a synthetic ice profile based on four parameters:
;salinity profile type: C, S, ?, I
;thickness (number of layers)
;surface temperature
;average salinity

;writes it to a file

pro create_syn_prof, fname, stype, s_ave, temp, n, dd
  tsea=271.		;sea temperature
  if n_elements(dd) eq 0 then dd=0.05		;thickness of each layer
  x=findgen(n)*dd	;depth

  denst=0.9		;density at the top
  densb=0.9		;density at the bottom

  ;aspect_ratio=1/10.
  aspect_ratio=1.

  ;salinity profile:
  s=create_salinity_prof(x, s_ave, stype)

  dbub=fltarr(n)+1.2	;constant bubble diameter
  dgr=fltarr(n)+1.2	;constant grain diameter

  ;linear temperature profile:
  temp=findgen(n)*(tsea-temp)/(n-1)+temp

  ;linear density profile
  dens=findgen(n)*(densb-denst)/(n-1)+denst

  ;write everything to the file
  openw, lun, fname, /get_lun
  ;write out angles:
  printf, lun, 7
  printf, lun, 88.065, 0.17132
  printf, lun, 80.248, 0.36076
  printf, lun, 67.624, 0.46790
  printf, lun, 51.735, 0.46790
  printf, lun, 33.840, 0.36076
  printf, lun, 14.933, 0.17132
  printf, lun, 0., 0.

  ;number of levels:
  printf, lun, n+2
  ;sea-water first:
  printf, lun, -x[n-1]-0.2-dd, tsea, 1., 0., 1., 0., 0., 34., 0., 0., 1

  ;print out each of the levels, one by one:
  for i=n-1, 0, -1 do begin
    printf, lun, -x[i]-dd, temp[i], dens[i], 0., 0.1, dgr[i], dbub[i], s[i], 0., aspect_ratio, n-i+1
  endfor
  ;print out the last layer twice:
  printf, lun, -x[0], temp[0], dens[0], 0., 0.1, dgr[0], dbub[0], s[0], 0., aspect_ratio, n+2

  free_lun, lun

  ;stop

end


