@create_salinity_prof

;calculates a representative salinity profile based on
;Eicken (1992) 'S' type, but designed so that:
;1. there is a constant offset that varies with salinity
;2. below 2 psu the profile is flat

function calc_sprof, s, n
  if n_elements(n) ne 1 then n=10

  ;verticle coordinate is arbitrary:
  x=findgen(n)/(n-1)

  s0prof=create_salinity_prof(x, 6., 'S')

  if s lt 2 then begin
    sprof=fltarr(n)+s
  endif else begin

    sprof1=(s0prof-2)*(s-2)+s/3.

    sprof=n*sprof1*(s-2)/total(sprof1)+2

  endelse

  return, sprof

end

;little test routine:
s=findgen(21)/2
ps_start, "sprof.ps", /portrait
plot, findgen(10), xrange=[0, 20], yrange=[-1, 0], /nodata, $
		xtitle="salinity", ytitle="z", charsize=1.5
for i=0, 20 do begin
  sprof=calc_sprof(s[i])
  print, total(sprof)/10

  oplot, sprof, -findgen(10)/9
endfor

ps_end

end

