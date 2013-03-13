
n=0L
openr, 1, "icetopo.txt"
readf, 1, n
data=fltarr(3, n)
readf, 1, data
close, 1
s3=reform(data[0, *])
h1=reform(data[1, *])
z1=reform(data[2, *])

n=0L
openr, 1, "pissfit.txt"
readf, 1, n
data=fltarr(2, n)
readf, 1, data
s1=reform(data[0, *])
h=reform(data[1, *])
readf, 1, n
data=fltarr(2, n)
readf, 1, data
s2=reform(data[0, *])
z=reform(data[1, *])
close, 1

xrange=[13050, 13100]

;ps_start, "piss_traj.ps"

round_sym, 0.5
plot, s1, h, yrange=[-3, 3], xrange=xrange, xstyle=1, ystyle=1
oplot, s2, z
oplot, s3, h1, psym=8
oplot, s3, z1, psym=8

openr, 1, "piss_traj.txt"

line=""

x=0.
y=0.

newflag=1

ntraj=0

while eof(1) eq 0 do begin
  readf, 1, line
  if line ne "" then begin
    reads, line, x, y
    if newflag eq 1 then begin
      ;plot, s1, h, yrange=[-3, 3], xrange=xrange
      ;oplot, s2, z
      newflag=0
      ntraj=ntraj+1
    endif else begin
      oplot, [lastx, x], [lasty, y]
    endelse
    lastx=x
    lasty=y
  endif else begin
    newflag=1
    ;stop
  endelse

  ;if ntraj mod 100 eq 0 then stop

endwhile

close, 1

;ps_end

end
