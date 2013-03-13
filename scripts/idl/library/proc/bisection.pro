
function bisection, f, x1, x2, tol, params=params
  f1=call_function(f, x1, params=params)
  f2=call_function(f, x2, params=params)

  if (f1 gt 0 or f2 lt 0) then begin
    message, "Root must be bracketed"
  endif

  while abs(x1-x2) gt tol do begin
    ;print, f1, f2
    xmid=(x1+x2)/2.
    fmid=call_function(f, xmid, params=params)

    ;stop

    if fmid lt 0 then begin
      x1=xmid
      f1=fmid
    endif else if fmid ge 0 then begin
      x2=xmid
      f2=fmid
    endif
  endwhile

  return, (x1+x2)/2.

end


