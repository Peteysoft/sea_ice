;calculates the emissivity of the 'U' component 
;from a sloping surface

;cross product of a vector:
function cross, v1, v2
  return, [v1[1]*v2[2]-v1[2]*v2[1], v1[0]*v2[2]-v1[2]*v2[0], $
		v1[0]*v2[1]-v1[1]*v2[2]]

end

function dot, v1, v2
  return, total(v1*v2)
end

