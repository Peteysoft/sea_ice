;reads the averaged data from Lars:
pro read_averaged, tbv, tbh, d, dstd, filename=filename, $
		tbh_std=tbh_std, tbv_std=tbv_std, irange=irange, $
		track=track

  if n_elements(filename) ne 1 then filename="/freax/storage/home/pmills/smos/data/police_campaign_averaged.txt"

  openr, lun, filename, /get_lun

  data=fltarr(13, 32)
  readf, lun, data

  tbv=reform(data[4, *])
  tbh=reform(data[6, *])

  tbh_std=reform(data[5, *])
  tbv_std=reform(data[7, *])

  d=reform(data[8, *])
  dstd=reform(data[9, *])

  track=reform(data[10, *])
  irange=data[11:12, *]

  free_lun, lun

end

