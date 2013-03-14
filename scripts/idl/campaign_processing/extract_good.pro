;extract the good data from the averaged stuff:

restore, "collave.idlsav"

ind=where(tbs_std[0, *] lt 10 and tbs_std[1, *] lt 10 and za_std lt 1 and $
		(fno_bin eq 7 and pindex1 gt 20000 and pindex2 lt 22000) eq 0 and $
		(fno_bin eq 6 and pindex1 gt 13000 and pindex2 lt 17500) eq 0 and $
		(fno_bin eq 1 and pindex2 lt 25000) eq 0 and $
		fno_bin ne 3 and fno_bin ne 2 and $
		(fno_bin eq 0 and pindex1 gt 20000) eq 0)

dbin=dbin[ind]
fno_bin=fno_bin[ind]
lat_bin=lat_bin[ind]
lon_bin=lon_bin[ind]
nr_bin=nr_bin[ind]
pindex1=pindex1[ind]
pindex2=pindex2[ind]
stdbin=stdbin[ind]
tbs_bin=tbs_bin[*, ind]
ts_bin=ts_bin[ind]
za_bin=za_bin[ind]
za_std=za_std[ind]
tbs_std=tbs_std[*, ind]
tbin=tbin[ind]

save, filename="collave_good.idlsav", dbin, fno_bin, lat_bin, lon_bin, nr_bin, $
		pindex1, pindex2, stdbin, tbs_bin, ts_bin, za_bin, za_std, tbs_std, $
		rname, tbin

end

