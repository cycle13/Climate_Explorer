pro nino_indicators_Kate2015

outdir = '/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/'
climato = '6190'
regione=[120,-20,-70,20]

hadisst_filename = '/project/hadobs1/OBS/marine/HadISST/anoms/HadISST1.1_sst_1870on_1dg_anm6190.pp'

if climato eq '6190' then begin
   allclim = ppa(hadisst_filename,'f.lbyr ge 1961 and f.lbyr le 1990',region=regione,/all)
   climatology = allclim(0:11)
   for m = 1,12 do begin
      climatology(m-1) = pp_avg(allclim(where(allclim.lbmon eq m)),/mdtol)
   endfor
   climindc = '1961-1990'
endif

if climato eq '8110' then begin
   allclim = ppa(hadisst_filename,'f.lbyr ge 1981 and f.lbyr le 2010',region=regione,/all)
   climatology = allclim(0:11)
   for m = 1,12 do begin
      climatology(m-1) = pp_avg(allclim(where(allclim.lbmon eq m)),/mdtol)
   endfor
   climindc = '1981-2010'
endif

if climato eq '7100' then begin
   allclim = ppa(hadisst_filename,'f.lbyr ge 1971 and f.lbyr le 2000',region=regione,/all)
   climatology = allclim(0:11)
   for m = 1,12 do begin
      climatology(m-1) = pp_avg(allclim(where(allclim.lbmon eq m)),/mdtol)
   endfor
   climindc = '1971-2000'
endif

if climato eq '7605' then begin
   allclim = ppa(hadisst_filename,'f.lbyr ge 1976 and f.lbyr le 2005',region=regione,/all)
   climatology = allclim(0:11)
   for m = 1,12 do begin
      climatology(m-1) = pp_avg(allclim(where(allclim.lbmon eq m)),/mdtol)
   endfor
   climindc = '1976-2005'
endif


anom = ppa('/project/hadobs1/OBS/marine/HadISST/anoms/HadISST1.1_sst_1870on_1dg_anm6190.pp','f.lbyr ge 1870',region=[120,-20,-70,20],/all)

next_month = anom(n_elements(anom)-1).lbmon+1
next_year = anom(n_elements(anom)-1).lbyr

if next_month eq 13 then begin
    next_month = 1
    next_year = next_year+1
endif

snext_month = nicenumber(next_month)
if next_month lt 10 then snext_month = '0'+snext_month
snext_year = nicenumber(next_year)

filename = '/project/hadobs1/OBS/marine/HadISST/temp/HadISST1_forAATSRval+insitu_'+snext_year+snext_month+'_anm6190.pp'

if canopen(filename) then begin
    latest = ppa(filename,region=[120,-20,-70,20],/all)
    all = [anom,latest]
endif else begin
    all = anom
endelse

;reanomalise
for i = 0,n_elements(all)-1 do begin
   all(i) = pp_ff('a-b',all(i),climatology(all(i).lbmon-1),/quiet)
endfor

nino3 = tsarea_avg(all,[-150,-5,-90,5])
nino4 = tsarea_avg(all,[160,-5,-150,5])
nino34 = tsarea_avg(all,[-170,-5,-120,5])
nino12 = tsarea_avg(all,[-90,-10,-80,0])

;stop

; Write out
openw,1,outdir+'HadISST1.1_Nino3.4ts_18502015_6190_NOV2015.txt'
openw,2,outdir+'HadISST1.1_Nino3ts_18502015_6190_NOV2015.txt'
openw,3,outdir+'HadISST1.1_Nino4ts_18502015_6190_NOV2015.txt'
openw,4,outdir+'HadISST1.1_Nino1.2ts_18502015_6190_NOV2015.txt'
yearsie=1870
monsie=0
monarr=['01','02','03','04','05','06','07','08','09','10','11','12']
for mm = 0,n_elements(nino3.data)-1 do begin
  printf,1,strcompress(yearsie,/remove_all),monarr(monsie),nino34.data(mm),format='(a4,a2,x,f7.3)'
  printf,2,strcompress(yearsie,/remove_all),monarr(monsie),nino3.data(mm),format='(a4,a2,x,f7.3)'
  printf,3,strcompress(yearsie,/remove_all),monarr(monsie),nino4.data(mm),format='(a4,a2,x,f7.3)'
  printf,4,strcompress(yearsie,/remove_all),monarr(monsie),nino12.data(mm),format='(a4,a2,x,f7.3)'
  monsie=monsie+1
  if (monsie EQ 12) then begin
    monsie=0
    yearsie=yearsie+1
  endif
endfor

close,1
close,2
close,3
close,4
end
