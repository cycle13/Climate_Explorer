pro run_annualanommaps


; STUFF FOR LAND
version='4.0.0.2017f'
nowmon='JAN'
nowyear='2018'
thenmon='JAN'
thenyear='2018'
homogtype='IDPHA' ; 'PHA','PHADPD','IDPHA'
param='RH'	;q,e,RH,T,Tw,Td,DPD
climst=1981			;1981, 1976
climed=2010			;2010, 2005
climie=strmid(strcompress(climst,/remove_all),2,2)+strmid(strcompress(climed,/remove_all),2,2)

; NOTE!!! Code can renormalise to new period but this is not necessary with HadISDH from 2017 onwards

;; STUFF FOR BLEND
;version='1.0.0.2016p'
;nowmon='JAN'
;nowyear='2017'
;thenmon='JAN'
;thenyear='2017'
;homogtype='FULL' ; 'FULL','FULLship'
;param='q'	;q,e,RH,T,Tw,Td,DPD
;climst=1981			;1981, 1976
;climed=2010			;2010, 2005
;climie=strmid(strcompress(climst,/remove_all),2,2)+strmid(strcompress(climed,/remove_all),2,2)

;; STUFF FOR MARINE
;version='1.0.0.2016p'
;nowmon='JAN'
;nowyear='2017'
;thenmon='JAN'
;thenyear='2017'
;;homogtype='OBSclim2BClocal' ; 'OBSclim2BClocal', 'OBSclim2BClocalship'
;param='q'	;q,e,RH,T,Tw,Td,DPD
;climst=1981			;1981, 1976
;climed=2010			;2010, 2005
;climie=strmid(strcompress(climst,/remove_all),2,2)+strmid(strcompress(climed,/remove_all),2,2)

for yy=1973,2017 DO BEGIN
  dir='/data/local/hadkw/HADCRUH2/UPDATE2017/'
  ystr=string(yy,format='(i4)')
  print,yy,ystr

  ; LAND
  filein=dir+'STATISTICS/GRIDS/HadISDH.land'+param+'.'+version+'_FLATgrid'+homogtype+'5by5_anoms8110_'+thenmon+thenyear+'_cf.nc'
  plotout=dir+'IMAGES/MAPS/HadISDH.land'+param+'.'+version+'_FLATgrid'+homogtype+'5by5_anoms8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
  fileout=dir+'STATISTICS/GRIDS/HadISDH.land'+param+'.'+version+'_FLATgrid'+homogtype+'5by5_anoms8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr

;  ; MARINE
;  filein=dir+'STATISTICS/GRIDS/HadISDH.marine'+param+'.'+version+'_'+homogtype+'_anoms8110_'+thenmon+thenyear+'_cf.nc'
;  plotout=dir+'IMAGES/MAPS/HadISDH.marine'+param+'.'+version+'_'+homogtype+'_anoms8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  fileout=dir+'STATISTICS/GRIDS/HadISDH.marine'+param+'.'+version+'_'+homogtype+'_anoms8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr

;  ; BLEND
;  filein=dir+'STATISTICS/GRIDS/HadISDH.blend'+param+'.'+version+'_'+homogtype+'_anoms8110_'+thenmon+thenyear+'_cf.nc'
;  plotout=dir+'IMAGES/MAPS/HadISDH.blend'+param+'.'+version+'_'+homogtype+'_anoms8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  fileout=dir+'STATISTICS/GRIDS/HadISDH.blend'+param+'.'+version+'_'+homogtype+'_anoms8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr

;  filein=dir+'STATISTICS/GRIDS/BLEND_NOCSv2.0_HadISDH.land'+param+'.'+version+'_'+nowmon+nowyear+'.nc'
;  plotout=dir+'IMAGES/MAPS/BLEND_NOCSv2.0_HadISDH.land'+param+'.'+version+'_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  fileout=dir+'STATISTICS/GRIDS/BLEND_NOCSv2.0_HadISDH.land'+param+'.'+version+'_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  plotout=dir+'IMAGES/MAPS/BLENDMASK_NOCSv2.0_HadISDH.land'+param+'.'+version+'_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  fileout=dir+'STATISTICS/GRIDS/BLENDMASK_NOCSv2.0_HadISDH.land'+param+'.'+version+'_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr

;  filein=dir+'STATISTICS/NOCSv2.0_oceanW_8110_'+nowmon+nowyear+'.nc'
;  plotout=dir+'IMAGES/NOCSv2.0_oceanW_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  fileout=dir+'STATISTICS/NOCSv2.0_oceanW_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  plotout=dir+'IMAGES/NOCSv2.0_MASKSEDoceanW_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  fileout=dir+'STATISTICS/NOCSv2.0_MASKEDoceanW_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  plotout=dir+'IMAGES/NOCSv2.0_QoceanW_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  fileout=dir+'STATISTICS/NOCSv2.0_QoceanW_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr

;  filein=dir+'OTHERDATA/waswind_v1_0_1.monthly.nc'
;  plotout=dir+'IMAGES/waswind_oceanW_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  fileout=dir+'OTHERDATA/waswind_oceanW_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr

;  filein=dir+'OTHERDATA/ERAINTERIM_q_5by519812010clim7913anomskate_19732014.nc'
;  plotout=dir+'IMAGES/ERAINTERIM_q_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  fileout=dir+'OTHERDATA/ERAINTERIM_q_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr

;  filein=dir+'OTHERDATA/ERAINTERIM_RH_5by519812010clim7913anomskate_19732014.nc'
;  plotout=dir+'IMAGES/ERAINTERIM_RH_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  fileout=dir+'OTHERDATA/ERAINTERIM_RH_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr

;  filein=dir+'OTHERDATA/ERAINTERIM_T_5by519812010clim7913anomskate_19732014.nc'
;  plotout=dir+'IMAGES/ERAINTERIM_T_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  fileout=dir+'OTHERDATA/ERAINTERIM_T_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr

;  filein=dir+'OTHERDATA/HadCRUT.4.3.0.0.median.nc'
;  plotout=dir+'IMAGES/MAPS/ANOMS8110/HadCRUI4300_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr
;  fileout=dir+'OTHERDATA/STATISTICS/HadCRUT4300_8110_'+nowmon+nowyear+'_annualanom'+climie+'_'+ystr


  plot_annualanommaps_FEB2015,string(yy,format='(i4)'),filein,fileout,plotout,climst,climed,param		;string(yy,format='(a4)')
;  stop
  fileouteps=plotout+'.eps'
  fileoutpng=plotout+'.png'
  spawn,'convert -rotate -90 '+fileouteps+' '+fileoutpng+'
endfor

end
