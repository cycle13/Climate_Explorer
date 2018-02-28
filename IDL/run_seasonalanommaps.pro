pro run_seasonalanommaps

version='4.0.0.2017f'
nowmon='JAN'
nowyear='2018'
thenmon='JAN'
thenyear='2018'
homogtype='IDPHA'
param='T'
climst=1981
climed=2010
climie=strmid(strcompress(climst,/remove_all),2,2)+strmid(strcompress(climed,/remove_all),2,2)

for yy=1973,2017 DO BEGIN
  dir='/data/local/hadkw/HADCRUH2/UPDATE2017/'
  ystr=string(yy,format='(i4)')
  print,yy,ystr
  filein=dir+'STATISTICS/GRIDS/HadISDH.land'+param+'.'+version+'_FLATgrid'+homogtype+'5by5_anoms8110_'+thenmon+thenyear+'_cf.nc'
  plotout=dir+'IMAGES/MAPS/ANOMS8110/HadISDH.land'+param+'.'+version+'_FLATgrid'+homogtype+'5by5_anoms8110_'+nowmon+nowyear+'_seasonalanom'+climie+'_'+ystr
  fileout=dir+'STATISTICS/GRIDS/ANOMS8110/HadISDH.land'+param+'.'+version+'_FLATgrid'+homogtype+'5by5_anoms8110_'+nowmon+nowyear+'_seasonalanom'+climie+'_'+ystr

;  filein=dir+'STATISTICS/GRIDS/HadISDH.land'+param+'.'+version+'_FLATgrid'+homogtype+'5by5_'+thenmon+thenyear+'_cf.nc'
;  plotout=dir+'IMAGES/MAPS/ANOMS7605/HadISDH.land'+param+'.'+version+'_FLATgrid'+homogtype+'5by5_'+nowmon+nowyear+'_seasonalanom'+climie+'_'+ystr
;  fileout=dir+'STATISTICS/GRIDS/ANOMS7605/HadISDH.land'+param+'.'+version+'_FLATgrid'+homogtype+'5by5_'+nowmon+nowyear+'_seasonalanom'+climie+'_'+ystr

;  filein=dir+'STATISTICS/GRIDS/BLEND_NOCSv2.0_HadISDH.land'+param+'.'+version+'_'+nowmon+nowyear+'.nc'
;;  plotout=dir+'IMAGES/MAPS/ANOMS8110/BLEND_NOCSv2.0_HadISDH.land'+param+'.'+version+'_8110_'+nowmon+nowyear+'_seasonalanom'+climie+'_'+ystr
;;  fileout=dir+'STATISTICS/GRIDS/ANOMS8110/BLEND_NOCSv2.0_HadISDH.land'+param+'.'+version+'_8110_'+nowmon+nowyear+'_seasonalanom'+climie+'_'+ystr
;  plotout=dir+'IMAGES/MAPS/ANOMS8110/BLENDMASK_NOCSv2.0_HadISDH.land'+param+'.'+version+'_'+nowmon+nowyear+'_seasonalanom'+climie+'_'+ystr
;  fileout=dir+'STATISTICS/GRIDS/ANOMS8110/BLENDMASK_NOCSv2.0_HadISDH.land'+param+'.'+version+'_'+nowmon+nowyear+'_seasonalanom'+climie+'_'+ystr

  plot_seasonalanommaps_FEB2014,string(yy,format='(i4)'),filein,fileout,plotout,climst,climed,param		;string(yy,format='(a4)')
;  stop
  fileouteps=plotout+'.eps'
  fileoutpng=plotout+'.png'
  spawn,'convert -rotate -90 '+fileouteps+' '+fileoutpng+'
endfor

end
