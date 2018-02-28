pro plot_climatology_MAR2013,param,typee

; param = q, RH, Td, Tw, T, e, DPD
; typee = 'Annual','Seasonal'
nowmon='SEP'
nowyear='2014'

;climst=1976
;climed=2005
climst=1981
climed=2010

Domain = 'MARINE' ; 'LAND', 'MARINE'

IF (typee EQ 'Seasonal') THEN bits='climsSEAS' ELSE bits='climsANN'

CASE param OF
  'q': BEGIN
;    filein='/data/local/hadkw/HADCRUH2/UPDATE2013/STATISTICS/HadISDH.landq.2.0.0.2013p_FLATgridIDPHA5by5_JAN2014.nc'
;    fileout='/data/local/hadkw/HADCRUH2/UPDATE2013/IMAGES/HadISDH.landq.2.0.0.2013p_FLATgridIDPHA5by5_'+nowmon+nowyear+'_'+bits+'7605.eps'
    filein='/project/hadobs2/hadisdh/marine/ICOADS.2.5.1/GRIDS3/ERAclimNBC_5x5_monthly_from_daily_both_relax.nc'
    fileout='/data/local/hadkw/HADCRUH2/MARINE/IMAGES/ERAclimNBC_5x5_monthly_from_daily_both_relax_q_'+nowmon+nowyear+'_'+bits+'8110.eps'
    titlees='Surface Specific Humidity'
    units='(g kg!E-1!N)'    
    unitees='g/kg'
;    kcolsarrST=[255,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
;    levsarrST=[-2e+30,0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,14.0,16.0,18.0,20.0,21.0,30.0]
    kcolsarrST=[255,1,3,5,7,9,10,12,14,16,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,30.0]
    variname='specific humidity'
    letteree='a)'
  END
  'RH': BEGIN
;    filein='/data/local/hadkw/HADCRUH2/UPDATE2013/STATISTICS/HadISDH.landRH.2.0.0.2013p_FLATgridIDPHA5by5_JAN2014.nc'
;    fileout='/data/local/hadkw/HADCRUH2/UPDATE2013/IMAGES/HadISDH.landRH.2.0.0.2013p_FLATgridIDPHA5by5_'+nowmon+nowyear+'_'+bits+'7605.eps'
    filein='/project/hadobs2/hadisdh/marine/ICOADS.2.5.1/GRIDS3/ERAclimNBC_5x5_monthly_from_daily_both_relax.nc'
    fileout='/data/local/hadkw/HADCRUH2/MARINE/IMAGES/ERAclimNBC_5x5_monthly_from_daily_both_relax_RH_'+nowmon+nowyear+'_'+bits+'8110.eps'
    titlees='Surface Relative Humidity'
    units='(%rh)'   
    unitees='%rh' 
;    kcolsarrST=[255,1,2,3,4,5,6,7,9,10,11,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
;    levsarrST=[-2e+30,-100.0,5.0,10.,15.,20.0,25.,30.,40.0,50.,60.0,70.,75.,80.,85.,90.,95.0,110]
    kcolsarrST=[255,1,3,5,7,9,10,12,14,16,18]	; 14 colours + black to ensure no botching
;    levsarrST=[-2e+30,-100.0,10.,20.0,30.,40.0,50.,60.0,70.,80.,90.,110]
    levsarrST=[-2e+30,-100.0,55.,60.,65.0,70.,75.0,80.,85.0,90.,95.,110]
    variname='relative humidity'
    letteree='b)'
  END
  'e': BEGIN
;    filein='/data/local/hadkw/HADCRUH2/UPDATE2013/STATISTICS/HadISDH.lande.2.0.0.2013p_FLATgridIDPHA5by5_JAN2014.nc'
;    fileout='/data/local/hadkw/HADCRUH2/UPDATE2013/IMAGES/HadISDH.lande.2.0.0.2013p_FLATgridIDPHA5by5_'+nowmon+nowyear+'_'+bits+'7605.eps'
    filein='/project/hadobs2/hadisdh/marine/ICOADS.2.5.1/GRIDS3/ERAclimNBC_5x5_monthly_from_daily_both_relax.nc'
    fileout='/data/local/hadkw/HADCRUH2/MARINE/IMAGES/ERAclimNBC_5x5_monthly_from_daily_both_relax_e_'+nowmon+nowyear+'_'+bits+'8110.eps'
    titlees='Surface Vapour Pressure'
    units='(hPa)'  
    unitees='hPa'  
    kcolsarrST=[255,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,100.0]
    variname='vapour pressure'
    letteree='b)'
  END
  'Td': BEGIN
;    filein='/data/local/hadkw/HADCRUH2/UPDATE2013/STATISTICS/HadISDH.landTd.2.0.0.2013p_FLATgridPHADPD5by5_JAN2014.nc'
;    fileout='/data/local/hadkw/HADCRUH2/UPDATE2013/IMAGES/HadISDH.landTd.2.0.0.2013p_FLATgridPHADPD5by5_'+nowmon+nowyear+'_'+bits+'7605.eps'
    filein='/project/hadobs2/hadisdh/marine/ICOADS.2.5.1/GRIDS3/ERAclimNBC_5x5_monthly_from_daily_both_relax.nc'
    fileout='/data/local/hadkw/HADCRUH2/MARINE/IMAGES/ERAclimNBC_5x5_monthly_from_daily_both_relax_Td_'+nowmon+nowyear+'_'+bits+'8110.eps'
    titlees='Surface Dew Point Temperature'
    units='(!Eo!NC)'
    unitees='deg C'    
;    kcolsarrST=[255,1,2,3,4,5,6,7,9,10,11,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
;    levsarrST=[-2e+30,-40.0,-35.,-30.,-25.0,-20.,-15.0,-10.,-5., 0.0, 5.,10.,15.0,20.,25.0,30.,35.,100.]
    kcolsarrST=[255,1,3,5,7,9,10,12,14,16,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-40.0,-10.,-5., 0.0, 5.,10.,15.0,20.,25.0,30.,100.]
    variname='dew point temperature'
    letteree='c)' 
  END
  'Tw': BEGIN
;    filein='/data/local/hadkw/HADCRUH2/UPDATE2013/STATISTICS/HadISDH.landTw.2.0.0.2013p_FLATgridIDPHA5by5_JAN2014.nc'
;    fileout='/data/local/hadkw/HADCRUH2/UPDATE2013/IMAGES/HadISDH.landTw.2.0.0.2013p_FLATgridIDPHA5by5_'+nowmon+nowyear+'_'+bits+'7605.eps'
    filein='/project/hadobs2/hadisdh/marine/ICOADS.2.5.1/GRIDS3/ERAclimNBC_5x5_monthly_from_daily_both_relax.nc'
    fileout='/data/local/hadkw/HADCRUH2/MARINE/IMAGES/ERAclimNBC_5x5_monthly_from_daily_both_relax_Tw_'+nowmon+nowyear+'_'+bits+'8110.eps'
    titlees='Surface Wet Bulb Temperature'
    units='(!Eo!NC)'
    unitees='deg C'    
    kcolsarrST=[255,1,2,3,4,5,6,7,9,10,11,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-40.0,-35.,-30.,-25.0,-20.,-15.0,-10.,-5., 0.0,5.,10.,15.0,20.,25.0,30.,35.,100.]
    variname='wet bulb temperature'
    letteree='a)'
  END
  'T': BEGIN
;    filein='/data/local/hadkw/HADCRUH2/UPDATE2013/STATISTICS/HadISDH.landT.2.0.0.2013p_FLATgridIDPHA5by5_JAN2014.nc'
;    fileout='/data/local/hadkw/HADCRUH2/UPDATE2013/IMAGES/HadISDH.landT.2.0.0.2013p_FLATgridIDPHA5by5_'+nowmon+nowyear+'_'+bits+'7605.eps'
    filein='/project/hadobs2/hadisdh/marine/ICOADS.2.5.1/GRIDS3/ERAclimNBC_5x5_monthly_from_daily_both_relax.nc'
    fileout='/data/local/hadkw/HADCRUH2/MARINE/IMAGES/ERAclimNBC_5x5_monthly_from_daily_both_relax_T_'+nowmon+nowyear+'_'+bits+'8110.eps'
    titlees='Surface Temperature'
    units='(!Eo!NC)'    
    unitees='deg C'    
;    kcolsarrST=[255,1,2,3,4,5,6,7,9,10,11,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
;    levsarrST=[-2e+30,-40.0,-35.,-30.,-25.0,-20.,-15.0,-10.,-5.,0.0, 5.,10.,15.0,20.,25.0,30.,35.,100.]
    kcolsarrST=[255,1,3,5,7,9,10,12,14,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-40.0,-10.,-5.,0.0, 5.,10.0,15.0,20.,25.,30.,100.]
    variname='temperature'
    letteree='d)'
  END
  'DPD': BEGIN
;    filein='/data/local/hadkw/HADCRUH2/UPDATE2013/STATISTICS/HadISDH.landDPD.2.0.0.2013p_FLATgridPHA5by5_JAN2014.nc'
;    fileout='/data/local/hadkw/HADCRUH2/UPDATE2013/IMAGES/HadISDH.landDPD.2.0.0.2013p_FLATgridPHA5by5_'+nowmon+nowyear+'_'+bits+'7605.eps'
    filein='/project/hadobs2/hadisdh/marine/ICOADS.2.5.1/GRIDS3/ERAclimNBC_5x5_monthly_from_daily_both_relax.nc'
    fileout='/data/local/hadkw/HADCRUH2/MARINE/IMAGES/ERAclimNBC_5x5_monthly_from_daily_both_relax_DPD_'+nowmon+nowyear+'_'+bits+'8110.eps'
    titlees='Surface Dew Point Depression'
    units='(!Eo!NC)'    
    unitees='deg C'    
    kcolsarrST=[255,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,0,1.,2.0,3.,4.0,5.,6.0,7.,8.0,9.,10.0,11.,12.,13.,14.0,15.,16.0,17.,100.0]
    variname='dew point depression'
    letteree='c)'
  END

ENDCASE

mdi=-1e30
edyr=2015
styr=1973
nyrs=(edyr+1)-styr
nmons=nyrs*12
nlons=72
nlats=36
climtmp=make_array(nlons,nlats,12,/float,value=mdi)
maparrDJF=make_array(nlons,nlats,/float,value=mdi)
maparrMAM=make_array(nlons,nlats,/float,value=mdi)
maparrJJA=make_array(nlons,nlats,/float,value=mdi)
maparrSON=make_array(nlons,nlats,/float,value=mdi)
maparrANN=make_array(nlons,nlats,/float,value=mdi)

latlg=5.
lonlg=5.
stlt=-90+(latlg/2.)
stln=-180+(lonlg/2.)
nlats=180/latlg
nlons=360/lonlg
nbox=LONG(nlats*nlons)

lats=(findgen(nlats)*latlg)+stlt
lons=(findgen(nlons)*lonlg)+stln

djf=[11,0,1]
mam=[2,3,4]
jja=[5,6,7]
son=[8,9,10]

inn=NCDF_OPEN(filein)
IF (Domain EQ 'LAND') THEN BEGIN
  IF (climst EQ '1976') THEN BEGIN
    CASE param OF
      'q': varid=NCDF_VARID(inn,'q_clims')
    ;  'q': varid=NCDF_VARID(inn,'blendmask_q_clims')
    ;  'q': varid=NCDF_VARID(inn,'blend_q_clims')
      'RH': varid=NCDF_VARID(inn,'rh_clims')
      'e': varid=NCDF_VARID(inn,'e_clims')
      'Td': varid=NCDF_VARID(inn,'td_clims')
      'Tw': varid=NCDF_VARID(inn,'tw_clims')
      'T': varid=NCDF_VARID(inn,'t_clims')
      'DPD': varid=NCDF_VARID(inn,'dpd_clims')
    ENDCASE
  ENDIF ELSE BEGIN
    CASE param OF
      'q': varid=NCDF_VARID(inn,'q_anoms')
    ;  'q': varid=NCDF_VARID(inn,'blendmask_q_anoms')
    ;  'q': varid=NCDF_VARID(inn,'blend_q_anoms')
      'RH': varid=NCDF_VARID(inn,'rh_anoms')
      'e': varid=NCDF_VARID(inn,'e_anoms')
      'Td': varid=NCDF_VARID(inn,'td_anoms')
      'Tw': varid=NCDF_VARID(inn,'tw_anoms')
      'T': varid=NCDF_VARID(inn,'t_anoms')
      'DPD': varid=NCDF_VARID(inn,'dpd_anoms')
    ENDCASE
  ENDELSE
  latid=NCDF_VARID(inn,'lat')
  lonid=NCDF_VARID(inn,'lon')
  NCDF_VARGET,inn,varid,tmp
  NCDF_VARGET,inn,latid,lats
  NCDF_VARGET,inn,lonid,lons
  NCDF_CLOSE,inn
ENDIF ELSE BEGIN
    CASE param OF
      'q': varid=NCDF_VARID(inn,'specific_humidity')
      'RH': varid=NCDF_VARID(inn,'relative_humidity')
      'e': varid=NCDF_VARID(inn,'vapor_pressure')
      'Td': varid=NCDF_VARID(inn,'dew_point_temperature')
      'Tw': varid=NCDF_VARID(inn,'wet_bulb_temperature')
      'T': varid=NCDF_VARID(inn,'marine_air_temperature')
      'DPD': varid=NCDF_VARID(inn,'dew_point_depression')
    ENDCASE
  latid=NCDF_VARID(inn,'latitude')
  lonid=NCDF_VARID(inn,'longitude')
  NCDF_VARGET,inn,varid,tmp
  NCDF_VARGET,inn,latid,lats
  NCDF_VARGET,inn,lonid,lons
  NCDF_CLOSE,inn
  ;tmp = reverse(tmp,2)
  bads = where(tmp LT -100,count)
  IF (count GT 0) THEN tmp(bads) = mdi
ENDELSE

IF (Domain EQ 'LAND') AND (climst NE 1976) THEN BEGIN	; rezero over new climatology period
  FOR i=0,nlons-1 DO BEGIN
    FOR j=0,nlats-1 DO BEGIN
      newtmp=reform(tmp(i,j,*),12,nyrs)
      FOR mm=0,11 DO BEGIN
        climperiod=newtmp(mm,climst-styr:climed-styr)
	gotsclim=WHERE(climperiod NE mdi,count)
	IF (count GE 15) THEN BEGIN
	  gotsobs=where(newtmp(mm,*) NE mdi,count)
	  newtmp(mm,gotsobs)=newtmp(mm,gotsobs)-MEAN(climperiod(gotsclim))
	ENDIF
      ENDFOR
      tmp(i,j,*)=REFORM(newtmp,nmons)
    ENDFOR
  ENDFOR
ENDIF

; make clims for each month
IF (Domain EQ 'MARINE') THEN BEGIN	; rezero over new climatology period
  FOR i=0,nlons-1 DO BEGIN
    FOR j=0,nlats-1 DO BEGIN
      newtmp=reform(tmp(i,j,*),12,nyrs)
      FOR mm=0,11 DO BEGIN
        climperiod=newtmp(mm,climst-styr:climed-styr)
	gotsclim=WHERE(climperiod NE mdi,count)
	IF (count GE 10) THEN BEGIN
	  gotsobs=where(newtmp(mm,*) NE mdi,count)
	  climtmp(i,j,mm)=MEAN(climperiod(gotsclim))
	ENDIF
      ENDFOR
    ENDFOR
  ENDFOR
ENDIF

; average over the year for each gridbox

FOR ltt=0,nlats-1 DO BEGIN
  FOR lnn=0,nlons-1 DO BEGIN
    IF (typee EQ 'Seasonal') THEN BEGIN
      subarr=climtmp(lnn,ltt,mam)
      gots=WHERE(subarr NE mdi,count)
      IF (count GE 2) THEN maparrMAM(lnn,ltt)=MEAN(subarr(gots)) ;2
      subarr=climtmp(lnn,ltt,jja)
      gots=WHERE(subarr NE mdi,count)
      IF (count GE 2) THEN maparrJJA(lnn,ltt)=MEAN(subarr(gots)) ;2
      subarr=climtmp(lnn,ltt,son)
      gots=WHERE(subarr NE mdi,count)
      IF (count GE 2) THEN maparrSON(lnn,ltt)=MEAN(subarr(gots)) ;2
      subarr=climtmp(lnn,ltt,djf)
      gots=WHERE(subarr NE mdi,count)
      IF (count GE 2) THEN maparrDJF(lnn,ltt)=MEAN(subarr(gots)) ;2
    ENDIF ELSE BEGIN
      subarr=climtmp(lnn,ltt,*)
      gots=WHERE(subarr NE mdi,count)
      IF (count GE 9) THEN maparrANN(lnn,ltt)=MEAN(subarr(gots)) ;9
    ENDELSE
  ENDFOR
ENDFOR

; colour settings - Uni ORegon blue to red
tvlct,0,0,0,0
tvlct,240,240,240,100

ncols=n_elements(kcolsarrST)
colsarrST=kcolsarrST(1:ncols-1)
nlevs=n_elements(levsarrST)-1
labsarrST=strarr(nlevs)
labsarrST(0)=''
labsarrST((nlevs-1))=''
labsarrST(1:nlevs-2)=string(levsarrST(2:nlevs-1),format='(i4)')

;------------------------------------------------------------------------------------
IF (typee EQ 'Seasonal') THEN BEGIN
  xpos1=[0.05,0.55,0.05,0.55]
  xpos2=[0.45,0.95,0.45,0.94]
  ypos1=[0.60,0.6,0.15,0.15]
  ypos2=[0.95,0.95,0.5,0.5]
ENDIF ELSE BEGIN
  xpos1=0.07
  xpos2=0.93
  ypos1=0.12
  ypos2=0.92
ENDELSE
;make the plot
set_plot,'PS'
device,filename=fileout,/ENCAPSUL,$
       xsize=26,ysize=20,/landscape,/color,/helvetica,/bold


!P.Font=0
!P.Thick=4
!X.Thick=4
!Y.Thick=4

tvlct,200,200,200,100
IF (param EQ 'q') OR (param EQ 'RH') OR (param EQ 'e') OR (param EQ 'Td') OR (param EQ 'Tw') THEN BEGIN
  tvlct,10,0,0,1
  tvlct,51,25,0,2
  tvlct,75,32,0,3
  tvlct,102,47,0,4
  tvlct,128,73,26,5
  tvlct,153,96,53,6
  tvlct,204,155,122,7
  tvlct,216,175,151,8
  tvlct,242,218,205,9
  tvlct,204,253,255,10
  tvlct,153,248,255,11
  tvlct,101,239,255,12
  tvlct,50,227,255,13
  tvlct,0,169,204,14
  tvlct,0,122,153,15
  tvlct,0,75,100,16
  tvlct,0,40,80,17
  tvlct,0,10,40,18
ENDIF ELSE IF (param EQ 'DPD') THEN BEGIN
;  tvlct,10,0,0,18
;  tvlct,20,1,0,17
;  tvlct,30,5,0,16
;  tvlct,51,25,0,15
;  tvlct,102,47,0,14
  tvlct,10,0,0,18
  tvlct,51,25,0,17
  tvlct,75,32,0,16
  tvlct,102,47,0,15
  tvlct,128,73,26,14
  tvlct,153,96,53,13
  tvlct,204,155,122,12
  tvlct,216,175,151,11
  tvlct,242,218,205,10
  tvlct,204,253,255,9
  tvlct,153,248,255,8
  tvlct,101,239,255,7
  tvlct,50,227,255,6
  tvlct,0,169,204,5
  tvlct,0,122,153,4
  tvlct,0,75,100,3
  tvlct,0,40,80,2
  tvlct,0,10,40,1
ENDIF ELSE IF (param EQ 'T') THEN BEGIN
    tvlct,0,0,110,1
    tvlct,36,0,216,2
    tvlct,24,28,247,3
    tvlct,40,87,255,4
    tvlct,61,135,255,5
    tvlct,86,176,255,6
    tvlct,117,211,255,7
    tvlct,153,234,255,8
    tvlct,188,249,255,9
    tvlct,255,241,188,10
    tvlct,255,214,153,11
    tvlct,255,172,117,12
    tvlct,255,120,86,13
    tvlct,255,61,61,14
    tvlct,247,39,53,15
    tvlct,216,21,47,16
    tvlct,165,0,33,17
    tvlct,70,0,0,18
ENDIF

IF (typee EQ 'Seasonal') THEN titlees=['DJF','MAM','JJA','SON'] ELSE titlees=variname


IF (typee EQ 'Seasonal') THEN BEGIN
  !P.Position=[xpos1(0),ypos1(0),xpos2(0),ypos2(0)]

  Map_set,/continents,color=0,/noborder,/noerase,/robinson,/grid, glinestyle=2
  XYOUTS,xpos1(0)+((xpos2(0)-xpos1(0))/2.),ypos2(0)+0.02,titlees(0),/NORMAL,alignment=0.5,charsize=1.2,color=0
  XYOUTS,xpos1(0)-0.01,ypos2(0)-0.02,'a)',/NORMAL,alignment=0.5,charsize=1.2,color=0

  boxfill,maparrDJF,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
  Map_set,/continents,color=0,/noborder,/noerase,/robinson,/grid, glinestyle=2

  PLOTS,[-179.9,0],[89.9,89.9],color=0
  PLOTS,[0,179.9],[89.9,89.9],color=0
  PLOTS,[-179.9,0],[-89.9,-89.9],color=0
  PLOTS,[0,179.9],[-89.9,-89.9],color=0

  !P.Position=[xpos1(1),ypos1(1),xpos2(1),ypos2(1)]

  Map_set,/continents,color=0,/noborder,/noerase,/robinson,/grid, glinestyle=2
  XYOUTS,xpos1(1)+((xpos2(1)-xpos1(1))/2.),ypos2(1)+0.01,titlees(1),/NORMAL,alignment=0.5,charsize=1.2,color=0
  XYOUTS,xpos1(1)-0.01,ypos2(1)-0.02,'b)',/NORMAL,alignment=0.5,charsize=1.2,color=0

  boxfill,maparrMAM,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
  Map_set,/continents,color=0,/noborder,/noerase,/robinson,/grid, glinestyle=2

  PLOTS,[-179.9,0],[89.9,89.9],color=0
  PLOTS,[0,179.9],[89.9,89.9],color=0
  PLOTS,[-179.9,0],[-89.9,-89.9],color=0
  PLOTS,[0,179.9],[-89.9,-89.9],color=0

  !P.Position=[xpos1(2),ypos1(2),xpos2(2),ypos2(2)]

  Map_set,/continents,color=0,/noborder,/noerase,/robinson,/grid, glinestyle=2
  XYOUTS,xpos1(2)+((xpos2(2)-xpos1(2))/2.),ypos2(2)+0.01,titlees(2),/NORMAL,alignment=0.5,charsize=1.2,color=0
  XYOUTS,xpos1(2)-0.01,ypos2(2)-0.02,'c)',/NORMAL,alignment=0.5,charsize=1.2,color=0

  boxfill,maparrJJA,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
       
  ; Overplot continents and islands again
  Map_set,/continents,color=0,/noborder,/noerase,/robinson,/grid, glinestyle=2

  PLOTS,[-179.9,0],[89.9,89.9],color=0
  PLOTS,[0,179.9],[89.9,89.9],color=0
  PLOTS,[-179.9,0],[-89.9,-89.9],color=0
  PLOTS,[0,179.9],[-89.9,-89.9],color=0

  !P.Position=[xpos1(3),ypos1(3),xpos2(3),ypos2(3)]

  Map_set,/continents,color=0,/noborder,/noerase,/robinson,/grid, glinestyle=2
  XYOUTS,xpos1(3)+((xpos2(3)-xpos1(3))/2.),ypos2(3)+0.01,titlees(3),/NORMAL,alignment=0.5,charsize=1.2,color=0
  XYOUTS,xpos1(3)-0.01,ypos2(3)-0.02,'d)',/NORMAL,alignment=0.5,charsize=1.2,color=0

  boxfill,maparrSON,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
  Map_set,/continents,color=0,/noborder,/noerase,/robinson,/grid, glinestyle=2

  PLOTS,[-179.9,0],[89.9,89.9],color=0
  PLOTS,[0,179.9],[89.9,89.9],color=0
  PLOTS,[-179.9,0],[-89.9,-89.9],color=0
  PLOTS,[0,179.9],[-89.9,-89.9],color=0

ENDIF ELSE BEGIN
  !P.Position=[xpos1(0),ypos1(0),xpos2(0),ypos2(0)]

  XYOUTS,xpos1(0)+((xpos2(0)-xpos1(0))/2.),ypos2(0)+0.02,titlees(0),/NORMAL,alignment=0.5,charsize=2.5,color=0
  Map_set,/continents,color=0,/noborder,/noerase,/robinson,/grid, glinestyle=2
  boxfill,maparrANN,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
  Map_set,/continents,color=0,/noborder,/noerase,/robinson,/grid, glinestyle=2

  PLOTS,[-179.9,0],[89.9,89.9],color=0
  PLOTS,[0,179.9],[89.9,89.9],color=0
  PLOTS,[-179.9,0],[-89.9,-89.9],color=0
  PLOTS,[0,179.9],[-89.9,-89.9],color=0
ENDELSE

MAKE_KEY,0.05,0.08,0.9,0.03,-0.01,-0.035,/NORMAL,COLORS=colsarrST,labels=labsarrST,$
         charsize=1.8,charthick=4,bcolor=0  ;,orientation=1

IF (typee EQ 'Seasonal') THEN BEGIN
  XYOUTS,0.5,0.015,'Seasonal Climatology over '+strcompress(climst,/remove_all)+'-'$
                +strcompress(climed,/remove_all)+' '+units,/normal,color=0,$
		charsize=1.2,alignment=0.5 
ENDIF ELSE BEGIN
  XYOUTS,0.5,0.015,'annual climatology over '+strcompress(climst,/remove_all)+'-'$
                +strcompress(climed,/remove_all)+' '+units,/normal,color=0,$
		charsize=1.8,alignment=0.5
  XYOUTS,0.05,0.9,letteree,/normal,color=0,charsize=2.5,alignment=0.5  
ENDELSE

DEVICE,/close


return

end

