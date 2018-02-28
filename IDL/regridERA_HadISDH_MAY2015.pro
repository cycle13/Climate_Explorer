pro regridERA_HadISDH_MAY2015

; regrid ERA monthly 1by1 data to monthly 5by5 data
; make anomalies for 1981-2010 and 1979-2005

;--------------------------------------------------------
indir='/data/local/hadkw/HADCRUH2/UPDATE2017/OTHERDATA/'
infil='_monthly_1by1_ERA-Interim_data_19792017.nc'
;infil='ERAINTERIM_evap_monthly00_12hr_1by1_197901201512.nc'
inlandmask='HadCRUT.4.3.0.0.land_fraction.nc'
outfil='_monthly_5by5_ERA-Interim_data_19792017.nc'
outfilA='_monthly_5by5_ERA-Interim_data_19792017_anoms1981-2010.nc'
outfilA7905='_monthly_5by5_ERA-Interim_data_19792017_anoms1979-2005.nc'

varlist=['q2m','rh2m','e2m','t2m','td2m','tw2m','dpd2m','evap']
unitslist=['g/kg','%rh','hPa','deg C','deg C','deg C','deg C','cm w.e']
mdi=-1e+30

monarr=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
styr=1979
edyr=2017
climst1=1981
climed1=2010
climst2=1979
climed2=2005
nyrs=(edyr+1)-styr
nmons=nyrs*12
int_mons=indgen(nmons)
dayssince=findgen(nmons)
stpoint=JULDAY(1,1,1979)

; array of months in days since Jan 1st 1979 including Jan 2015 (nmons+1)
monthies=TIMEGEN(nmons+1,START=(stpoint),UNITS="Months")-stpoint ; an array of floats describing each month in days since Jan 1st 1973
monthies=ROUND(monthies*100.)/100.	; nicer numbers

; Now need to use the mid-point of each month.
FOR mm=0,nmons-1 DO BEGIN
  dayssince(mm)=monthies(mm)+(monthies(mm+1)-monthies(mm))/2.
ENDFOR

latlg=5.
lonlg=5.
stlt=-90+(latlg/2.)
stln=-180+(lonlg/2.)
nlats=180/latlg
nlons=360/lonlg
nbox=LONG(nlats*nlons)

lats=(findgen(nlats)*latlg)+stlt
lons=(findgen(nlons)*lonlg)+stln

absarr=fltarr(nlons,nlats,nmons)
;-----------------------------------------------------------
FOR loo=0,6 DO BEGIN ; loop through the variables 7,7 for just evap
  ; Read in 1by1 file
  IF (varlist[loo] NE 'evap') THEN BEGIN
    inn=NCDF_OPEN(indir+varlist[loo]+infil) 
    varid=NCDF_VARID(inn,varlist[loo])
    NCDF_VARGET,inn,varid,hiabs
  ENDIF ELSE BEGIN
    inn=NCDF_OPEN(indir+infil)
    varid=NCDF_VARID(inn,'e')
    NCDF_VARGET,inn,varid,oldhiabs
    NCDF_ATTGET,inn,varid,'scale_factor',sf
    NCDF_ATTGET,inn,varid,'add_offset',off
    oldhiabs=(oldhiabs*sf)+off
    ; This is actually two vals per month at 00hrs and 12hrs so need to be combined
    hiabs=fltarr(360,181,nmons)
    FOR lnn=0,359 DO BEGIN
      FOR ltt=0,180 DO BEGIN
        hiabs(lnn,ltt,*)=REFORM(mean(REFORM(oldhiabs(lnn,ltt,*),2,nmons),dimension=1),nmons)
      ENDFOR
    ENDFOR  
  ENDELSE
  NCDF_CLOSE,inn
  
  ; Regrid to 5by5 by simple averaging
  ; ERAI is 181 lats 90.5 to -90.5 (GB centres 90 to -90) therefore we need to make it 0.5 degrees and then average over box
  ; ERAI is 360 longs -0.5 to 359.5 (GB centres 0 to 359) therefore we need to make it 0.5 degrees and average over box
  hlt=0
  hln=0
  FOR ltt=0,nlats-1 DO BEGIN
    FOR lnn=0,nlons-1 DO BEGIN
;      print,hlt,hln
      FOR mm=0,nmons-1 DO BEGIN
	subhi=fltarr(12,12) 
	FOR i=0,5 DO BEGIN
	  FOR j=0,5 DO BEGIN
	    scol=(i*2)
	    jrow=(j*2)
	    IF (hln EQ 360) THEN hln=0
	    subhi(scol:scol+1,jrow:jrow+1)=hiabs(hln,hlt,mm)
	    hlt=hlt+1
          ENDFOR
	  hlt=hlt-6
	  hln=hln+1
	ENDFOR
	hln=hln-6
	absarr(lnn,ltt,mm)=MEAN(subhi(1:11,1:11))
      ENDFOR
      hln=hln+5
;      stop
    ENDFOR
;    stop
    hln=0
    hlt=hlt+5
  ENDFOR
;  stop
  ; now flip and shift so its -90 to 90, -180 to 180
  FOR mm=0,nmons-1 DO BEGIN
    absarr(*,*,mm)=REVERSE(absarr(*,*,mm),2)
    absarr(*,*,mm)=SHIFT(absarr(*,*,mm),nlons/2,0)
  ENDFOR
  
  ; Output to netCDF
  wilma=NCDF_CREATE(indir+varlist[loo]+outfil,/clobber)
  tid=NCDF_DIMDEF(wilma,'time',nmons)
  charid=NCDF_DIMDEF(wilma, 'Character', 3)
  latid=NCDF_DIMDEF(wilma,'latitude',nlats)
  lonid=NCDF_DIMDEF(wilma,'longitude',nlons)
  
  timesvar=NCDF_VARDEF(wilma,'time',[tid],/SHORT)
  latsvar=NCDF_VARDEF(wilma,'latitude',[latid],/FLOAT)
  lonsvar=NCDF_VARDEF(wilma,'longitude',[lonid],/FLOAT)
  absvar=NCDF_VARDEF(wilma,'actuals',[lonid,latid,tid],/FLOAT)

  NCDF_ATTPUT,wilma,'time','standard_name','time'
  NCDF_ATTPUT,wilma,'time','long_name','time'
  NCDF_ATTPUT,wilma,'time','units','days since 1979-1-1 00:00:00'
  NCDF_ATTPUT,wilma,'time','axis','T'
  NCDF_ATTPUT,wilma,'time','calendar','gregorian'
  NCDF_ATTPUT,wilma,'time','start_year',styr
  NCDF_ATTPUT,wilma,'time','end_year',edyr
  NCDF_ATTPUT,wilma,'time','start_month',1
  NCDF_ATTPUT,wilma,'time','end_month',12

  NCDF_ATTPUT,wilma,'latitude','standard_name','latitude'
  NCDF_ATTPUT,wilma,'latitude','long_name','latitude'
  NCDF_ATTPUT,wilma,'latitude','units','degrees_north'
  NCDF_ATTPUT,wilma,'latitude','point_spacing','even'
  NCDF_ATTPUT,wilma,'latitude','axis','X'

  NCDF_ATTPUT,wilma,'longitude','standard_name','longitude'
  NCDF_ATTPUT,wilma,'longitude','long_name','longitude'
  NCDF_ATTPUT,wilma,'longitude','units','degrees_east'
  NCDF_ATTPUT,wilma,'longitude','point_spacing','even'
  NCDF_ATTPUT,wilma,'longitude','axis','X'

  NCDF_ATTPUT,wilma,absvar,'long_name','Monthly mean absolutes'
  NCDF_ATTPUT,wilma,absvar,'units',unitslist[loo]
  NCDF_ATTPUT,wilma,absvar,'axis','T'
  valid=WHERE(absarr NE -1.E+30, tc)
  IF tc GE 1 THEN BEGIN
    min_t=MIN(absarr(valid))
    max_t=MAX(absarr(valid))
    NCDF_ATTPUT,wilma,absvar,'valid_min',min_t(0)
    NCDF_ATTPUT,wilma,absvar,'valid_max',max_t(0)
  ENDIF
  NCDF_ATTPUT,wilma,absvar,'missing_value',-1.e+30
  NCDF_ATTPUT,wilma,absvar,'_FillValue',-1.e+30

  current_time=SYSTIME()

  NCDF_ATTPUT,wilma,/GLOBAL,'file_created',STRING(current_time)
  NCDF_ATTPUT,wilma,/GLOBAL,'description',"ERA-Interim monthly mean land surface "+varlist[loo]+" reanalysis product from 1979 onwards. "+$
                                         "regridded to 5by5 degree boxes by simple averaging (no smoothing."
  NCDF_ATTPUT,wilma,/GLOBAL,'title',"ERA-Interim monthly mean land surface "+varlist[loo]+" climate monitoring product from 1979 onwards."
  NCDF_ATTPUT,wilma,/GLOBAL,'institution',"ECMWF (regridded at Met Office Hadley Centre by Kate Willett"
  NCDF_ATTPUT,wilma,/GLOBAL,'history',"Updated "+STRING(current_time)
  NCDF_ATTPUT,wilma,/GLOBAL,'source',"http://apps.ecmwf.int/datasets/data/interim-full-daily/ to ERAMONTHLY_6hrly_1by1_month chunks"
  NCDF_ATTPUT,wilma,/GLOBAL,'comment'," "
  NCDF_ATTPUT,wilma,/GLOBAL,'reference',"NA"
  NCDF_ATTPUT,wilma,/GLOBAL,'version',"Last Download Feb 2018"
  NCDF_ATTPUT,wilma,/GLOBAL,'Conventions',"CF-1.0"

  NCDF_CONTROL,wilma,/ENDEF
  NCDF_VARPUT, wilma,timesvar,dayssince
  NCDF_VARPUT, wilma,latsvar, lats
  NCDF_VARPUT, wilma,lonsvar, lons
  NCDF_VARPUT, wilma,absvar,absarr

  NCDF_CLOSE,wilma

; now also create anomalies 81-10
; have a land masked and non-land masked version

  anomsall=make_array(nlons,nlats,nmons,/float,value=mdi)
  anomsland=make_array(nlons,nlats,nmons,/float,value=mdi)
  anomssea=make_array(nlons,nlats,nmons,/float,value=mdi)
  FOR lnn=0,nlons-1 DO BEGIN
    FOR ltt=0,nlats-1 DO BEGIN
      subarr=reform(absarr(lnn,ltt,*),12,nyrs)
      FOR mm=0,11 DO BEGIN
        ; NO MISSING DATA IN ERA!
	subarr(mm,*)=subarr(mm,*)-MEAN(subarr(mm,climst1-styr:climed1-styr))
      ENDFOR
      anomsall(lnn,ltt,*)=reform(subarr,nmons)  
    ENDFOR
  ENDFOR
  
  ; read in land mask and make masked arrays
  inn=NCDF_OPEN(indir+inlandmask)
  varid=NCDF_VARID(inn,'land_area_fraction')
  NCDF_VARGET,inn,varid,pct_land
  NCDF_CLOSE,inn
  pct_land=pct_land(*,*,0)
  
  FOR lnn=0,nlons-1 DO BEGIN
    FOR ltt=0,nlats-1 DO BEGIN
      IF (pct_land(lnn,ltt) GT 0) THEN anomsland(lnn,ltt,*)=anomsall(lnn,ltt,*)
      IF (pct_land(lnn,ltt) LT 0.9) THEN BEGIN
        print,lnn,ltt,pct_land(lnn,ltt)
	anomssea(lnn,ltt,*)=anomsall(lnn,ltt,*)
      ENDIF
    ENDFOR
  ENDFOR  

  ; output to netCDF
  wilma=NCDF_CREATE(indir+varlist[loo]+outfilA,/clobber)
  tid=NCDF_DIMDEF(wilma,'time',nmons)
  latid=NCDF_DIMDEF(wilma,'latitude',nlats)
  lonid=NCDF_DIMDEF(wilma,'longitude',nlons)
  
  timesvar=NCDF_VARDEF(wilma,'time',[tid],/SHORT)
  latsvar=NCDF_VARDEF(wilma,'latitude',[latid],/FLOAT)
  lonsvar=NCDF_VARDEF(wilma,'longitude',[lonid],/FLOAT)
  allvar=NCDF_VARDEF(wilma,'anomalies',[lonid,latid,tid],/FLOAT)
  landvar=NCDF_VARDEF(wilma,'anomalies_land',[lonid,latid,tid],/FLOAT)
  seavar=NCDF_VARDEF(wilma,'anomalies_sea',[lonid,latid,tid],/FLOAT)

  NCDF_ATTPUT,wilma,'time','standard_name','time'
  NCDF_ATTPUT,wilma,'time','long_name','time'
  NCDF_ATTPUT,wilma,'time','units','days since 1979-1-1 00:00:00'
  NCDF_ATTPUT,wilma,'time','axis','T'
  NCDF_ATTPUT,wilma,'time','calendar','gregorian'
  NCDF_ATTPUT,wilma,'time','start_year',styr
  NCDF_ATTPUT,wilma,'time','end_year',edyr
  NCDF_ATTPUT,wilma,'time','start_month',1
  NCDF_ATTPUT,wilma,'time','end_month',12

  NCDF_ATTPUT,wilma,'latitude','standard_name','latitude'
  NCDF_ATTPUT,wilma,'latitude','long_name','latitude'
  NCDF_ATTPUT,wilma,'latitude','units','degrees_north'
  NCDF_ATTPUT,wilma,'latitude','point_spacing','even'
  NCDF_ATTPUT,wilma,'latitude','axis','X'

  NCDF_ATTPUT,wilma,'longitude','standard_name','longitude'
  NCDF_ATTPUT,wilma,'longitude','long_name','longitude'
  NCDF_ATTPUT,wilma,'longitude','units','degrees_east'
  NCDF_ATTPUT,wilma,'longitude','point_spacing','even'
  NCDF_ATTPUT,wilma,'longitude','axis','X'

  NCDF_ATTPUT,wilma,allvar,'long_name','Monthly mean anomalies relative to 1981-2010'
  NCDF_ATTPUT,wilma,allvar,'units',unitslist[loo]
  NCDF_ATTPUT,wilma,allvar,'axis','T'
  valid=WHERE(anomsall NE -1.E+30, tc)
  IF tc GE 1 THEN BEGIN
    min_t=MIN(anomsall(valid))
    max_t=MAX(anomsall(valid))
    NCDF_ATTPUT,wilma,allvar,'valid_min',min_t(0)
    NCDF_ATTPUT,wilma,allvar,'valid_max',max_t(0)
  ENDIF
  NCDF_ATTPUT,wilma,allvar,'missing_value',-1.e+30
  NCDF_ATTPUT,wilma,allvar,'_FillValue',-1.e+30
  NCDF_ATTPUT,wilma,allvar,'reference_period','1981 to 2010'

  NCDF_ATTPUT,wilma,landvar,'long_name','Land only monthly mean anomalies relative to 1981-2010'
  NCDF_ATTPUT,wilma,landvar,'units',unitslist[loo]
  NCDF_ATTPUT,wilma,landvar,'axis','T'
  valid=WHERE(anomsland NE -1.E+30, tc)
  IF tc GE 1 THEN BEGIN
    min_t=MIN(anomsland(valid))
    max_t=MAX(anomsland(valid))
    NCDF_ATTPUT,wilma,landvar,'valid_min',min_t(0)
    NCDF_ATTPUT,wilma,landvar,'valid_max',max_t(0)
  ENDIF
  NCDF_ATTPUT,wilma,landvar,'missing_value',-1.e+30
  NCDF_ATTPUT,wilma,landvar,'_FillValue',-1.e+30
  NCDF_ATTPUT,wilma,landvar,'reference_period','1981 to 2010'

  NCDF_ATTPUT,wilma,seavar,'long_name','Sea only monthly mean anomalies relative to 1981-2010'
  NCDF_ATTPUT,wilma,seavar,'units',unitslist[loo]
  NCDF_ATTPUT,wilma,seavar,'axis','T'
  valid=WHERE(anomssea NE -1.E+30, tc)
  IF tc GE 1 THEN BEGIN
    min_t=MIN(anomssea(valid))
    max_t=MAX(anomssea(valid))
    NCDF_ATTPUT,wilma,seavar,'valid_min',min_t(0)
    NCDF_ATTPUT,wilma,seavar,'valid_max',max_t(0)
  ENDIF
  NCDF_ATTPUT,wilma,seavar,'missing_value',-1.e+30
  NCDF_ATTPUT,wilma,seavar,'_FillValue',-1.e+30
  NCDF_ATTPUT,wilma,seavar,'reference_period','1981 to 2010'

  current_time=SYSTIME()

  NCDF_ATTPUT,wilma,/GLOBAL,'file_created',STRING(current_time)
  NCDF_ATTPUT,wilma,/GLOBAL,'description',"ERA-Interim monthly mean land surface "+varlist[loo]+" reanalysis product from 1979 onwards. "+$
                                         "regridded to 5by5 degree boxes by simple averaging (no smoothing."
  NCDF_ATTPUT,wilma,/GLOBAL,'title',"ERA-Interim monthly mean land surface "+varlist[loo]+" climate monitoring product from 1979 onwards."
  NCDF_ATTPUT,wilma,/GLOBAL,'institution',"ECMWF (regridded at Met Office Hadley Centre by Kate Willett"
  NCDF_ATTPUT,wilma,/GLOBAL,'history',"Updated "+STRING(current_time)
  NCDF_ATTPUT,wilma,/GLOBAL,'source',"http://apps.ecmwf.int/datasets/data/interim-full-daily/ to ERAMONTHLY_6hrly_1by1_ month chunks"
  NCDF_ATTPUT,wilma,/GLOBAL,'comment'," "
  NCDF_ATTPUT,wilma,/GLOBAL,'reference',"NA"
  NCDF_ATTPUT,wilma,/GLOBAL,'version',"Last Download Feb 2018"
  NCDF_ATTPUT,wilma,/GLOBAL,'Conventions',"CF-1.0"

  NCDF_CONTROL,wilma,/ENDEF
  NCDF_VARPUT, wilma,timesvar,dayssince
  NCDF_VARPUT, wilma,latsvar, lats
  NCDF_VARPUT, wilma,lonsvar, lons
  NCDF_VARPUT, wilma,allvar,anomsall
  NCDF_VARPUT, wilma,landvar,anomsland
  NCDF_VARPUT, wilma,seavar,anomssea

  NCDF_CLOSE,wilma

; now also create anomalies 79-05 (NOT IDEAL - BETTER TO MOVE 81-10 TO MATCH HADISDH OVER 79-05 PERIOD?
; have a land masked and non-land masked version

  anomsall=make_array(nlons,nlats,nmons,/float,value=mdi)
  anomsland=make_array(nlons,nlats,nmons,/float,value=mdi)
  anomssea=make_array(nlons,nlats,nmons,/float,value=mdi)
  FOR lnn=0,nlons-1 DO BEGIN
    FOR ltt=0,nlats-1 DO BEGIN
      subarr=reform(absarr(lnn,ltt,*),12,nyrs)
      FOR mm=0,11 DO BEGIN
        ; NO MISSING DATA IN ERA!
	subarr(mm,*)=subarr(mm,*)-MEAN(subarr(mm,climst2-styr:climed2-styr))
      ENDFOR
      anomsall(lnn,ltt,*)=reform(subarr,nmons)  
    ENDFOR
  ENDFOR
  
  ; already read in land mask and make masked arrays
  
  FOR lnn=0,nlons-1 DO BEGIN
    FOR ltt=0,nlats-1 DO BEGIN
      IF (pct_land(lnn,ltt) GT 0) THEN anomsland(lnn,ltt,*)=anomsall(lnn,ltt,*)
      IF (pct_land(lnn,ltt) LT 0.9) THEN BEGIN
        print,lnn,ltt,pct_land(lnn,ltt)
	anomssea(lnn,ltt,*)=anomsall(lnn,ltt,*)
      ENDIF
    ENDFOR
  ENDFOR  

  ; output to netCDF
  wilma=NCDF_CREATE(indir+varlist[loo]+outfilA7905,/clobber)
  tid=NCDF_DIMDEF(wilma,'time',nmons)
  latid=NCDF_DIMDEF(wilma,'latitude',nlats)
  lonid=NCDF_DIMDEF(wilma,'longitude',nlons)
  
  timesvar=NCDF_VARDEF(wilma,'time',[tid],/SHORT)
  latsvar=NCDF_VARDEF(wilma,'latitude',[latid],/FLOAT)
  lonsvar=NCDF_VARDEF(wilma,'longitude',[lonid],/FLOAT)
  allvar=NCDF_VARDEF(wilma,'anomalies',[lonid,latid,tid],/FLOAT)
  landvar=NCDF_VARDEF(wilma,'anomalies_land',[lonid,latid,tid],/FLOAT)
  seavar=NCDF_VARDEF(wilma,'anomalies_sea',[lonid,latid,tid],/FLOAT)

  NCDF_ATTPUT,wilma,'time','standard_name','time'
  NCDF_ATTPUT,wilma,'time','long_name','time'
  NCDF_ATTPUT,wilma,'time','units','days since 1979-1-1 00:00:00'
  NCDF_ATTPUT,wilma,'time','axis','T'
  NCDF_ATTPUT,wilma,'time','calendar','gregorian'
  NCDF_ATTPUT,wilma,'time','start_year',styr
  NCDF_ATTPUT,wilma,'time','end_year',edyr
  NCDF_ATTPUT,wilma,'time','start_month',1
  NCDF_ATTPUT,wilma,'time','end_month',12

  NCDF_ATTPUT,wilma,'latitude','standard_name','latitude'
  NCDF_ATTPUT,wilma,'latitude','long_name','latitude'
  NCDF_ATTPUT,wilma,'latitude','units','degrees_north'
  NCDF_ATTPUT,wilma,'latitude','point_spacing','even'
  NCDF_ATTPUT,wilma,'latitude','axis','X'

  NCDF_ATTPUT,wilma,'longitude','standard_name','longitude'
  NCDF_ATTPUT,wilma,'longitude','long_name','longitude'
  NCDF_ATTPUT,wilma,'longitude','units','degrees_east'
  NCDF_ATTPUT,wilma,'longitude','point_spacing','even'
  NCDF_ATTPUT,wilma,'longitude','axis','X'

  NCDF_ATTPUT,wilma,allvar,'long_name','Monthly mean anomalies relative to 1979-2005'
  NCDF_ATTPUT,wilma,allvar,'units',unitslist[loo]
  NCDF_ATTPUT,wilma,allvar,'axis','T'
  valid=WHERE(anomsall NE -1.E+30, tc)
  IF tc GE 1 THEN BEGIN
    min_t=MIN(anomsall(valid))
    max_t=MAX(anomsall(valid))
    NCDF_ATTPUT,wilma,allvar,'valid_min',min_t(0)
    NCDF_ATTPUT,wilma,allvar,'valid_max',max_t(0)
  ENDIF
  NCDF_ATTPUT,wilma,allvar,'missing_value',-1.e+30
  NCDF_ATTPUT,wilma,allvar,'_FillValue',-1.e+30
  NCDF_ATTPUT,wilma,allvar,'reference_period','1979 to 2005'

  NCDF_ATTPUT,wilma,landvar,'long_name','Land only monthly mean anomalies relative to 1979-2005'
  NCDF_ATTPUT,wilma,landvar,'units',unitslist[loo]
  NCDF_ATTPUT,wilma,landvar,'axis','T'
  valid=WHERE(anomsland NE -1.E+30, tc)
  IF tc GE 1 THEN BEGIN
    min_t=MIN(anomsland(valid))
    max_t=MAX(anomsland(valid))
    NCDF_ATTPUT,wilma,landvar,'valid_min',min_t(0)
    NCDF_ATTPUT,wilma,landvar,'valid_max',max_t(0)
  ENDIF
  NCDF_ATTPUT,wilma,landvar,'missing_value',-1.e+30
  NCDF_ATTPUT,wilma,landvar,'_FillValue',-1.e+30
  NCDF_ATTPUT,wilma,landvar,'reference_period','1979 to 2005'

  NCDF_ATTPUT,wilma,seavar,'long_name','Sea only monthly mean anomalies relative to 1979-2005'
  NCDF_ATTPUT,wilma,seavar,'units',unitslist[loo]
  NCDF_ATTPUT,wilma,seavar,'axis','T'
  valid=WHERE(anomssea NE -1.E+30, tc)
  IF tc GE 1 THEN BEGIN
    min_t=MIN(anomssea(valid))
    max_t=MAX(anomssea(valid))
    NCDF_ATTPUT,wilma,seavar,'valid_min',min_t(0)
    NCDF_ATTPUT,wilma,seavar,'valid_max',max_t(0)
  ENDIF
  NCDF_ATTPUT,wilma,seavar,'missing_value',-1.e+30
  NCDF_ATTPUT,wilma,seavar,'_FillValue',-1.e+30
  NCDF_ATTPUT,wilma,seavar,'reference_period','1979 to 2005'

  current_time=SYSTIME()

  NCDF_ATTPUT,wilma,/GLOBAL,'file_created',STRING(current_time)
  NCDF_ATTPUT,wilma,/GLOBAL,'description',"ERA-Interim monthly mean land surface "+varlist[loo]+" reanalysis product from 1979 onwards. "+$
                                         "regridded to 5by5 degree boxes by simple averaging (no smoothing."
  NCDF_ATTPUT,wilma,/GLOBAL,'title',"ERA-Interim monthly mean land surface "+varlist[loo]+" climate monitoring product from 1979 onwards."
  NCDF_ATTPUT,wilma,/GLOBAL,'institution',"ECMWF (regridded at Met Office Hadley Centre by Kate Willett"
  NCDF_ATTPUT,wilma,/GLOBAL,'history',"Updated "+STRING(current_time)
  NCDF_ATTPUT,wilma,/GLOBAL,'source',"http://apps.ecmwf.int/datasets/data/interim-full-daily/ to ERAMONTHLY_6hrly_1by1_ month chunks"
  NCDF_ATTPUT,wilma,/GLOBAL,'comment'," "
  NCDF_ATTPUT,wilma,/GLOBAL,'reference',"NA"
  NCDF_ATTPUT,wilma,/GLOBAL,'version',"Last Download Feb 2018"
  NCDF_ATTPUT,wilma,/GLOBAL,'Conventions',"CF-1.0"

  NCDF_CONTROL,wilma,/ENDEF
  NCDF_VARPUT, wilma,timesvar,dayssince
  NCDF_VARPUT, wilma,latsvar, lats
  NCDF_VARPUT, wilma,lonsvar, lons
  NCDF_VARPUT, wilma,allvar,anomsall
  NCDF_VARPUT, wilma,landvar,anomsland
  NCDF_VARPUT, wilma,seavar,anomssea

  NCDF_CLOSE,wilma
  
ENDFOR


END
