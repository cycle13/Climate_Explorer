pro sortout_CRUTEM_fields_MAY2015

; program to read in CRUTEM anomalies and station/sampling uncertainty fields
; reads in HadISDH.landT anomalies too
; outputs a single netCDF in HadISDH stylee with anomaliies and uncertainties 
; BUT the anomalies have been renormalised over the HadISDH period (1976-2005) (no need to do anything to uncertainties!)
; AND there is another file where CRUTEM is masked to the coverage of HadISDH
; This is so Robert can compute HadISDH like regional uncertainties
; Looks like CRUTEM uncs are 1sigma so I've multiplied by two to match with HadISDH 2sigma unc.

;----------------------------------------------------------------------------
inCanom='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/CRUTEM.4.3.0.0.anomalies.nc'
inCsamp='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/CRUTEM.4.3.0.0.sampling_error.nc'
inCstat='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/CRUTEM.4.3.0.0.station_error.nc'

inHanom='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf.nc'

outfil='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/CRUTEM.4.3.0.0.anomalies_uncertainties19762005clim_19732014.nc'
outfilMASK='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/CRUTEM.4.3.0.0.MASKEDanomalies_uncertainties19762005clim_19732014.nc'

mdi=-1e30

stH=1973
edH=2014
stC=1850

nyrs=(edH+1)-stH
nmons=nyrs*12

dayssince=findgen(nmons)
stpoint=JULDAY(1,1,1973)

; array of months in days since Jan 1st 1979 including Jan 2015 (nmons+1)
monthies=TIMEGEN(nmons+1,START=(stpoint),UNITS="Months")-stpoint ; an array of floats describing each month in days since Jan 1st 1973
monthies=ROUND(monthies*100.)/100.	; nicer numbers

; Now need to use the mid-point of each month.
FOR mm=0,nmons-1 DO BEGIN
  dayssince(mm)=monthies(mm)+(monthies(mm+1)-monthies(mm))/2.
ENDFOR

climst=1976
climed=2005

nlats=36
nlons=72

cruanomsFULL=0
cruanomsFULLSAMP=0
cruanomsFULLSTAT=0
cruanomsSHORT=make_array(nlons,nlats,nmons,/float,value=mdi)
cruanomsSHORTMASK=make_array(nlons,nlats,nmons,/float,value=mdi)
cruanomsSHORTSAMP=make_array(nlons,nlats,nmons,/float,value=mdi)
cruanomsSHORTSTAT=make_array(nlons,nlats,nmons,/float,value=mdi)
cruanomsSHORTSAMPMASK=make_array(nlons,nlats,nmons,/float,value=mdi)
cruanomsSHORTSTATMASK=make_array(nlons,nlats,nmons,/float,value=mdi)

hadisdh=0
;------------------------------------------------------------------------------
; read in HadISDH and lons, lats and times
inn=NCDF_OPEN(inHanom)
varid=NCDF_VARID(inn,'t_anoms')
latid=NCDF_VARID(inn,'latitude')
lonid=NCDF_VARID(inn,'longitude')
timid=NCDF_VARID(inn,'time')
NCDF_VARGET,inn,varid,hadisdh
NCDF_VARGET,inn,latid,lats
NCDF_VARGET,inn,lonid,lons
NCDF_VARGET,inn,timid,times
NCDF_CLOSE,inn

; read in CRUTEMS
inn=NCDF_OPEN(inCanom)
varid=NCDF_VARID(inn,'temperature_anomaly')
NCDF_VARGET,inn,varid,cruanomsFULL
NCDF_CLOSE,inn

inn=NCDF_OPEN(inCsamp)
varid=NCDF_VARID(inn,'standard_error')
NCDF_VARGET,inn,varid,cruanomsFULLSAMP
NCDF_CLOSE,inn

inn=NCDF_OPEN(inCstat)
varid=NCDF_VARID(inn,'standard_error')
NCDF_VARGET,inn,varid,cruanomsFULLSTAT
NCDF_CLOSE,inn

; shorten to HadISDH time points only
cruanomsSHORT(*,*,*)=cruanomsFULL(*,*,(stH-stC)*12:(((edH+1)-stC)*12)-1)
cruanomsSHORTSAMP(*,*,*)=cruanomsFULLSAMP(*,*,(stH-stC)*12:(((edH+1)-stC)*12)-1)
cruanomsSHORTSTAT(*,*,*)=cruanomsFULLSTAT(*,*,(stH-stC)*12:(((edH+1)-stC)*12)-1)

cruanomsSHORTMASK(*,*,*)=cruanomsSHORT
cruanomsSHORTSAMPMASK(*,*,*)=cruanomsSHORTSAMP
cruanomsSHORTSTATMASK(*,*,*)=cruanomsSHORTSTAT

; renormalise to new clim period
FOR ltt=0,nlats-1 DO BEGIN
  FOR lnn=0,nlons-1 DO BEGIN
    subarr=REFORM(cruanomsSHORT(lnn,ltt,*),12,nyrs)    
    FOR mm=0,11 DO BEGIN
      climsie=subarr(mm,climst-stH:(climed+1)-stH)
      gotclims=where(climsie GT mdi,count)
      IF (count GT 15) THEN BEGIN
        gots=where(subarr(mm,*) GT mdi,count)
	subarr(mm,gots)=subarr(mm,gots)-MEAN(climsie(gotclims))
      ENDIF ELSE BEGIN
        subarr(mm,*)=mdi
      ENDELSE
    ENDFOR
    cruanomsSHORT(lnn,ltt,*)=REFORM(subarr,nmons)
  ENDFOR
ENDFOR

; Mask to HadISDH
FOR ltt=0,nlats-1 DO BEGIN
  FOR lnn=0,nlons-1 DO BEGIN
    gots=WHERE(hadisdh(lnn,ltt,*) EQ mdi,count)
    IF (count GT 0) THEN BEGIN
      ;print,'GOT SOME MISSINGS!'
      cruanomsSHORTMASK(lnn,ltt,gots)=mdi
      cruanomsSHORTSAMPMASK(lnn,ltt,gots)=mdi
      cruanomsSHORTSTATMASK(lnn,ltt,gots)=mdi      
    ENDIF 
  ENDFOR
ENDFOR

; multple CRUTEM uncertainties by two to give 2sigma uncertainty
gots=where(cruanomsSHORTSAMP GT mdi)
cruanomsSHORTSAMP(gots)=cruanomsSHORTSAMP(gots)*2.
gots=where(cruanomsSHORTSTAT GT mdi)
cruanomsSHORTSTAT(gots)=cruanomsSHORTSTAT(gots)*2.
gots=where(cruanomsSHORTSAMPMASK GT mdi)
cruanomsSHORTSAMPMASK(gots)=cruanomsSHORTSAMPMASK(gots)*2.
gots=where(cruanomsSHORTSTATMASK GT mdi)
cruanomsSHORTSTATMASK(gots)=cruanomsSHORTSTATMASK(gots)*2.

; Output NetCDFs
; Output to netCDF
wilma=NCDF_CREATE(outfil,/clobber)
tid=NCDF_DIMDEF(wilma,'time',nmons)
latid=NCDF_DIMDEF(wilma,'latitude',nlats)
lonid=NCDF_DIMDEF(wilma,'longitude',nlons)

timesvar=NCDF_VARDEF(wilma,'time',[tid],/SHORT)
latsvar=NCDF_VARDEF(wilma,'latitude',[latid],/FLOAT)
lonsvar=NCDF_VARDEF(wilma,'longitude',[lonid],/FLOAT)
absvar=NCDF_VARDEF(wilma,'t_anoms',[lonid,latid,tid],/FLOAT)
sampvar=NCDF_VARDEF(wilma,'t_samplingerr',[lonid,latid,tid],/FLOAT)
statvar=NCDF_VARDEF(wilma,'t_stationerr',[lonid,latid,tid],/FLOAT)

NCDF_ATTPUT,wilma,'time','standard_name','time'
NCDF_ATTPUT,wilma,'time','long_name','time'
NCDF_ATTPUT,wilma,'time','units','days since 1979-1-1 00:00:00'
NCDF_ATTPUT,wilma,'time','axis','T'
NCDF_ATTPUT,wilma,'time','calendar','gregorian'
NCDF_ATTPUT,wilma,'time','start_year',stH
NCDF_ATTPUT,wilma,'time','end_year',edH
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

NCDF_ATTPUT,wilma,absvar,'long_name','Monthly mean anomalies relative to 1976 to 2005'
NCDF_ATTPUT,wilma,absvar,'units','deg C'
NCDF_ATTPUT,wilma,absvar,'axis','T'
valid=WHERE(cruanomsSHORT NE -1.E+30, tc)
IF tc GE 1 THEN BEGIN
  min_t=MIN(cruanomsSHORT(valid))
  max_t=MAX(cruanomsSHORT(valid))
  NCDF_ATTPUT,wilma,absvar,'valid_min',min_t(0)
  NCDF_ATTPUT,wilma,absvar,'valid_max',max_t(0)
ENDIF
NCDF_ATTPUT,wilma,absvar,'missing_value',-1.e+30
NCDF_ATTPUT,wilma,absvar,'_FillValue',-1.e+30

NCDF_ATTPUT,wilma,sampvar,'long_name','2 sigma sampling uncertainty'
NCDF_ATTPUT,wilma,sampvar,'units','deg C'
NCDF_ATTPUT,wilma,sampvar,'axis','T'
valid=WHERE(cruanomsSHORTSAMP NE -1.E+30, tc)
IF tc GE 1 THEN BEGIN
  min_t=MIN(cruanomsSHORTSAMP(valid))
  max_t=MAX(cruanomsSHORTSAMP(valid))
  NCDF_ATTPUT,wilma,sampvar,'valid_min',min_t(0)
  NCDF_ATTPUT,wilma,sampvar,'valid_max',max_t(0)
ENDIF
NCDF_ATTPUT,wilma,sampvar,'missing_value',-1.e+30
NCDF_ATTPUT,wilma,sampvar,'_FillValue',-1.e+30

NCDF_ATTPUT,wilma,statvar,'long_name','2 sigma station uncertainty'
NCDF_ATTPUT,wilma,statvar,'units','deg C'
NCDF_ATTPUT,wilma,statvar,'axis','T'
valid=WHERE(cruanomsSHORTSTAT NE -1.E+30, tc)
IF tc GE 1 THEN BEGIN
  min_t=MIN(cruanomsSHORTSTAT(valid))
  max_t=MAX(cruanomsSHORTSTAT(valid))
  NCDF_ATTPUT,wilma,statvar,'valid_min',min_t(0)
  NCDF_ATTPUT,wilma,statvar,'valid_max',max_t(0)
ENDIF
NCDF_ATTPUT,wilma,statvar,'missing_value',-1.e+30
NCDF_ATTPUT,wilma,statvar,'_FillValue',-1.e+30

current_time=SYSTIME()

NCDF_ATTPUT,wilma,/GLOBAL,'file_created',STRING(current_time)
NCDF_ATTPUT,wilma,/GLOBAL,'description',"CRUTEM4.3.0.0 monthly mean anomalies renormalised to HadISDH period (1976-2005) with uncertainties"
NCDF_ATTPUT,wilma,/GLOBAL,'title',"CRUTEM4.3.0.0. surface temperature from 1973 onwards"
NCDF_ATTPUT,wilma,/GLOBAL,'institution',"MOHC (reformed at Met Office Hadley Centre by Kate Willett)"
NCDF_ATTPUT,wilma,/GLOBAL,'history',"Updated "+STRING(current_time)
NCDF_ATTPUT,wilma,/GLOBAL,'source',"www.metoffice.gov/hadobs/crutem4/"
NCDF_ATTPUT,wilma,/GLOBAL,'comment'," "
NCDF_ATTPUT,wilma,/GLOBAL,'reference',"NA"
NCDF_ATTPUT,wilma,/GLOBAL,'version',"Last Download Feb 2015"
NCDF_ATTPUT,wilma,/GLOBAL,'Conventions',"CF-1.0"

NCDF_CONTROL,wilma,/ENDEF
NCDF_VARPUT, wilma,timesvar,dayssince
NCDF_VARPUT, wilma,latsvar, lats
NCDF_VARPUT, wilma,lonsvar, lons
NCDF_VARPUT, wilma,absvar,cruanomsSHORT
NCDF_VARPUT, wilma,sampvar,cruanomsSHORTSAMP
NCDF_VARPUT, wilma,statvar,cruanomsSHORTSTAT

NCDF_CLOSE,wilma


wilma=NCDF_CREATE(outfilMASK,/clobber)
tid=NCDF_DIMDEF(wilma,'time',nmons)
latid=NCDF_DIMDEF(wilma,'latitude',nlats)
lonid=NCDF_DIMDEF(wilma,'longitude',nlons)

timesvar=NCDF_VARDEF(wilma,'time',[tid],/SHORT)
latsvar=NCDF_VARDEF(wilma,'latitude',[latid],/FLOAT)
lonsvar=NCDF_VARDEF(wilma,'longitude',[lonid],/FLOAT)
absvar=NCDF_VARDEF(wilma,'t_anoms',[lonid,latid,tid],/FLOAT)
sampvar=NCDF_VARDEF(wilma,'t_samplingerr',[lonid,latid,tid],/FLOAT)
statvar=NCDF_VARDEF(wilma,'t_stationerr',[lonid,latid,tid],/FLOAT)

NCDF_ATTPUT,wilma,'time','standard_name','time'
NCDF_ATTPUT,wilma,'time','long_name','time'
NCDF_ATTPUT,wilma,'time','units','days since 1979-1-1 00:00:00'
NCDF_ATTPUT,wilma,'time','axis','T'
NCDF_ATTPUT,wilma,'time','calendar','gregorian'
NCDF_ATTPUT,wilma,'time','start_year',stH
NCDF_ATTPUT,wilma,'time','end_year',edH
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

NCDF_ATTPUT,wilma,absvar,'long_name','Monthly mean anomalies relative to 1976 to 2005 masked to HadISDH coverage'
NCDF_ATTPUT,wilma,absvar,'units','deg C'
NCDF_ATTPUT,wilma,absvar,'axis','T'
valid=WHERE(cruanomsSHORTMASK NE -1.E+30, tc)
IF tc GE 1 THEN BEGIN
  min_t=MIN(cruanomsSHORTMASK(valid))
  max_t=MAX(cruanomsSHORTMASK(valid))
  NCDF_ATTPUT,wilma,absvar,'valid_min',min_t(0)
  NCDF_ATTPUT,wilma,absvar,'valid_max',max_t(0)
ENDIF
NCDF_ATTPUT,wilma,absvar,'missing_value',-1.e+30
NCDF_ATTPUT,wilma,absvar,'_FillValue',-1.e+30

NCDF_ATTPUT,wilma,sampvar,'long_name','2 sigma sampling uncertainty'
NCDF_ATTPUT,wilma,sampvar,'units','deg C'
NCDF_ATTPUT,wilma,sampvar,'axis','T'
valid=WHERE(cruanomsSHORTSAMPMASK NE -1.E+30, tc)
IF tc GE 1 THEN BEGIN
  min_t=MIN(cruanomsSHORTSAMPMASK(valid))
  max_t=MAX(cruanomsSHORTSAMPMASK(valid))
  NCDF_ATTPUT,wilma,sampvar,'valid_min',min_t(0)
  NCDF_ATTPUT,wilma,sampvar,'valid_max',max_t(0)
ENDIF
NCDF_ATTPUT,wilma,sampvar,'missing_value',-1.e+30
NCDF_ATTPUT,wilma,sampvar,'_FillValue',-1.e+30

NCDF_ATTPUT,wilma,statvar,'long_name','2 sigma station uncertainty'
NCDF_ATTPUT,wilma,statvar,'units','deg C'
NCDF_ATTPUT,wilma,statvar,'axis','T'
valid=WHERE(cruanomsSHORTSTATMASK NE -1.E+30, tc)
IF tc GE 1 THEN BEGIN
  min_t=MIN(cruanomsSHORTSTATMASK(valid))
  max_t=MAX(cruanomsSHORTSTATMASK(valid))
  NCDF_ATTPUT,wilma,statvar,'valid_min',min_t(0)
  NCDF_ATTPUT,wilma,statvar,'valid_max',max_t(0)
ENDIF
NCDF_ATTPUT,wilma,statvar,'missing_value',-1.e+30
NCDF_ATTPUT,wilma,statvar,'_FillValue',-1.e+30

current_time=SYSTIME()

NCDF_ATTPUT,wilma,/GLOBAL,'file_created',STRING(current_time)
NCDF_ATTPUT,wilma,/GLOBAL,'description',"CRUTEM4.3.0.0 monthly mean anomalies renormalised to HadISDH period (1976-2005) with uncertainties and masked to HadISDH coverage"
NCDF_ATTPUT,wilma,/GLOBAL,'title',"CRUTEM4.3.0.0. surface temperature from 1973 onwards masked to HadISDH coverage"
NCDF_ATTPUT,wilma,/GLOBAL,'institution',"MOHC (reformed at Met Office Hadley Centre by Kate Willett)"
NCDF_ATTPUT,wilma,/GLOBAL,'history',"Updated "+STRING(current_time)
NCDF_ATTPUT,wilma,/GLOBAL,'source',"www.metoffice.gov/hadobs/crutem4/"
NCDF_ATTPUT,wilma,/GLOBAL,'comment'," "
NCDF_ATTPUT,wilma,/GLOBAL,'reference',"NA"
NCDF_ATTPUT,wilma,/GLOBAL,'version',"Last Download Feb 2015"
NCDF_ATTPUT,wilma,/GLOBAL,'Conventions',"CF-1.0"

NCDF_CONTROL,wilma,/ENDEF
NCDF_VARPUT, wilma,timesvar,dayssince
NCDF_VARPUT, wilma,latsvar, lats
NCDF_VARPUT, wilma,lonsvar, lons
NCDF_VARPUT, wilma,absvar,cruanomsSHORTMASK
NCDF_VARPUT, wilma,sampvar,cruanomsSHORTSAMPMASK
NCDF_VARPUT, wilma,statvar,cruanomsSHORTSTATMASK

NCDF_CLOSE,wilma

end
