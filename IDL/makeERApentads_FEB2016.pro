pro makeERApentads_FEB2016,param

; param is the variable you want to work with
; it can be: 'q2m','e2m','tw2m','t2m','td2m','rh2m','dpd2m','p2m'

; Read in 6 hourly ERA files
; Convert to other variables
; Average to pentad
; create pentad climatology
; create pentad anomalies

; MARINE OPTION:
; regrid to lie on 180 (89.5 to -89.5) latitudes rather than 181 (90 to -90)
;  ; regrid to 0.5 degree grids then re-average over 1 degree boxes
; shift to lie -180 to 180 longitude rather than 0 to 360
; Save

; USE calc_evap
; USE testleap

;************************************************

; File and vars set up

MARINE=1 ; 0 if normal, 1 if regrid and shift

mdi=-1e30

styr=1979	; 1979
edyr=2015	; 2013
clmst=1981
clmed=2010
nyrs=(edyr+1)-styr
nptds=nyrs*73
ndecs=5
declengths=[10,10,10,6,1]
int_ptds=indgen(nptds)
nleaps=9	;1980,1984,1988,1992,1996,2000,2004,2008,2012
;ndays=LONG((nyrs*365)+nleaps)
;nhrs=LONG(ndays*4)	;51136

nlons=360
nlats=181
marlats=180

;indir='/data/local/hadkw/HADCRUH2/UPDATE2015/OTHERDATA/'
;outdir='/data/local/hadkw/HADCRUH2/UPDATE2015/OTHERDATA/'
indir='/project/hadobs2/hadisdh/land/UPDATE2015/OTHERDATA/'
outdir='/project/hadobs2/hadisdh/land/UPDATE2015/OTHERDATA/'
infilT='ERAINTERIM_drybulbT_6hr_1by1_'
infilTd='ERAINTERIM_dewpointT_6hr_1by1_'
infilP='ERAINTERIM_sp_6hr_1by1_'
indecs=['1979010119881231','1989010119981231','1999010120081231','2009010120141231','2015010120151231']

;myparams=['q2m','e2m','tw2m','t2m','td2m','rh2m','dpd2m','p2m']
;IF (where(myparams EQ param) LT 0) THEN stop,"No param set"

if (MARINE EQ 0) THEN $
    outfile=param+'_pentad_1by1_ERA-Interim_data_'+strcompress(styr,/remove_all)+strcompress(edyr,/remove_all)+'.nc' ELSE $
    outfile=param+'_pentad_1by1marine_ERA-Interim_data_'+strcompress(styr,/remove_all)+strcompress(edyr,/remove_all)+'.nc'
    
pentad=make_array(nlons,nlats,nptds,type=FLOAT,value=mdi)

;************************************************
; Loop through each year, read in correct file
county=0	; counter through the years
; Make the pentad day times on the fly
moo=999
startee=JULDAY(1,1,styr)	; midday on 1st Jan 1979
FOR yy=0,nyrs-1 DO BEGIN
  isleap=0
  IF (testleap(yy+styr) EQ 0.0) THEN isleap=1 
  FOR pp=0,72 DO BEGIN
    ; time point creation
    IF (moo(0) EQ 999) THEN BEGIN; not yet initialised
      moo=TIMEGEN(1,UNIT='DAYS',start=JULDAY(1,1,styr))
      moolast=moo
    ENDIF ELSE BEGIN
      IF (isleap EQ 0) OR (pp NE 12) THEN BEGIN
        moolast=moolast+5
        moo=[moo,moolast]
      ENDIF ELSE BEGIN ; the start point after the 6 day leap year pentad
        moolast=moolast+6
        moo=[moo,moolast]
      ENDELSE
    ENDELSE
    print,yy+styr,moolast(0),pp,isleap
  ENDFOR
ENDFOR

; make times zero'd to 1979-1-1 that are half way through pentad
goomoo=(moo-startee)+2.5

FOR dd=0,ndecs-1 DO BEGIN
  sttim=0
  edtim=-1
  FOR yy=0,declengths(dd)-1 DO BEGIN
    
    ; get the correct pointers
    yrr=county+styr
    
    IF (yrr GT edyr) THEN break	; done all years!
    
    IF (testleap(yrr) EQ 0.0) THEN BEGIN
      HrTot=366*4
      sttim=edtim+1
      edtim=edtim+366*4
      pointarr=((indgen(73)*5)+5)*4  ; a pointer for the last point (hr) in each pentad (will need to -1 for indexing)
      pointarr(11:72)=pointarr(11:72)+4 ; added extra for the leap year
    ENDIF ELSE BEGIN
      HrTot=365*4
      sttim=edtim+1
      edtim=edtim+365*4
      pointarr=((indgen(73)*5)+5)*4  ; a pointer for the last point (hr) in each pentad (will need to -1 for indexing)        
    ENDELSE
    
    print,yrr,sttim,edtim
    
    IF (param NE 'p2m') THEN BEGIN
      ; get the correct filename and read in the year
      inn=NCDF_OPEN(indir+infilT+indecs(dd)+'.nc')
      varid=NCDF_VARID(inn,'t2m')
      latid=NCDF_VARID(inn,'latitude')
      lonid=NCDF_VARID(inn,'longitude')
    ; Only read in the data we need to save memory
      NCDF_VARGET,inn,varid,t,COUNT=[360,181,HrTot],STRIDE=[1,1,1],OFFSET=[0,0,sttim]
      NCDF_VARGET,inn,latid,lats
      NCDF_VARGET,inn,lonid,lons
      NCDF_ATTGET,inn,varid,'scale_factor',scaleet
      NCDF_ATTGET,inn,varid,'add_offset',offseteet
      NCDF_CLOSE,inn
      t=(((t*scaleet)+offseteet)-273.15)
      print,'Read T'
;    stop
    ENDIF
    
    IF ((param NE 't2m') AND (param NE 'p2m')) THEN BEGIN
      inn=NCDF_OPEN(indir+infilTd+indecs(dd)+'.nc')
      varid=NCDF_VARID(inn,'d2m')
      latid=NCDF_VARID(inn,'latitude')
      lonid=NCDF_VARID(inn,'longitude')
      NCDF_VARGET,inn,varid,d,COUNT=[360,181,HrTot],STRIDE=[1,1,1],OFFSET=[0,0,sttim]
      NCDF_VARGET,inn,latid,lats
      NCDF_VARGET,inn,lonid,lons
      NCDF_ATTGET,inn,varid,'scale_factor',scaleed
      NCDF_ATTGET,inn,varid,'add_offset',offseteed
      NCDF_CLOSE,inn
      d=(((d*scaleed)+offseteed)-273.15)
      print,'Read Td'
    ENDIF
    
    IF ((param NE 't2m') AND (param NE 'td2m')) THEN BEGIN
      inn=NCDF_OPEN(indir+infilP+indecs(dd)+'.nc')
      varid=NCDF_VARID(inn,'sp')
      latid=NCDF_VARID(inn,'latitude')
      lonid=NCDF_VARID(inn,'longitude')
      NCDF_VARGET,inn,varid,p,COUNT=[360,181,HrTot],STRIDE=[1,1,1],OFFSET=[0,0,sttim]
      NCDF_VARGET,inn,latid,lats
      NCDF_VARGET,inn,lonid,lons
      NCDF_ATTGET,inn,varid,'scale_factor',scaleep
      NCDF_ATTGET,inn,varid,'add_offset',offseteep
      NCDF_CLOSE,inn
      p=(((p*scaleep)+offseteep)/100.)
      print,'Read P'
    ENDIF  
    
    FOR ltt=0,nlats-1 DO BEGIN
      FOR lnn=0,nlons-1 DO BEGIN
        IF ((param NE 'td2m') AND (param NE 'p2m')) THEN tf=t(lnn,ltt,*)
        IF ((param NE 't2m') AND (param NE 'p2m')) THEN df=d(lnn,ltt,*)
        IF ((param NE 't2m') AND (param NE 'td2m')) THEN pf=p(lnn,ltt,*)
        IF ((param EQ 'e2m') OR (param EQ 'q2m') OR (param EQ 'rh2m') OR (param EQ 'tw2m')) THEN BEGIN
	  e=calc_evap(df,pf)
          esat=calc_evap(tf,pf)
          tw=calc_wetbulb(e,pf,df,tf)
    ; NOW REDO FOR ICE BULBS - BIT OF ERROR AS Twet CALC FROM Ewater not Eice
          icicles=WHERE(tw LE 0.,count)
          IF (count GT 0) THEN BEGIN
            e(icicles)=calc_evap_ice(df(icicles),pf(icicles))
            esat(icicles)=calc_evap_ice(tf(icicles),pf(icicles))
            tw(icicles)=calc_wetbulb(e(icicles),pf(icicles),df(icicles),tf(icicles))
          ENDIF
          IF (param EQ 'q2m') THEN q=calc_qhum(e,pf)
          IF (param EQ 'rh2m') THEN rh=(e/esat)*100.
        ENDIF
	IF (param EQ 'dpd2m') THEN dpd=tf-df
	
        pointst=0
        pointed=-1
        FOR pp=0,72 DO BEGIN
	  pcount=((yrr-styr)*73)+pp
          pointst=pointed+1
          pointed=pointarr(pp)-1
          IF (param EQ 'e2m') THEN pentad(lnn,ltt,pcount)=mean(e(pointst:pointed))
          IF (param EQ 'q2m') THEN pentad(lnn,ltt,pcount)=mean(q(pointst:pointed))
          IF (param EQ 'rh2m') THEN pentad(lnn,ltt,pcount)=mean(rh(pointst:pointed))
          IF (param EQ 'td2m') THEN pentad(lnn,ltt,pcount)=mean(df(pointst:pointed))
          IF (param EQ 't2m') THEN pentad(lnn,ltt,pcount)=mean(tf(pointst:pointed))
          IF (param EQ 'tw2m') THEN pentad(lnn,ltt,pcount)=mean(tw(pointst:pointed))
          IF (param EQ 'dpd2m') THEN pentad(lnn,ltt,pcount)=mean(dpd(pointst:pointed))
          IF (param EQ 'p2m') THEN pentad(lnn,ltt,pcount)=mean(pf(pointst:pointed))
;          print,pp,pointst,pointed
        ENDFOR
      ENDFOR
;      print,ltt,lnn
    ENDFOR

    county=county+1
  ENDFOR
ENDFOR  

; IF its output for marine data then regrid and shift
IF (MARINE EQ 1) THEN BEGIN
  oldpentad=pentad
  pentad=make_array(nlons,marlats,nptds,type=FLOAT,value=mdi)
  FOR lnn=0,nlons-1 DO BEGIN
    FOR pt=0,nptds-1 DO BEGIN
      subarr=transpose(oldpentad(lnn,*,pt)) ; this is 90 to -90, so 181 boxes, needs to be transpose to make it all one row
      ; expand to repeat for half deg
      subarr=REFORM(transpose(REFORM([subarr,subarr],nlats,2)),nlats*2)
      ; get rid of redundant 90 to 90.5 and -90 to -90.5 boxes
      subarr=subarr(1:(nlats*2)-2)
      ; reform to two columns and lots of rows, and then compute means over each pair of columns
      pentad(lnn,*,pt)=mean(REFORM(subarr,2,marlats),dimension=1)    
    ENDFOR
  ENDFOR
  lons=findgen(nlons)-179.5
  lats=(findgen(marlats)*(-1))+89.5
  ; shift lons
  pentad=shift(pentad,180) ; should shift each lon,lat field 180 lons over and not touch the lats of pentad order
ENDIF

; make pentad climatologies and standard deviations and anomalies
IF (MARINE EQ 1) THEN nlats=marlats ; now 180 rather than 181
clims=make_array(nlons,nlats,73,type=FLOAT,value=mdi)
stdevs=make_array(nlons,nlats,73,type=FLOAT,value=mdi)
anoms=make_array(nlons,nlats,nptds,type=FLOAT,value=mdi)

FOR ltt=0,nlats-1 DO BEGIN
  FOR lnn=0,nlons-1 DO BEGIN
    subarr=REFORM(pentad(lnn,ltt,*),73,nyrs)
    FOR pp=0,72 DO BEGIN
      clims(lnn,ltt,pp) = mean(subarr(pp,clmst-styr:clmed-styr))
      stdevs(lnn,ltt,pp) = stdev(subarr(pp,clmst-styr:clmed-styr))
      subarr(pp,*)=subarr(pp,*)-clims(lnn,ltt,pp)    
    ENDFOR
    anoms(lnn,ltt,*)=REFORM(subarr,nptds)
;    stop
  ENDFOR
ENDFOR

;************************************************
; Save as netCDF

filename=outdir+outfile
file_out=NCDF_CREATE(filename,/clobber)
time_id=NCDF_DIMDEF(file_out,'time',nptds)
timvarid=NCDF_VARDEF(file_out,'time',[time_id],/DOUBLE)
pt_id=NCDF_DIMDEF(file_out,'pentad_time',73)
ptvarid=NCDF_VARDEF(file_out,'pentad_time',[pt_id],/DOUBLE)
lat_id=NCDF_DIMDEF(file_out,'latitude',nlats)
latvarid=NCDF_VARDEF(file_out,'latitude',[lat_id],/DOUBLE)
lon_id=NCDF_DIMDEF(file_out,'longitude',nlons)
lonvarid=NCDF_VARDEF(file_out,'longitude',[lon_id],/DOUBLE)

IF (param EQ 'e2m') THEN BEGIN
  varid=NCDF_VARDEF(file_out,'e2m',[lon_id,lat_id,time_id],/FLOAT)
  climid=NCDF_VARDEF(file_out,'e2m_clims',[lon_id,lat_id,pt_id],/FLOAT)
  sdid=NCDF_VARDEF(file_out,'e2m_stdevs',[lon_id,lat_id,pt_id],/FLOAT)
  anomid=NCDF_VARDEF(file_out,'e2m_anoms',[lon_id,lat_id,time_id],/FLOAT)

  NCDF_ATTPUT,file_out,timvarid,'units','days since 1979-1-1 00:00:00'
  NCDF_ATTPUT,file_out,timvarid,'long_name','Time'
  NCDF_ATTPUT,file_out,ptvarid,'units','pentads'
  NCDF_ATTPUT,file_out,ptvarid,'long_name','pentads of the year'
  NCDF_ATTPUT,file_out,varid,'units','hPa'
  NCDF_ATTPUT,file_out,varid,'long_name','2 metre vapour pressure from 6hrly T and Td'
  NCDF_ATTPUT,file_out,varid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,climid,'units','hPa'
  NCDF_ATTPUT,file_out,climid,'long_name','2 metre vapour pressure from 6hrly T and Td climatology 1981-2010'
  NCDF_ATTPUT,file_out,climid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,sdid,'units','hPa'
  NCDF_ATTPUT,file_out,sdid,'long_name','2 metre vapour pressure from 6hrly T and Td pentad standard deviation 1981-2010'
  NCDF_ATTPUT,file_out,sdid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,anomid,'units','hPa'
  NCDF_ATTPUT,file_out,anomid,'long_name','2 metre vapour pressure from 6hrly T and Td climate anomalies'
  NCDF_ATTPUT,file_out,anomid,'missing_value',-1e30
ENDIF 
IF (param EQ 'p2m') THEN BEGIN
  varid=NCDF_VARDEF(file_out,'p2m',[lon_id,lat_id,time_id],/FLOAT)
  climid=NCDF_VARDEF(file_out,'p2m_clims',[lon_id,lat_id,pt_id],/FLOAT)
  sdid=NCDF_VARDEF(file_out,'p2m_stdevs',[lon_id,lat_id,pt_id],/FLOAT)
  anomid=NCDF_VARDEF(file_out,'p2m_anoms',[lon_id,lat_id,time_id],/FLOAT)

  NCDF_ATTPUT,file_out,timvarid,'units','days since 1979-1-1 00:00:00'
  NCDF_ATTPUT,file_out,timvarid,'long_name','Time'
  NCDF_ATTPUT,file_out,ptvarid,'units','pentads'
  NCDF_ATTPUT,file_out,ptvarid,'long_name','pentads of the year'
  NCDF_ATTPUT,file_out,varid,'units','hPa'
  NCDF_ATTPUT,file_out,varid,'long_name','2 metre surface pressure from 6hrly'
  NCDF_ATTPUT,file_out,varid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,climid,'units','hPa'
  NCDF_ATTPUT,file_out,climid,'long_name','2 metre surface pressure from 6hrly climatology 1981-2010'
  NCDF_ATTPUT,file_out,climid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,sdid,'units','hPa'
  NCDF_ATTPUT,file_out,sdid,'long_name','2 metre surface pressure from 6hrly pentad standard deviation 1981-2010'
  NCDF_ATTPUT,file_out,sdid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,anomid,'units','hPa'
  NCDF_ATTPUT,file_out,anomid,'long_name','2 metre surface pressure from 6hrly climate anomalies'
  NCDF_ATTPUT,file_out,anomid,'missing_value',-1e30
ENDIF 
IF (param EQ 'tw2m') THEN BEGIN
  varid=NCDF_VARDEF(file_out,'tw2m',[lon_id,lat_id,time_id],/FLOAT)
  climid=NCDF_VARDEF(file_out,'tw2m_clims',[lon_id,lat_id,pt_id],/FLOAT)
  sdid=NCDF_VARDEF(file_out,'tw2m_stdevs',[lon_id,lat_id,pt_id],/FLOAT)
  anomid=NCDF_VARDEF(file_out,'tw2m_anoms',[lon_id,lat_id,time_id],/FLOAT)

  NCDF_ATTPUT,file_out,timvarid,'units','days since 1979-1-1 00:00:00'
  NCDF_ATTPUT,file_out,timvarid,'long_name','Time'
  NCDF_ATTPUT,file_out,ptvarid,'units','pentads'
  NCDF_ATTPUT,file_out,ptvarid,'long_name','pentads of the year'
  NCDF_ATTPUT,file_out,varid,'units','degrees C'
  NCDF_ATTPUT,file_out,varid,'long_name','2 metre wet bulb temperature from 6hrly T and Td'
  NCDF_ATTPUT,file_out,varid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,climid,'units','degrees C'
  NCDF_ATTPUT,file_out,climid,'long_name','2 metre wet bulb temperature from 6hrly T and Td climatology 1981-2010'
  NCDF_ATTPUT,file_out,climid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,sdid,'units','degrees C'
  NCDF_ATTPUT,file_out,sdid,'long_name','2 metre wet bulb temperature from 6hrly T and Td pentad standard deviation 1981-2010'
  NCDF_ATTPUT,file_out,sdid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,anomid,'units','degrees C'
  NCDF_ATTPUT,file_out,anomid,'long_name','2 metre wet bulb temperature from 6hrly T and Td climate anomalies'
  NCDF_ATTPUT,file_out,anomid,'missing_value',-1e30
ENDIF
IF (param EQ 'q2m') THEN BEGIN
  varid=NCDF_VARDEF(file_out,'q2m',[lon_id,lat_id,time_id],/FLOAT)
  climid=NCDF_VARDEF(file_out,'q2m_clims',[lon_id,lat_id,pt_id],/FLOAT)
  sdid=NCDF_VARDEF(file_out,'q2m_stdevs',[lon_id,lat_id,pt_id],/FLOAT)
  anomid=NCDF_VARDEF(file_out,'q2m_anoms',[lon_id,lat_id,time_id],/FLOAT)

  NCDF_ATTPUT,file_out,timvarid,'units','days since 1979-1-1 00:00:00'
  NCDF_ATTPUT,file_out,timvarid,'long_name','Time'
  NCDF_ATTPUT,file_out,ptvarid,'units','pentads'
  NCDF_ATTPUT,file_out,ptvarid,'long_name','pentads of the year'
  NCDF_ATTPUT,file_out,varid,'units','g/kg'
  NCDF_ATTPUT,file_out,varid,'long_name','2 metre specific humidity from 6hrly T and Td'
  NCDF_ATTPUT,file_out,varid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,climid,'units','g/kg'
  NCDF_ATTPUT,file_out,climid,'long_name','2 metre specific humidity from 6hrly T and Td climatology 1981-2010'
  NCDF_ATTPUT,file_out,climid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,sdid,'units','g/kg'
  NCDF_ATTPUT,file_out,sdid,'long_name','2 metre specific humidity from 6hrly T and Td pentad standard deviation 1981-2010'
  NCDF_ATTPUT,file_out,sdid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,anomid,'units','g/kg'
  NCDF_ATTPUT,file_out,anomid,'long_name','2 metre specific humidity from 6hrly T and Td climate anomalies'
  NCDF_ATTPUT,file_out,anomid,'missing_value',-1e30
ENDIF
IF (param EQ 'rh2m') THEN BEGIN
  varid=NCDF_VARDEF(file_out,'rh2m',[lon_id,lat_id,time_id],/FLOAT)
  climid=NCDF_VARDEF(file_out,'rh2m_clims',[lon_id,lat_id,pt_id],/FLOAT)
  sdid=NCDF_VARDEF(file_out,'rh2m_stdevs',[lon_id,lat_id,pt_id],/FLOAT)
  anomid=NCDF_VARDEF(file_out,'rh2m_anoms',[lon_id,lat_id,time_id],/FLOAT)

  NCDF_ATTPUT,file_out,timvarid,'units','days since 1979-1-1 00:00:00'
  NCDF_ATTPUT,file_out,timvarid,'long_name','Time'
  NCDF_ATTPUT,file_out,ptvarid,'units','pentads'
  NCDF_ATTPUT,file_out,ptvarid,'long_name','pentads of the year'
  NCDF_ATTPUT,file_out,varid,'units','%rh'
  NCDF_ATTPUT,file_out,varid,'long_name','2 metre relative humidity from 6hrly T and Td'
  NCDF_ATTPUT,file_out,varid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,climid,'units','%rh'
  NCDF_ATTPUT,file_out,climid,'long_name','2 metre relative humidity from 6hrly T and Td climatology 1981-2010'
  NCDF_ATTPUT,file_out,climid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,sdid,'units','%rh'
  NCDF_ATTPUT,file_out,sdid,'long_name','2 metre relative humidity from 6hrly T and Td pentad standard deviation 1981-2010'
  NCDF_ATTPUT,file_out,sdid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,anomid,'units','%rh'
  NCDF_ATTPUT,file_out,anomid,'long_name','2 metre relative humidity from 6hrly T and Td climate anomalies'
  NCDF_ATTPUT,file_out,anomid,'missing_value',-1e30
ENDIF
IF (param EQ 't2m') THEN BEGIN
  varid=NCDF_VARDEF(file_out,'t2m',[lon_id,lat_id,time_id],/FLOAT)
  climid=NCDF_VARDEF(file_out,'t2m_clims',[lon_id,lat_id,pt_id],/FLOAT)
  sdid=NCDF_VARDEF(file_out,'t2m_stdevs',[lon_id,lat_id,pt_id],/FLOAT)
  anomid=NCDF_VARDEF(file_out,'t2m_anoms',[lon_id,lat_id,time_id],/FLOAT)

  NCDF_ATTPUT,file_out,timvarid,'units','days since 1979-1-1 00:00:00'
  NCDF_ATTPUT,file_out,timvarid,'long_name','Time'
  NCDF_ATTPUT,file_out,ptvarid,'units','pentads'
  NCDF_ATTPUT,file_out,ptvarid,'long_name','pentads of the year'
  NCDF_ATTPUT,file_out,varid,'units','degrees C'
  NCDF_ATTPUT,file_out,varid,'long_name','2 metre dry bulb temperature from 6hrly T and Td'
  NCDF_ATTPUT,file_out,varid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,climid,'units','degrees C'
  NCDF_ATTPUT,file_out,climid,'long_name','2 metre dry bulb temperature from 6hrly T and Td climatology 1981-2010'
  NCDF_ATTPUT,file_out,climid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,sdid,'units','degrees C'
  NCDF_ATTPUT,file_out,sdid,'long_name','2 metre dry bulb temperature from 6hrly T and Td pentad standard deviation 1981-2010'
  NCDF_ATTPUT,file_out,sdid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,anomid,'units','degrees C'
  NCDF_ATTPUT,file_out,anomid,'long_name','2 metre dry bulb temperature from 6hrly T and Td climate anomalies'
  NCDF_ATTPUT,file_out,anomid,'missing_value',-1e30
ENDIF
IF (param EQ 'td2m') THEN BEGIN
  varid=NCDF_VARDEF(file_out,'td2m',[lon_id,lat_id,time_id],/FLOAT)
  climid=NCDF_VARDEF(file_out,'td2m_clims',[lon_id,lat_id,pt_id],/FLOAT)
  sdid=NCDF_VARDEF(file_out,'td2m_stdevs',[lon_id,lat_id,pt_id],/FLOAT)
  anomid=NCDF_VARDEF(file_out,'td2m_anoms',[lon_id,lat_id,time_id],/FLOAT)

  NCDF_ATTPUT,file_out,timvarid,'units','days since 1979-1-1 00:00:00'
  NCDF_ATTPUT,file_out,timvarid,'long_name','Time'
  NCDF_ATTPUT,file_out,ptvarid,'units','pentads'
  NCDF_ATTPUT,file_out,ptvarid,'long_name','pentads of the year'
  NCDF_ATTPUT,file_out,varid,'units','degrees C'
  NCDF_ATTPUT,file_out,varid,'long_name','2 metre dew point temperature from 6hrly T and Td'
  NCDF_ATTPUT,file_out,varid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,climid,'units','degrees C'
  NCDF_ATTPUT,file_out,climid,'long_name','2 metre dew point temperature from 6hrly T and Td climatology 1981-2010'
  NCDF_ATTPUT,file_out,climid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,sdid,'units','degrees C'
  NCDF_ATTPUT,file_out,sdid,'long_name','2 metre dew point temperature from 6hrly T and Td pentad standard deviation 1981-2010'
  NCDF_ATTPUT,file_out,sdid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,anomid,'units','degrees C'
  NCDF_ATTPUT,file_out,anomid,'long_name','2 metre dew point temperature from 6hrly T and Td climate anomalies'
  NCDF_ATTPUT,file_out,anomid,'missing_value',-1e30
ENDIF
IF (param EQ 'dpd2m') THEN BEGIN
  varid=NCDF_VARDEF(file_out,'dpd2m',[lon_id,lat_id,time_id],/FLOAT)
  climid=NCDF_VARDEF(file_out,'dpd2m_clims',[lon_id,lat_id,pt_id],/FLOAT)
  sdid=NCDF_VARDEF(file_out,'dpd2m_stdevs',[lon_id,lat_id,pt_id],/FLOAT)
  anomid=NCDF_VARDEF(file_out,'dpd2m_anoms',[lon_id,lat_id,time_id],/FLOAT)

  NCDF_ATTPUT,file_out,timvarid,'units','days since 1979-1-1 00:00:00'
  NCDF_ATTPUT,file_out,timvarid,'long_name','Time'
  NCDF_ATTPUT,file_out,ptvarid,'units','pentads'
  NCDF_ATTPUT,file_out,ptvarid,'long_name','pentads of the year'
  NCDF_ATTPUT,file_out,varid,'units','degrees C'
  NCDF_ATTPUT,file_out,varid,'long_name','2 metre dew point depression from 6hrly T and Td'
  NCDF_ATTPUT,file_out,varid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,climid,'units','degrees C'
  NCDF_ATTPUT,file_out,climid,'long_name','2 metre dew point depression from 6hrly T and Td climatology 1981-2010'
  NCDF_ATTPUT,file_out,climid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,sdid,'units','degrees C'
  NCDF_ATTPUT,file_out,sdid,'long_name','2 metre dew point depression from 6hrly T and Td pentad standard deviation 1981-2010'
  NCDF_ATTPUT,file_out,sdid,'missing_value',-1e30
  NCDF_ATTPUT,file_out,anomid,'units','degrees C'
  NCDF_ATTPUT,file_out,anomid,'long_name','2 metre dew point depression from 6hrly T and Td climate anomalies'
  NCDF_ATTPUT,file_out,anomid,'missing_value',-1e30
ENDIF

NCDF_ATTPUT,file_out,latvarid,'units','degrees_north'
NCDF_ATTPUT,file_out,latvarid,'long_name','latitude'
NCDF_ATTPUT,file_out,lonvarid,'units','degrees_east'
NCDF_ATTPUT,file_out,lonvarid,'long_name','longitude'

NCDF_CONTROL,file_out,/ENDEF

NCDF_VARPUT,file_out,varid,pentad
NCDF_VARPUT,file_out,climid,clims
NCDF_VARPUT,file_out,sdid,stdevs
NCDF_VARPUT,file_out,anomid,anoms
NCDF_VARPUT,file_out,latvarid,lats
NCDF_VARPUT,file_out,lonvarid,lons
NCDF_VARPUT,file_out,time_id,goomoo
NCDF_CLOSE,file_out

;************************************************

;stop
end
