pro makeOTHER_area_avg_ts

; to calculate area average - NO missing data tolerance - need a land/sea mask for that

; read in files
; make sure they are SW to NE
; convert missing data
; mask for land/stations present gridboxes only
; extract appropriate time chunk
; create anomalies 1981-2010
; regrid to HadISDH grids - apply cosine(lat) weigting NOT DONE YET
; create a HadISDH masked version too
; Calculate global, hemisperic and tropical average series
; save to netCDF

; written by Kate Willett
; last modified JUL 2012

;---------------------------------------------------
; set up directories and filenames
mdi=-1e+30

param='GHCNM'	;'CRUTEM4','CRUTS3_21','GHCNM','ERAINTERIM','GISS','BERKELEY'
param2='T'	;'DPD','Td','T','e','q','RH'
nowmon='MAY'
nowyear='2015'

styr=1973	; 1973
edyr=2014	; 2013
nyrs=(edyr+1)-styr
nmons=nyrs*12
int_mons=indgen(nmons)
climst=1976
climed=2005

indir='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/'
indirH='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/GRIDS/'
outdir='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TIMESERIES/'

CASE param OF
  'BERKELEY': BEGIN
    infile='Complete_TAVG_LatLong1'	;temperature_anomaly
    latlg=1.					;latitude
    lonlg=1.					;longitude
    latst=-89.5
    lonst=-179.5
    stlt=-89.5
    stln=-179.5
    actst=1750
    acted=2014 ; even though this is April 2015 it works if you say 2014 here as it cuts to nmons
    actmoned=12 ; even though this is April 2015 it works if you say 2014 here as it cuts to nmons
    actmdi='NaN'
    scaled=0
    anom=1
    type='nc'
  END
  'GISS': BEGIN
    infile='gistemp250'	;temperature_anomaly
    latlg=2.					;latitude
    lonlg=2.					;longitude
    latst=-89.
    lonst=-179.
    stlt=-89.
    stln=-179.
    actst=1880
    acted=2015
    actmoned=4
    actmdi=327.67 ; once scaled by 0.01
    scaled=1
    anom=1
    type='nc'
  END
  'CRUTEM4': BEGIN
    infile='CRUTEM.4.3.0.0.anomalies'	;temperature_anomaly
    latlg=5.					;latitude
    lonlg=5.					;longitude
    latst=-87.5
    lonst=-177.5
    stlt=-87.5
    stln=-177.5
    actst=1850
    acted=2014
    actmoned=12
    actmdi=-1.e+30
    scaled=0
    anom=1
    type='nc'
  END
  'CRUTS3_21': BEGIN
    infile='cru_ts3.21.1901.2012.vap.dat'	;vap
    innobs='cru_ts3.21.1901.2012.vap.stn'	;stn
    latlg=0.5					;lat
    lonlg=0.5					;lon
    latst=-89.75
    lonst=-179.75
    stlt=-89.75
    stln=-179.75
    actst=1901
    acted=2012
    actmoned=12
    actmdi=9.969e+36
    scaled=0
    anom=0
    type='nc'
  END
  'ERAINTERIM': BEGIN
    inlandmask='new_coverpercentjul08'
    IF (param2 EQ 'T') OR (param2 EQ 'Td') OR (param2 EQ 'DPD') THEN infile='ERA-INTERIM_t_td_19792013mm'	;t2m
    IF (param2 EQ 'Tw') OR (param2 EQ 'e') THEN BEGIN
;      infile='ERA-INTERIM_t_td_19792013mm'	;t2m
;      infile2='ERA-INTERIM_p_19792013mm'
    ENDIF
    IF (param2 EQ 'e') THEN infile='e2m_monthly_ERA-Interim_data'
    IF (param2 EQ 'Tw') THEN infile='Tw2m_monthly_ERA-Interim_data'
    IF (param2 EQ 'q') THEN infile='q2m_monthly_ERA-Interim_data'
    IF (param2 EQ 'RH') THEN infile='rh2m_monthly_ERA-Interim_data'
    IF (param2 EQ 'T') OR (param2 EQ 'Td') OR (param2 EQ 'DPD') THEN BEGIN ;$
;       OR (param2 EQ 'Tw') OR (param2 EQ 'e') THEN BEGIN
      latlg=0.75 					;latitude, 1 extra box
      lonlg=0.75 			      ;longitude
      latst=90 
      lonst=0 
      stlt=-90
      stln=-180
      actst=1979
      acted=2014
      actmoned=12
      actmdi=-32767
      scaled=1
      anom=0
      type='nc' 
    ENDIF ELSE IF (param2 EQ 'Tw') OR (param2 EQ 'e') OR (param2 EQ 'q')OR (param2 EQ 'RH') THEN BEGIN
      latlg=1.
      lonlg=1.				      ;longitude
      latst=90
      lonst=0
      stlt=-90
      stln=-180
      actst=1979
      acted=2014
      actmoned=12
      actmdi=-32767
      scaled=0
      anom=0
      type='nc'
    ENDIF ELSE BEGIN
      latlg=1.
      lonlg=1.				      ;longitude
      latst=89.5
      lonst=-179.5
      stlt=-89.5
      stln=-179.5
      actst=1979
      acted=2014
      actmoned=12
      actmdi=-32767
      scaled=0
      anom=0
      type='txt'
    ENDELSE
    
  END
  'GHCNM': BEGIN
    infile='grid-mntp-1880-current-v3.2.2.dat'
    latlg=5.					
    lonlg=5.					
    latst=87.5
    lonst=-177.5
    stlt=-87.5
    stln=-177.5
    actst=1880
    acted=2015
    actmoned=12
    actmdi=-99.99	;-99.99
    scaled=1		; /100.
    anom=1
    type='txt'
  END
ENDCASE
;ofileFULL=param+'_'+param2+'_areaTSkate_'+strcompress(styr,/remove_all)+strcompress(edyr,/remove_all)	
;ofileMASK=param+'_'+param2+'_HadISDHMASKareaTSkate_'+strcompress(styr,/remove_all)+strcompress(edyr,/remove_all)	
;ofileGRIDS=param+'_'+param2+'_5by519812010clim7913anomskate_'+strcompress(styr,/remove_all)+strcompress(edyr,/remove_all)	
ofileFULL=param+'_'+param2+'_areaTS_'+strcompress(styr,/remove_all)+strcompress(edyr,/remove_all)	
ofileMASK=param+'_'+param2+'_HadISDHMASKareaTS_'+strcompress(styr,/remove_all)+strcompress(edyr,/remove_all)	
ofileGRIDS=param+'_'+param2+'_5by519762005clim_anoms_'+strcompress(styr,/remove_all)+strcompress(edyr,/remove_all)	

CASE param2 OF
  'RH': BEGIN
    unitees='% rh'
    inhadisdh='HadISDH.landRH.2.0.0.2013p_FLATgridIDPHA5by5_JAN2014.nc'
  END
  'e': BEGIN
    unitees='hPa'
    inhadisdh='HadISDH.lande.2.0.0.2013p_FLATgridIDPHA5by5_JAN2014.nc'
  END
  'q': BEGIN
    unitees='g/kg'
    inhadisdh='HadISDH.landq.2.0.0.2013p_FLATgridIDPHA5by5_JAN2014.nc'
  END
  'T': BEGIN
    unitees='deg C'
    inhadisdh='HadISDH.landT.2.0.1.2014p_FLATgridIDPHA5by5_JAN2015_cf.nc'
  END
  'Td': BEGIN
    unitees='deg C'
    inhadisdh='HadISDH.landTd.2.0.0.2013p_FLATgridPHADPD5by5_JAN2014.nc'
  END
  'Tw': BEGIN
    unitees='deg C'
    inhadisdh='HadISDH.landTw.2.0.0.2013p_FLATgridIDPHA5by5_JAN2014.nc'
  END
  'DPD': BEGIN
    unitees='deg C'
    inhadisdh='HadISDH.landDPD.2.0.0.2013p_FLATgridPHA5by5_JAN2014.nc'
  END
ENDCASE

nlats=180/latlg
IF (param EQ 'ERAINTERIM') AND ((param2 EQ 'T') OR (param2 EQ 'Td') $
   OR (param2 EQ 'DPD') OR (param2 EQ 'Tw') OR (param2 EQ 'e') OR (param2 EQ 'q') OR (param2 EQ 'RH')) THEN nlats=nlats+1
nlons=360/lonlg
nbox=LONG(nlats*nlons)

lats=(findgen(nlats)*latlg)+stlt
lons=(findgen(nlons)*lonlg)+stln

valsarr=make_array(nlons,nlats,nmons,/float,value=mdi)
anomalies=make_array(nlons,nlats,nmons,/float,value=mdi)

hstlt=-87.5
hstln=-177.5
hlatlg=5.
hlonlg=5.
hnlats=180/hlatlg
hnlons=360/hlonlg
hlats=(findgen(hnlats)*hlatlg)+hstlt
hlons=(findgen(hnlons)*hlonlg)+hstln

; making hi res lats at 1/8th of a degree res.
hrstlt=-89.9375
hrlatlg=0.125
hrnlats=180/hrlatlg
hrlats=(findgen(hrnlats)*hrlatlg)+hrstlt
hireslatsarr=fltarr(360/hrlatlg,hrnlats)
FOR ltt=0,hrnlats-1 DO BEGIN
  hireslatsarr(*,ltt)=hrlats(ltt)
ENDFOR
;stop

newanomalies=make_array(hnlons,hnlats,nmons,/float,value=mdi)

print,param,param2,actst,acted

;----------------------------------------------------
; read in files

IF (type EQ 'nc') THEN BEGIN
  filee=NCDF_OPEN(indir+infile+'.nc')
  CASE param OF
    'BERKELEY': varid=NCDF_VARID(filee,'temperature')
    'GISS': varid=NCDF_VARID(filee,'tempanomaly')
    'CRUTEM4': varid=NCDF_VARID(filee,'temperature_anomaly')
    'CRUTS3_21': varid=NCDF_VARID(filee,'vap')
    'ERAINTERIM': BEGIN
      IF (param2 EQ 'T') THEN varid=NCDF_VARID(filee,'t2m')
      IF (param2 EQ 'Td') THEN varid=NCDF_VARID(filee,'d2m')
      IF (param2 EQ 'q') THEN varid=NCDF_VARID(filee,'q2m')
      IF (param2 EQ 'RH') THEN varid=NCDF_VARID(filee,'rh2m')
      IF (param2 EQ 'DPD') THEN BEGIN ;OR (param2 EQ 'Tw') OR (param2 EQ 'e') THEN BEGIN
        varid=NCDF_VARID(filee,'t2m')
        varid2=NCDF_VARID(filee,'d2m')
      ENDIF
      IF (param2 EQ 'Tw') THEN varid=NCDF_VARID(filee,'tw2m')
      IF (param2 EQ 'e') THEN varid=NCDF_VARID(filee,'e2m')      
    END
  ENDCASE
  NCDF_VARGET,filee,varid,var
  IF (scaled EQ 1) THEN BEGIN
    NCDF_ATTGET,filee,varid,'scale_factor',scalee
    IF (param EQ 'GISS') THEN BEGIN
      var=var*scalee
;      stop
    ENDIF ELSE BEGIN
      NCDF_ATTGET,filee,varid,'add_offset',offsetee
;    stop,'CHECK scale and offset'
      var=((var*scalee)+offsetee)-273.15
    ENDELSE
  ENDIF
  bads=WHERE(var LT -150,count)
  IF (count GT 0) THEN var(bads)=mdi
  bads=WHERE(var GT 150,count)
  IF (count GT 0) THEN var(bads)=mdi
  IF (param EQ 'ERAINTERIM') AND (param2 EQ 'DPD') THEN BEGIN; ((param2 EQ 'DPD') THEN BEGIN $
;     OR (param2 EQ 'Tw') OR (param2 EQ 'e')) THEN BEGIN
    NCDF_VARGET,filee,varid2,var2
    IF (scaled EQ 1) THEN BEGIN
      NCDF_ATTGET,filee,varid2,'scale_factor',scalee2
      NCDF_ATTGET,filee,varid2,'add_offset',offsetee2
      var2=((var2*scalee2)+offsetee2)-273.15
    ENDIF
    bads=WHERE(var2 LT -100,count)
    IF (count GT 0) THEN var2(bads)=mdi
  ENDIF
  
  IF (param EQ 'GISS') THEN NCDF_ATTGET,filee,varid,'_FillValue',actmdi ELSE NCDF_ATTGET,filee,varid,'missing_value',actmdi
  
  NCDF_CLOSE,filee

;  IF (param EQ 'ERAINTERIM') AND ((param2 EQ 'Tw') OR (param2 EQ 'e')) THEN BEGIN
;    filee=NCDF_OPEN(indir+infile2+'.nc')
;    varid=NCDF_VARID(filee,'sp')
;    NCDF_VARGET,filee,varid,varP
;    NCDF_ATTGET,filee,varid,'scale_factor',scalee
;    NCDF_ATTGET,filee,varid,'add_offset',offsetee
;    stop
;    varP=((varP*scalee)+offsetee)/100.
;    stop
;    NCDF_CLOSE,filee
;  ENDIF
  
  IF (param EQ 'ERAINTERIM') AND (param2 EQ 'DPD') THEN BEGIN
    goods=WHERE(var > mdi,count)
    var(goods)=var(goods)-var2(goods)
  ENDIF ;ELSE IF (param EQ 'ERAINTERIM') AND (param2 EQ 'Tw') THEN BEGIN
;    ; Td, P to calc e
;    evar=calc_evap(var2,varP)
;    ; T, Td, P and e to calc Tw
;    twvar=calc_wetbulb(evar,varP,var2,var)
;    ; where Tw <= 0 then eice calculated
;    lows=WHERE(twvar LE 0.,count)
;    IF (count GT 0) THEN evar(lows)=calc_evap_ice(var2(lows),varP(lows))
;    ; Tw recalculated from T, Td, e/eice and P
;    var=calc_wetbulb(evar,varP,var2,var)
;  ENDIF ELSE IF (param EQ 'ERAINTERIM') AND (param2 EQ 'e') THEN BEGIN
;    ; Td, P to calc e
;    evar=calc_evap(var2,varP)
;    ; T, Td, P and e to calc Tw
;    twvar=calc_wetbulb(evar,varP,var2,var)
;    ; where Tw <= 0 then eice calculated
;    lows=WHERE(twvar LE 0.,count)
;    IF (count GT 0) THEN evar(lows)=calc_evap_ice(var2(lows),varP(lows))
;    var=evar
;  ENDIF

  IF (param EQ 'CRUTS3_21') THEN BEGIN
    filee=NCDF_OPEN(indir+innobs+'.nc')
    varid=NCDF_VARID(filee,'stn')
    NCDF_VARGET,filee,varid,nobs
    NCDF_CLOSE,filee
    bads=WHERE(nobs NE actmdi AND nobs LT 1,count)
    IF (count GT 0) THEN var(bads)=actmdi
  ENDIF
  
ENDIF ELSE IF (type EQ 'txt') THEN BEGIN
  anmons=((acted+1)-actst)*12
  var=make_array(nlons,nlats,anmons,/float,value=mdi)
  openr,1,indir+infile
  IF (param EQ 'GHCNM') THEN BEGIN
    counter=0
    yrcount=0
    moncount=0
    totcount=0
    WHILE NOT EOF(1) DO BEGIN
      dummy=' '
      tmparr=intarr(72)
      IF (counter EQ 0) THEN readf,1,dummy,format='(a)' ELSE readf,1,tmparr,format='(72i6)'
      IF (counter NE 0) THEN var(*,counter-1,totcount)=tmparr/100.
      counter=counter+1
      IF (counter EQ 37) THEN BEGIN
        moncount=moncount+1
        totcount=totcount+1
        IF (moncount EQ 12) THEN BEGIN
          moncount=0
	  yrcount=yrcount+1
        ENDIF
        counter=0
      ENDIF    
    ENDWHILE
  ENDIF ELSE IF (param EQ 'ERAINTERIM') THEN BEGIN
    dummy=' '
    rowcount=0
    latcount=0
    totcount=0
    stcl=0
    edcl=5
    WHILE NOT EOF(1) DO BEGIN
      dummy=' '
      tmparr=intarr(6)
      IF (rowcount NE 0) THEN BEGIN
        readf,1,tmparr,format='(6f8.4)'
        var(stcl:edcl,latcount,totcount)=tmparr
	stcl=stcl+6
	edcl=edcl+6
	IF (edcl EQ 365) THEN BEGIN
	  stcl=0
	  edcl=5
	  latcount=latcount+1
	  rowcount=-1	;add 1 during rest of loop to bring to 0
	  IF (latcount EQ 180) THEN BEGIN
	    latcount=0
	    totcount=totcount+1
	  ENDIF
	ENDIF
      ENDIF ELSE BEGIN
        readf,1,dummy,format='(a)'
        IF (latcount EQ 0) THEN BEGIN
          readf,1,dummy,format='(a)'
          readf,1,dummy,format='(a)'
	ENDIF
      ENDELSE
      rowcount=rowcount+1
    ENDWHILE
  ENDIF
  close,1
ENDIF

;stop

IF (actst LE styr) THEN othstmon=(styr-actst)*12 ELSE othstmon=0 
IF (actst GT styr) THEN mystmon=(actst-styr)*12 ELSE mystmon=0 
IF (acted LT edyr) THEN myedmon=(((acted)-styr)*12)-1+actmoned ELSE myedmon=nmons-1
IF (acted GT edyr) THEN othedmon=(((edyr+1)-actst)*12)-1 ELSE othedmon=othstmon+nmons-1    ;n_elements(var(0,0,*))-1

;stop

valsarr(*,*,mystmon:myedmon)=var(*,*,othstmon:othedmon)

IF (param EQ 'BERKELEY') THEN bads=WHERE(Finite(valsarr) EQ 0,count) ELSE bads=WHERE(valsarr EQ actmdi,count)
IF (count GT 0) THEN valsarr(bads)=mdi

; reverse arrays where lats begin in the North
IF (latst GT 0) THEN BEGIN
  FOR mm=0,nmons-1 DO BEGIN
    valsarr(*,*,mm)=REVERSE(valsarr(*,*,mm),2)
  ENDFOR
ENDIF

; shift lons not beginning at -180 
IF (lonst GT -175.) THEN BEGIN
  FOR mm=0,nmons-1 DO BEGIN
    valsarr(*,*,mm)=SHIFT(valsarr(*,*,mm),(nlons/2),0)
  ENDFOR
ENDIF

    
; make anomalies from the monthlies/re normalise to clim period
FOR ltt=0,nlats-1 DO BEGIN
  FOR lnn=0,nlons-1 DO BEGIN
    subarr=REFORM(valsarr(lnn,ltt,*),12,nyrs)
    FOR mm=0,11 DO BEGIN
      gots=WHERE(subarr(mm,*) NE mdi,count)
      climsub=subarr(mm,climst-styr:climed-styr)
      gotsC=WHERE(climsub NE mdi,countC)
      IF (countC GE 15) AND (count GE 15) THEN subarr(mm,gots)=subarr(mm,gots)-MEAN(climsub(gotsC)) ELSE subarr(mm,*)=mdi
    ENDFOR
    anomalies(lnn,ltt,*)=REFORM(subarr,nmons)
  ENDFOR
ENDFOR

; regrid to HadISDH grids
; ERA T/Td is 0.75 centred on 90 so 0.375 is half a gridbox
; ERA e/Tq is 1. centred on 90 so 0.5 is half a gridbox
IF (latlg NE 5.) THEN BEGIN
  FOR mm=0,nmons-1 DO BEGIN
    print,'Month',mm
    ; make hi res 0.125 by 0.125 grids
    nlonshi=360*8
    nlatshi=180*8
    hires=make_array(nlonshi,nlatshi,/float,value=mdi)
    movelat=latlg/0.125
    movelon=lonlg/0.125
    sltt=0
    IF (param EQ 'ERAINTERIM') AND ((param2 EQ 'T') OR (param2 EQ 'Td') $
       OR (param2 EQ 'DPD') OR (param2 EQ 'Tw') OR (param2 EQ 'e') OR (param2 EQ 'q') OR (param2 EQ 'RH')) THEN eltt=(movelat/2.)-1 ELSE eltt=movelat-1
    FOR ltt=0,nlats-1 DO BEGIN ; ERA extra box
      slnn=0
      IF (param EQ 'ERAINTERIM') AND ((param2 EQ 'T') OR (param2 EQ 'Td') $
         OR (param2 EQ 'DPD') OR (param2 EQ 'Tw') OR (param2 EQ 'e') OR (param2 EQ 'q') OR (param2 EQ 'RH')) THEN elnn=(movelon/2.)-1 ELSE elnn=movelon-1    
      FOR lnn=0,nlons-1 DO BEGIN
;        print,sltt,eltt,slnn,elnn
        hires(slnn:elnn,sltt:eltt)=anomalies(lnn,ltt,mm)
	slnn=elnn+1
	elnn=elnn+movelon
	IF (lnn EQ nlons-1) AND (param EQ 'ERAINTERIM') AND $
	   ((param2 EQ 'T') OR (param2 EQ 'Td') OR (param2 EQ 'DPD') $
	   OR (param2 EQ 'Tw') OR (param2 EQ 'e') OR (param2 EQ 'q') OR (param2 EQ 'RH')) THEN BEGIN
	  elnn=elnn-(movelon/2.)
          hires(slnn:elnn,sltt:eltt)=anomalies(0,ltt,mm)
	ENDIF
      ENDFOR
;      stop
      sltt=eltt+1
      eltt=eltt+movelat
      IF (ltt EQ nlats-2) AND (param EQ 'ERAINTERIM') AND $
          ((param2 EQ 'T') OR (param2 EQ 'Td') OR (param2 EQ 'DPD') $
	  OR (param2 EQ 'Tw') OR (param2 EQ 'e') OR (param2 EQ 'q') OR (param2 EQ 'RH')) $
	   THEN eltt=eltt-(movelat/2.)
    ENDFOR 
;    stop
    
    movelat=hlatlg/0.125
    movelon=hlonlg/0.125
    sltt=0
    eltt=movelat-1
    FOR ltt=0,hnlats-1 DO BEGIN
      slnn=0
      elnn=movelon-1        
      FOR lnn=0,hnlons-1 DO BEGIN
        subarr=hires(slnn:elnn,sltt:eltt)
	sublats=hireslatsarr(slnn:elnn,sltt:eltt)
	gots=WHERE(subarr NE mdi,count)
        IF (count GT 0) THEN BEGIN ; Don't have a minimum requirement because we say one station is sufficient to create a gridbox average
;        IF (count GT 800) THEN BEGIN ; total boxes is 40*40 so 1600 - made a minimum requirement to prevent too much extrapolation
	  totsarr=TOTAL(subarr(gots)*cos(sublats(gots)*!dtor))
	  totslats=TOTAL(cos(sublats(gots)*!dtor))
	  newanomalies(lnn,ltt,mm)=totsarr/totslats
;	  stop
	ENDIF
	slnn=elnn+1
	elnn=elnn+movelon
      ENDFOR
      sltt=eltt+1
      eltt=eltt+movelat
    ENDFOR
;    stop
  ENDFOR
ENDIF ELSE newanomalies=anomalies

anomalies=newanomalies ; reset anomalies to newanomalies, newanomalies gets masked to HadISDH
;---------------------------------------------------------------

;stop

;make mask - set anything greater than 70 deg lat to mdi
inn=NCDF_OPEN(indirH+inhadisdh)
CASE param2 OF
  'RH': varid=NCDF_VARID(inn,'rh_anoms')
  'q': varid=NCDF_VARID(inn,'q_anoms')
  'e': varid=NCDF_VARID(inn,'e_anoms')
  'T': varid=NCDF_VARID(inn,'t_anoms')
  'Td': varid=NCDF_VARID(inn,'td_anoms')
  'Tw': varid=NCDF_VARID(inn,'tw_anoms')
  'DPD': varid=NCDF_VARID(inn,'dpd_anoms')
ENDCASE
NCDF_VARGET,inn,varid,hadisdh
NCDF_CLOSE,inn

IF (param EQ 'ERAINTERIM') THEN BEGIN
  ; get land mask for ERA
  inn=NCDF_OPEN(indir+inlandmask+'.nc')
  varid=NCDF_VARID(inn,'pct_land')
  NCDF_VARGET,inn,varid,pctland
  NCDF_CLOSE,inn

  pctland=REVERSE(pctland,2)
  
  ; Mask to get rid of marine only data and then to HadISDH coverage
  anomalies=newanomalies
  FOR mm=0,nmons-1 DO BEGIN
    if (mm GE n_elements(hadisdh(0,0,*))) THEN subhad=hadisdh(*,*,mm-12) ELSE subhad=hadisdh(*,*,mm) ; copes with ERA when its longer than HadISDH
    suboth=newanomalies(*,*,mm)
    bads=WHERE(subhad LT -900 and pctland EQ 0.,count)
    print,'BADS: ',count
    IF (count GT 0) THEN suboth(bads)=mdi
    anomalies(*,*,mm)=suboth ; anomalies = unmasked, newanomalies = masked
    bads=WHERE(subhad LT -900,count)
    IF (count GT 0) THEN suboth(bads)=mdi
    newanomalies(*,*,mm)=suboth
  ENDFOR
ENDIF ELSE BEGIN
  bads=WHERE(hadisdh LT -900.,count)
  IF (count GT 0) THEN newanomalies(bads)=mdi
ENDELSE
 
global_mask = anomalies(*,*,0)
global_mask(*,*) = 1
nhem_mask=global_mask
shem_mask=global_mask
trop_mask=global_mask

stlt=-87.5
stln=-177.5
latlg=5.
lonlg=5.
nlats=180/latlg
nlons=360/lonlg
lats=(findgen(nlats)*latlg)+stlt
lons=(findgen(nlons)*lonlg)+stln


for deg=0,nlats-1 DO BEGIN
  if (lats(deg) LT -70.) OR (lats(deg) GT 70.) THEN BEGIN
    global_mask(*,deg) = mdi
  endif
  if (lats(deg) LT 20.) OR (lats(deg) GT 70.) THEN BEGIN
    nhem_mask(*,deg) = mdi
  endif
  if (lats(deg) LT -70.) OR (lats(deg) GT -20.) THEN BEGIN
    shem_mask(*,deg) = mdi
  endif
  if (lats(deg) LT -20.) OR (lats(deg) GT 20.) THEN BEGIN
    trop_mask(*,deg) = mdi
  endif  
endfor


field_values=make_array(nlons,nlats,nmons,/float)	;all values or a subsection
num_val=n_elements(field_values(0,0,*))

global_mask_3d=fltarr(nlons,nlats,num_val)
nhem_mask_3d=fltarr(nlons,nlats,num_val)
shem_mask_3d=fltarr(nlons,nlats,num_val)
trop_mask_3d=fltarr(nlons,nlats,num_val)

for add=0,num_val-1 do begin
  global_mask_3d(*,*,add) = global_mask(*,*)
  nhem_mask_3d(*,*,add)=nhem_mask(*,*)
  shem_mask_3d(*,*,add)=shem_mask(*,*)
  trop_mask_3d(*,*,add)=trop_mask(*,*)
endfor


field_values(*,*,*)=anomalies(*,*,0:nmons-1)
print,n_elements(field_values(0,0,*))
num_val=n_elements(field_values(0,0,*))
glob_avg_ts=fltarr(num_val)
glob_avg_ts=globalmean(field_values,lats,mdi,mask=global_mask_3d)

print,n_elements(glob_avg_ts)
print,max(glob_avg_ts)
print,min(glob_avg_ts)
print,glob_avg_ts(0:10)

nhem_avg_ts=fltarr(num_val)
nhem_avg_ts=globalmean(field_values,lats,mdi,mask=nhem_mask_3d)

print,n_elements(nhem_avg_ts)
print,max(nhem_avg_ts)
print,min(nhem_avg_ts)
print,nhem_avg_ts(0:10)

shem_avg_ts=fltarr(num_val)
shem_avg_ts=globalmean(field_values,lats,mdi,mask=shem_mask_3d)

print,n_elements(shem_avg_ts)
print,max(shem_avg_ts)
print,min(shem_avg_ts)
print,shem_avg_ts(0:10)

trop_avg_ts=fltarr(num_val)
trop_avg_ts=globalmean(field_values,lats,mdi,mask=trop_mask_3d)

print,n_elements(trop_avg_ts)
print,max(trop_avg_ts)
print,min(trop_avg_ts)
print,trop_avg_ts(0:10)
  
;save to file;
;
filename=outdir+ofileFULL+'.nc'
file_out=NCDF_CREATE(filename,/clobber)
time_id=NCDF_DIMDEF(file_out,'time',nmons)
timvarid=NCDF_VARDEF(file_out,'time',[time_id],/DOUBLE)

glob_varid=NCDF_VARDEF(file_out,'glob_anoms',[time_id],/FLOAT)
trop_varid=NCDF_VARDEF(file_out,'trop_anoms',[time_id],/FLOAT)
nhem_varid=NCDF_VARDEF(file_out,'nhem_anoms',[time_id],/FLOAT)
shem_varid=NCDF_VARDEF(file_out,'shem_anoms',[time_id],/FLOAT)

NCDF_ATTPUT,file_out,timvarid,'units','months since 1973'
NCDF_ATTPUT,file_out,timvarid,'long_name','Time'
NCDF_ATTPUT,file_out,glob_varid,'units',unitees
NCDF_ATTPUT,file_out,glob_varid,'long_name','Globally Averaged Anomalies 70S-70N'
NCDF_ATTPUT,file_out,glob_varid,'missing_value',mdi

NCDF_ATTPUT,file_out,trop_varid,'units',unitees
NCDF_ATTPUT,file_out,trop_varid,'long_name','Tropics Averaged Anomalies 20S-20N'
NCDF_ATTPUT,file_out,trop_varid,'missing_value',mdi

NCDF_ATTPUT,file_out,nhem_varid,'units',unitees
NCDF_ATTPUT,file_out,nhem_varid,'long_name','Northern Hemisphere Averaged Anomalies 20N-70N'
NCDF_ATTPUT,file_out,nhem_varid,'missing_value',mdi

NCDF_ATTPUT,file_out,shem_varid,'units',unitees
NCDF_ATTPUT,file_out,shem_varid,'long_name','Southern Hemisphere Averaged Anomalies 70S-20S'
NCDF_ATTPUT,file_out,shem_varid,'missing_value',mdi
NCDF_CONTROL,file_out,/ENDEF

NCDF_VARPUT,file_out,glob_varid,glob_avg_ts
NCDF_VARPUT,file_out,trop_varid,trop_avg_ts
NCDF_VARPUT,file_out,nhem_varid,nhem_avg_ts
NCDF_VARPUT,file_out,shem_varid,shem_avg_ts
NCDF_VARPUT,file_out,time_id,int_mons
NCDF_CLOSE,file_out;

; mask with HadISDH

field_values(*,*,*)=newanomalies(*,*,0:nmons-1)
print,n_elements(field_values(0,0,*))
num_val=n_elements(field_values(0,0,*))
glob_avg_ts=fltarr(num_val)
glob_avg_ts=globalmean(field_values,lats,mdi,mask=global_mask_3d)

print,n_elements(glob_avg_ts)
print,max(glob_avg_ts)
print,min(glob_avg_ts)
print,glob_avg_ts(0:10)

nhem_avg_ts=fltarr(num_val)
nhem_avg_ts=globalmean(field_values,lats,mdi,mask=nhem_mask_3d)

print,n_elements(nhem_avg_ts)
print,max(nhem_avg_ts)
print,min(nhem_avg_ts)
print,nhem_avg_ts(0:10)

shem_avg_ts=fltarr(num_val)
shem_avg_ts=globalmean(field_values,lats,mdi,mask=shem_mask_3d)

print,n_elements(shem_avg_ts)
print,max(shem_avg_ts)
print,min(shem_avg_ts)
print,shem_avg_ts(0:10)

trop_avg_ts=fltarr(num_val)
trop_avg_ts=globalmean(field_values,lats,mdi,mask=trop_mask_3d)

print,n_elements(trop_avg_ts)
print,max(trop_avg_ts)
print,min(trop_avg_ts)
print,trop_avg_ts(0:10)
  
;save to file;
;
filename=outdir+ofileMASK+'.nc'
file_out=NCDF_CREATE(filename,/clobber)
time_id=NCDF_DIMDEF(file_out,'time',nmons)
timvarid=NCDF_VARDEF(file_out,'time',[time_id],/DOUBLE)

glob_varid=NCDF_VARDEF(file_out,'glob_anoms',[time_id],/FLOAT)
trop_varid=NCDF_VARDEF(file_out,'trop_anoms',[time_id],/FLOAT)
nhem_varid=NCDF_VARDEF(file_out,'nhem_anoms',[time_id],/FLOAT)
shem_varid=NCDF_VARDEF(file_out,'shem_anoms',[time_id],/FLOAT)

NCDF_ATTPUT,file_out,timvarid,'units','months since 1973'
NCDF_ATTPUT,file_out,timvarid,'long_name','Time'
NCDF_ATTPUT,file_out,glob_varid,'units',unitees
NCDF_ATTPUT,file_out,glob_varid,'long_name','Globally Averaged Anomalies 70S-70N'
NCDF_ATTPUT,file_out,glob_varid,'missing_value',mdi

NCDF_ATTPUT,file_out,trop_varid,'units',unitees
NCDF_ATTPUT,file_out,trop_varid,'long_name','Tropics Averaged Anomalies 20S-20N'
NCDF_ATTPUT,file_out,trop_varid,'missing_value',mdi

NCDF_ATTPUT,file_out,nhem_varid,'units',unitees
NCDF_ATTPUT,file_out,nhem_varid,'long_name','Northern Hemisphere Averaged Anomalies 20N-70N'
NCDF_ATTPUT,file_out,nhem_varid,'missing_value',mdi

NCDF_ATTPUT,file_out,shem_varid,'units',unitees
NCDF_ATTPUT,file_out,shem_varid,'long_name','Southern Hemisphere Averaged Anomalies 70S-20S'
NCDF_ATTPUT,file_out,shem_varid,'missing_value',mdi
NCDF_CONTROL,file_out,/ENDEF

NCDF_VARPUT,file_out,glob_varid,glob_avg_ts
NCDF_VARPUT,file_out,trop_varid,trop_avg_ts
NCDF_VARPUT,file_out,nhem_varid,nhem_avg_ts
NCDF_VARPUT,file_out,shem_varid,shem_avg_ts
NCDF_VARPUT,file_out,time_id,int_mons
NCDF_CLOSE,file_out;

IF (param EQ 'ERAINTERIM') OR (param EQ 'BERKELEY') OR (param EQ 'GISS') OR (param EQ 'GHCNM') THEN BEGIN
; output land masked arrays (5by5) to netCDF for later use

;save to file;
;
filename=indir+ofileGRIDS+'.nc'
file_out=NCDF_CREATE(filename,/clobber)
time_id=NCDF_DIMDEF(file_out,'time',nmons)
timvarid=NCDF_VARDEF(file_out,'time',[time_id],/DOUBLE)
lat_id=NCDF_DIMDEF(file_out,'lat',hnlats)
lon_id=NCDF_DIMDEF(file_out,'lon',hnlons)
lat_varid=NCDF_VARDEF(file_out,'lat',[lat_id],/DOUBLE)
lon_varid=NCDF_VARDEF(file_out,'lon',[lon_id],/DOUBLE)
varid=NCDF_VARDEF(file_out,'anomalies',[lon_id,lat_id,time_id],/DOUBLE)

NCDF_ATTPUT,file_out,timvarid,'units','months since 1973'
NCDF_ATTPUT,file_out,timvarid,'long_name','Time'
NCDF_ATTPUT,file_out,lat_varid,'units','degrees
NCDF_ATTPUT,file_out,lat_varid,'long_name','Gridbox centres in degrees south to north'
NCDF_ATTPUT,file_out,lon_varid,'units','degrees
NCDF_ATTPUT,file_out,lon_varid,'long_name','Gridbox centres in degrees west to east'

NCDF_ATTPUT,file_out,varid,'units',unitees
NCDF_ATTPUT,file_out,varid,'long_name','Climate anomalies from 1976-2005 land only'
;NCDF_ATTPUT,file_out,varid,'long_name','Climate anomalies from 1981-2010 land only'
NCDF_ATTPUT,file_out,varid,'missing_value',mdi
NCDF_CONTROL,file_out,/ENDEF

NCDF_VARPUT,file_out,varid,anomalies
NCDF_VARPUT,file_out,lat_varid,hlats
NCDF_VARPUT,file_out,lon_varid,hlons
NCDF_VARPUT,file_out,timvarid,int_mons
NCDF_CLOSE,file_out;

ENDIF


return

end

