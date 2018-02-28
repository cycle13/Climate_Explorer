; TIDL
; 
; Author: Kate Willett
; Created: 7 February 2017
; Last update: 7 February 2017
; Location: /data/local/hadkw/HADCRUH2/UPDATE2016/PROGS/IDL/	
; GitHub: https://github.com/Kate-Willett/Climate_Explorer					
; -----------------------
; CODE PURPOSE AND OUTPUT
; -----------------------
; This code reads in the HadISDH-marine data in its native format and
; outputs a HadISDH-land like format in addition to blending with HadISDH-land
; to product HadISDH-marine and HadISDH-blend
; 
; <references to related published material, e.g. that describes data set>
; 
; -----------------------
; LIST OF MODULES
; -----------------------
; <List of program modules required to run the code, or link to compiler/batch file>
; 
; -----------------------
; DATA
; -----------------------
; anomalies
; /project/hadobs2/hadisdh/marine/ICOADS.3.0.0/GRIDSOBSclim2BClocal/OBSclim2BClocal_5x5_monthly_renorm19812010_anomalies_from_daily_both_relax.nc
; actuals
; /project/hadobs2/hadisdh/marine/ICOADS.3.0.0/GRIDSOBSclim2BClocal/OBSclim2BClocal_5x5_monthly_from_daily_both_relax.nc
; climatology
; /project/hadobs2/hadisdh/marine/ICOADS.3.0.0/GRIDSOBSclim2BClocal/OBSclim2BClocal_5x5_monthly_climatology_from_daily_both_relax.nc
; standard deviation
; /project/hadobs2/hadisdh/marine/ICOADS.3.0.0/GRIDSOBSclim2BClocal/OBSclim2BClocal_5x5_monthly_stdev_from_daily_both_relax.nc
; or
; anomalies
; /project/hadobs2/hadisdh/marine/ICOADS.3.0.0/GRIDSOBSclim2BClocalship/OBSclim2BClocal_5x5_monthly_renorm19812010_anomalies_from_daily_both_relax.nc
; actuals
; /project/hadobs2/hadisdh/marine/ICOADS.3.0.0/GRIDSOBSclim2BClocalship/OBSclim2BClocal_5x5_monthly_from_daily_both_relax.nc
; climatology
; /project/hadobs2/hadisdh/marine/ICOADS.3.0.0/GRIDSOBSclim2BClocalship/OBSclim2BClocal_5x5_monthly_climatology_from_daily_both_relax.nc
; standard deviation
; /project/hadobs2/hadisdh/marine/ICOADS.3.0.0/GRIDSOBSclim2BClocalship/OBSclim2BClocal_5x5_monthly_stdev_from_daily_both_relax.nc
; may also be stored in /data/local/hadkw/HADCRUH2/UPDATE2016/OTHERDATA/
; The ship one I have changed its name to:
; OBSclim2BClocalship_5x5_monthly_renorm19812010_anomalies_from_daily_both_relax.nc
; AT A LATER DATE WE SHOULD ALSO READ IN ABS/CLIMS/STDEVS/UNCERTAINTIES to add to the package
; 
; -----------------------
; HOW TO RUN THE CODE
; -----------------------
; tidl
; .compile Blend_HadISDHLandOcean_APR2016.pro
; 
; -----------------------
; OUTPUT
; -----------------------
; /data/local/hadkw/HADCRUH2/UPDATE2016/STATISTICS/GRIDS/
; HadISDH-marineq.1.0.0.2016p_OBSclim2BClocal_anoms8110_JAN2017.nc
; HadISDH-blendq.1.0.0.2016p_FULL_anoms8110_JAN2017.nc (at the moment this involves renormalising the land data!)
; for the ship only version:
; HadISDH-marineq.1.0.0.2016p_OBSclim2BClocalship_anoms8110_JAN2017.nc
; HadISDH-blendq.1.0.0.2016p_FULLship_anoms8110_JAN2017.nc (at the moment this involves renormalising the land data!)
; 
; -----------------------
; VERSION/RELEASE NOTES
; -----------------------
;
; Version 1 (7 February 2017)
; ---------
;  
; Enhancements
;  
; Changes
;  
; Bug fixes
;  
; -----------------------
; OTHER INFORMATION
; -----------------------
;
;*******************************************************************

pro Blend_HadISDHLandOcean_APR2016

; Modified to read in from already created NOCS grids:
; extract_NOCSvars_FEB2014.pro
;	> UPDATE2014/OTHERDATA/NOCSv2.0_ocean*

;*******************************************************************
; EDITABLES:
mdi = -1e30

; What year?
styr = 1973 
edyr = 2016

; What climatology period? (land is generally 1976-2005 until I sort out a new version)
clstyr = 1981
cledyr = 2010
clst   = clstyr - styr
cled   = cledyr - styr
climchoice = strmid(strcompress(clstyr,/remove_all),2,2)+strmid(strcompress(cledyr,/remove_all),2,2)

; What versions?
LandVersion   = '3.0.0.2016p'
MarineVersion = '1.0.0.2016p'
BlendVersion  = '1.0.0.2016p'

; What Types?
LandType   = 'standard' ; this is actually 'FLATgrid*5by5' where IDPHA, PHADPD and PHA are used
;MarineType = 'OBSclim2BClocal' ; this is the QC'd and bias corrected version
;BlendType  = 'FULL' ; FULL=QC and bias corrected marine, homogenised land
MarineType = 'OBSclim2BClocalship' ; this is the QC'd and bias corrected version
BlendType  = 'FULLship' ; FULL=QC and bias corrected marine, homogenised land

; What working versions?
thenmon  = 'JAN'
thenyear = '2017'
nowmon   = 'JAN'
nowyear  = '2017'

; Work with anomalies or absolutes?
anomyes = 0	;0=anomalies, 1=absolutes

;*******************************************************************
; FILEPATHS AND VARIABLES 

updatedir  = 'UPDATE20'+strmid(strcompress(edyr,/remove_all),2,2)
workingdir = '/data/local/hadkw/HADCRUH2/'+updatedir

inMardir   = workingdir+'/OTHERDATA/'
inLanddir  = workingdir+'/STATISTICS/GRIDS/'
outdir     = workingdir+'/STATISTICS/GRIDS/'

inLandSeaMask = workingdir+'/OTHERDATA/new_coverpercentjul08.nc'

inMarfilANM   = inMardir+MarineType+'_5x5_monthly_renorm19812010_anomalies_from_daily_both_relax.nc'
inMarfilABS   = inMardir+MarineType+'_5x5_monthly_from_daily_both_relax.nc'
inMarfilCLM   = inMardir+MarineType+'_5x5_monthly_climatology_from_daily_both_relax.nc'
inMarfilSDV   = inMardir+MarineType+'_5x5_monthly_stdev_from_daily_both_relax.nc'

IF (LandType EQ 'standard') THEN BEGIN
  inLandq   = inLanddir+'HadISDH.landq.'+LandVersion+'_FLATgridIDPHA5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
  inLandRH  = inLanddir+'HadISDH.landRH.'+LandVersion+'_FLATgridIDPHA5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
  inLande   = inLanddir+'HadISDH.lande.'+LandVersion+'_FLATgridIDPHA5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
  inLandTw  = inLanddir+'HadISDH.landTw.'+LandVersion+'_FLATgridIDPHA5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
  inLandT   = inLanddir+'HadISDH.landT.'+LandVersion+'_FLATgridIDPHA5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
  inLandTd  = inLanddir+'HadISDH.landTd.'+LandVersion+'_FLATgridPHADPD5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
  inLandDPD = inLanddir+'HadISDH.landDPD.'+LandVersion+'_FLATgridPHA5by5_anoms7605_'+thenmon+thenyear+'_cf.nc'
ENDIF

outMarfilq   = outdir+'HadISDH.marineq.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outMarfile   = outdir+'HadISDH.marinee.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outMarfilTd  = outdir+'HadISDH.marineTd.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outMarfilT   = outdir+'HadISDH.marineT.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outMarfilTw  = outdir+'HadISDH.marineTw.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outMarfilRH  = outdir+'HadISDH.marineRH.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outMarfilDPD = outdir+'HadISDH.marineDPD.'+MarineVersion+'_'+MarineType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'

outBlendfilq   = outdir+'HadISDH.blendq.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outBlendfile   = outdir+'HadISDH.blende.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outBlendfilTd  = outdir+'HadISDH.blendTd.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outBlendfilT   = outdir+'HadISDH.blendT.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outBlendfilTw  = outdir+'HadISDH.blendTw.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outBlendfilRH  = outdir+'HadISDH.blendRH.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'
outBlendfilDPD = outdir+'HadISDH.blendDPD.'+BlendVersion+'_'+BlendType+'_anoms'+climchoice+'_'+nowmon+nowyear+'_cf.nc'

nyrs     = (edyr+1)-styr
nmons    = nyrs*12
int_mons = indgen(nmons)

nlons = 72
nlats = 36
lons  = (findgen(nlons)*5.)-177.5
lats  = (findgen(nlats)*(-5.))+87.5
lats  = REVERSE(lats)

;landarr   = make_array(nlons,nlats,nmons,/float,value=mdi)
;marinearr = make_array(nlons,nlats,nmons,/float,value=mdi)
blendarrANM  = make_array(nlons,nlats,nmons,/float,value=mdi)
blendarrABS  = make_array(nlons,nlats,nmons,/float,value=mdi)

;************************************************************************************
; read in land sea mask
inn   = NCDF_OPEN(inLandSeaMask)
varid = NCDF_VARID(inn,'pct_land')
NCDF_VARGET,inn,varid,pctl
NCDF_CLOSE,inn

;pctl = REVERSE(pctl,2)

; blend using coastal gridboxes as 25-75% weighting
; adjust pctl such that when value is not 0 or 100 it is at least 25 or max 75
; thus land and sea in coastal boxes will always be represented
lows = WHERE(pctl LT 25. AND pctl NE 0.,count)
IF (count GT 0) THEN pctl(lows) = 25.
highs = WHERE(pctl GT 75. AND pctl NE 100.,count)
IF (count GT 0) THEN pctl(highs) = 75.
;stop

; Loop through each variable
for loo = 0,6 DO BEGIN

  ; read in marine ANOMALIES
  inn = NCDF_OPEN(inMarFilANM) 
  CASE loo OF
    0: varid = NCDF_VARID(inn,'specific_humidity_anomalies')		
    1: varid = NCDF_VARID(inn,'vapor_pressure_anomalies')		
    2: varid = NCDF_VARID(inn,'dew_point_temperature_anomalies')	
    3: varid = NCDF_VARID(inn,'marine_air_temperature_anomalies')	
    4: varid = NCDF_VARID(inn,'wet_bulb_temperature_anomalies')		
    5: varid = NCDF_VARID(inn,'relative_humidity_anomalies')		
    6: varid = NCDF_VARID(inn,'dew_point_temperature_anomalies')	
  ENDCASE
  NCDF_VARGET,inn,varid,marinearrANM
  marinearrANM = reverse(marinearrANM,2) ; CHECK THIS IS STILL NECESSARY!
  bads = where(marinearrANM LT -100,count)
  if (count GT 0) THEN marinearrANM(bads) = mdi
  NCDF_CLOSE,inn
  ; read in marine ACTUALS
  inn = NCDF_OPEN(inMarFilABS) 
  CASE loo OF
    0: varid = NCDF_VARID(inn,'specific_humidity')		
    1: varid = NCDF_VARID(inn,'vapor_pressure')		
    2: varid = NCDF_VARID(inn,'dew_point_temperature')	
    3: varid = NCDF_VARID(inn,'marine_air_temperature')	
    4: varid = NCDF_VARID(inn,'wet_bulb_temperature')		
    5: varid = NCDF_VARID(inn,'relative_humidity')		
    6: varid = NCDF_VARID(inn,'dew_point_temperature')	
  ENDCASE
  NCDF_VARGET,inn,varid,marinearrABS
  marinearrABS = reverse(marinearrABS,2) ; CHECK THIS IS STILL NECESSARY!
  bads = where(marinearrABS LT -100,count)
  if (count GT 0) THEN marinearrABS(bads) = mdi
  NCDF_CLOSE,inn
  ; read in marine CLIMATOLOGY (1981-2010)
  inn = NCDF_OPEN(inMarFilCLM) 
  CASE loo OF
    0: varid = NCDF_VARID(inn,'specific_humidity')		
    1: varid = NCDF_VARID(inn,'vapor_pressure')		
    2: varid = NCDF_VARID(inn,'dew_point_temperature')	
    3: varid = NCDF_VARID(inn,'marine_air_temperature')	
    4: varid = NCDF_VARID(inn,'wet_bulb_temperature')		
    5: varid = NCDF_VARID(inn,'relative_humidity')		
    6: varid = NCDF_VARID(inn,'dew_point_temperature')	
  ENDCASE
  NCDF_VARGET,inn,varid,marinearrCLM
  marinearrCLM = reverse(marinearrCLM,2) ; CHECK THIS IS STILL NECESSARY!
  bads = where(marinearrCLM LT -100,count)
  if (count GT 0) THEN marinearrCLM(bads) = mdi
  NCDF_CLOSE,inn
  ; read in marine STDEVS
  inn = NCDF_OPEN(inMarFilSDV) 
  CASE loo OF
    0: varid = NCDF_VARID(inn,'specific_humidity')		
    1: varid = NCDF_VARID(inn,'vapor_pressure')		
    2: varid = NCDF_VARID(inn,'dew_point_temperature')	
    3: varid = NCDF_VARID(inn,'marine_air_temperature')	
    4: varid = NCDF_VARID(inn,'wet_bulb_temperature')		
    5: varid = NCDF_VARID(inn,'relative_humidity')		
    6: varid = NCDF_VARID(inn,'dew_point_temperature')	
  ENDCASE
  NCDF_VARGET,inn,varid,marinearrSDV
  marinearrSDV = reverse(marinearrSDV,2) ; CHECK THIS IS STILL NECESSARY!
  bads = where(marinearrSDV LT -100,count)
  if (count GT 0) THEN marinearrSDV(bads) = mdi
  NCDF_CLOSE,inn

  ;stop

  ; read in HadISDH data and renormalise to 1981-2010
  ; could also mask by where uncertainties are larger than the standard deviation
  CASE loo OF
    0: BEGIN
      inn=NCDF_OPEN(inLandq)
      Avarid=NCDF_VARID(inn,'q_anoms') 
      ABvarid=NCDF_VARID(inn,'q_abs') 
    END
    1: BEGIN
      inn=NCDF_OPEN(inLande)
      Avarid=NCDF_VARID(inn,'e_anoms') 
      ABvarid=NCDF_VARID(inn,'e_abs') 
    END
    2: BEGIN
      inn=NCDF_OPEN(inLandTd)
      Avarid=NCDF_VARID(inn,'td_anoms') 
      ABvarid=NCDF_VARID(inn,'td_abs') 
    END
    3: BEGIN
      inn=NCDF_OPEN(inLandT)
      Avarid=NCDF_VARID(inn,'t_anoms') 
      ABvarid=NCDF_VARID(inn,'t_abs') 
    END
    4: BEGIN
      inn=NCDF_OPEN(inLandTw)
      Avarid=NCDF_VARID(inn,'tw_anoms') 
      ABvarid=NCDF_VARID(inn,'tw_abs') 
    END
    5: BEGIN
      inn=NCDF_OPEN(inLandRH)
      Avarid=NCDF_VARID(inn,'rh_anoms') 
      ABvarid=NCDF_VARID(inn,'rh_abs') 
    END
    6: BEGIN
      inn=NCDF_OPEN(inLandDPD)
      Avarid=NCDF_VARID(inn,'dpd_anoms') 
      ABvarid=NCDF_VARID(inn,'dpd_abs') 
    END
    
  ENDCASE
  ;varid=NCDF_VARID(inn,'q_abs') 
  ;NCDF_VARGET,inn,Avarid,Alandarr
  NCDF_VARGET,inn,Avarid,landarrANM
  NCDF_VARGET,inn,ABvarid,landarrABS
  bads = where(landarrANM LT -100.,count)
  if (count GT 0) THEN landarrANM(bads) = mdi
  bads = where(landarrABS LT -100.,count)
  if (count GT 0) THEN landarrABS(bads) = mdi
  timvarid = NCDF_VARID(inn,'time') 
  NCDF_VARGET,inn,timvarid,landtims
  NCDF_CLOSE,inn
  ;stop

  ; renormalise the land array anomalies to chosen clim if no 7605 or marine array if not 8110?
  FOR i = 0,nlons-1 DO BEGIN
    FOR j = 0,nlats-1 DO BEGIN
      ; SORT OUT WHETHER TO RENORM LAND OR MARINE!!!
      IF (climchoice EQ '7605') THEN subarr = marinearrANM(i,j,*) ELSE subarr = landarrANM(i,j,*)
      gots = WHERE(subarr NE mdi,count)
      IF (count GT nmons/2) THEN BEGIN	; no point if there isn't much data
        monarr = REFORM(subarr,12,nyrs)
        FOR mm=0,11 DO BEGIN
          submon = monarr(mm,*)
	  gots = WHERE(submon NE mdi,count)
	  IF (count GT 20) THEN BEGIN
	    clims = submon(clst:cled)
	    gotclims = WHERE(clims NE mdi,count)
	    IF (count GT 15) THEN submon(gots) = submon(gots)-MEAN(clims(gotclims))
	  ENDIF ELSE submon(*,*) = mdi
	  monarr(mm,*) = submon
        ENDFOR
      IF (climchoice EQ '7605') THEN marinearrANM(i,j,*) = REFORM(monarr,nmons) ELSE landarrANM(i,j,*) = REFORM(monarr,nmons)
      ENDIF ELSE BEGIN
        IF (climchoice EQ '7605') THEN marinearrANM(i,j,*) = mdi ELSE landarrANM(i,j,*) = mdi
      ENDELSE
    ENDFOR
  ENDFOR
  ;stop

  ;Now do the blending
  FOR i=0,nlons-1 DO BEGIN
    FOR j=0,nlats-1 DO BEGIN
      FOR nn=0,nmons-1 DO BEGIN
        IF (pctl(i,j) EQ 100.) THEN BEGIN
          IF (marinearrANM(i,j,nn) EQ mdi) AND (landarrANM(i,j,nn) NE mdi) THEN blendarrANM(i,j,nn) = landarrANM(i,j,nn)
          IF (marinearrANM(i,j,nn) NE mdi) AND (landarrANM(i,j,nn) EQ mdi) THEN blendarrANM(i,j,nn) = marinearrANM(i,j,nn)
	  IF (marinearrANM(i,j,nn) NE mdi) AND (landarrANM(i,j,nn) NE mdi) THEN blendarrANM(i,j,nn) = (marinearrANM(i,j,nn)*0.25)+(landarrANM(i,j,nn)*0.75) 
          IF (marinearrABS(i,j,nn) EQ mdi) AND (landarrABS(i,j,nn) NE mdi) THEN blendarrABS(i,j,nn) = landarrABS(i,j,nn)
          IF (marinearrABS(i,j,nn) NE mdi) AND (landarrABS(i,j,nn) EQ mdi) THEN blendarrABS(i,j,nn) = marinearrABS(i,j,nn)
	  IF (marinearrABS(i,j,nn) NE mdi) AND (landarrABS(i,j,nn) NE mdi) THEN blendarrABS(i,j,nn) = (marinearrABS(i,j,nn)*0.25)+(landarrABS(i,j,nn)*0.75) 
        ENDIF ELSE IF (pctl(i,j) EQ 0.) THEN BEGIN
          IF (marinearrANM(i,j,nn) EQ mdi) AND (landarrANM(i,j,nn) NE mdi) THEN blendarrANM(i,j,nn) = landarrANM(i,j,nn)
          IF (marinearrANM(i,j,nn) NE mdi) AND (landarrANM(i,j,nn) EQ mdi) THEN blendarrANM(i,j,nn) = marinearrANM(i,j,nn)
	  IF (marinearrANM(i,j,nn) NE mdi) AND (landarrANM(i,j,nn) NE mdi) THEN blendarrANM(i,j,nn) = (marinearrANM(i,j,nn)*0.75)+(landarrANM(i,j,nn)*0.25) 
          IF (marinearrABS(i,j,nn) EQ mdi) AND (landarrABS(i,j,nn) NE mdi) THEN blendarrABS(i,j,nn) = landarrABS(i,j,nn)
          IF (marinearrABS(i,j,nn) NE mdi) AND (landarrABS(i,j,nn) EQ mdi) THEN blendarrABS(i,j,nn) = marinearrABS(i,j,nn)
	  IF (marinearrABS(i,j,nn) NE mdi) AND (landarrABS(i,j,nn) NE mdi) THEN blendarrABS(i,j,nn) = (marinearrABS(i,j,nn)*0.75)+(landarrABS(i,j,nn)*0.25) 
        ENDIF ELSE BEGIN
          IF (marinearrANM(i,j,nn) EQ mdi) AND (landarrANM(i,j,nn) NE mdi) THEN blendarrANM(i,j,nn) = landarrANM(i,j,nn)
          IF (marinearrANM(i,j,nn) NE mdi) AND (landarrANM(i,j,nn) EQ mdi) THEN blendarrANM(i,j,nn) = marinearrANM(i,j,nn)
	  IF (marinearrANM(i,j,nn) NE mdi) AND (landarrANM(i,j,nn) NE mdi) THEN blendarrANM(i,j,nn) = (marinearrANM(i,j,nn)*((100.-pctl(i,j))/100.))+(landarrANM(i,j,nn)*(pctl(i,j)/100.))  
          IF (marinearrABS(i,j,nn) EQ mdi) AND (landarrABS(i,j,nn) NE mdi) THEN blendarrABS(i,j,nn) = landarrABS(i,j,nn)
          IF (marinearrABS(i,j,nn) NE mdi) AND (landarrABS(i,j,nn) EQ mdi) THEN blendarrABS(i,j,nn) = marinearrABS(i,j,nn)
	  IF (marinearrABS(i,j,nn) NE mdi) AND (landarrABS(i,j,nn) NE mdi) THEN blendarrABS(i,j,nn) = (marinearrABS(i,j,nn)*((100.-pctl(i,j))/100.))+(landarrABS(i,j,nn)*(pctl(i,j)/100.))  
        ENDELSE
      ENDFOR  
    ENDFOR
  ENDFOR

  ; output gridded Blended product
  CASE loo OF
    0: wilma = NCDF_CREATE(outBlendfilq,/clobber)
    1: wilma = NCDF_CREATE(outBlendfile,/clobber)
    2: wilma = NCDF_CREATE(outBlendfilTd,/clobber)
    3: wilma = NCDF_CREATE(outBlendfilT,/clobber)
    4: wilma = NCDF_CREATE(outBlendfilTw,/clobber)
    5: wilma = NCDF_CREATE(outBlendfilRH,/clobber)
    6: wilma = NCDF_CREATE(outBlendfilDPD,/clobber)
  ENDCASE
  
  tid   = NCDF_DIMDEF(wilma,'time',nmons)
  clmid = NCDF_DIMDEF(wilma,'month',12)
  latid = NCDF_DIMDEF(wilma,'latitude',nlats)
  lonid = NCDF_DIMDEF(wilma,'longitude',nlons)
  
  timesvar = NCDF_VARDEF(wilma,'time',[tid],/SHORT)
  latsvar  = NCDF_VARDEF(wilma,'latitude',[latid],/FLOAT)
  lonsvar  = NCDF_VARDEF(wilma,'longitude',[lonid],/FLOAT)
  CASE loo OF
    0:   anomvar = NCDF_VARDEF(wilma,'q_anoms',[lonid,latid,tid],/FLOAT)
    1:   anomvar = NCDF_VARDEF(wilma,'e_anoms',[lonid,latid,tid],/FLOAT)
    2:   anomvar = NCDF_VARDEF(wilma,'td_anoms',[lonid,latid,tid],/FLOAT)
    3:   anomvar = NCDF_VARDEF(wilma,'t_anoms',[lonid,latid,tid],/FLOAT)
    4:   anomvar = NCDF_VARDEF(wilma,'tw_anoms',[lonid,latid,tid],/FLOAT)
    5:   anomvar = NCDF_VARDEF(wilma,'rh_anoms',[lonid,latid,tid],/FLOAT)
    6:   anomvar = NCDF_VARDEF(wilma,'dpd_anoms',[lonid,latid,tid],/FLOAT)
  ENDCASE
  CASE loo OF
    0:   absvar = NCDF_VARDEF(wilma,'q_abs',[lonid,latid,tid],/FLOAT)
    1:   absvar = NCDF_VARDEF(wilma,'e_abs',[lonid,latid,tid],/FLOAT)
    2:   absvar = NCDF_VARDEF(wilma,'td_abs',[lonid,latid,tid],/FLOAT)
    3:   absvar = NCDF_VARDEF(wilma,'t_abs',[lonid,latid,tid],/FLOAT)
    4:   absvar = NCDF_VARDEF(wilma,'tw_abs',[lonid,latid,tid],/FLOAT)
    5:   absvar = NCDF_VARDEF(wilma,'rh_abs',[lonid,latid,tid],/FLOAT)
    6:   absvar = NCDF_VARDEF(wilma,'dpd_abs',[lonid,latid,tid],/FLOAT)
  ENDCASE
  
  NCDF_ATTPUT,wilma,'time','long_name','time'
  NCDF_ATTPUT,wilma,'time','units','months beginning Jan 1973'
  NCDF_ATTPUT,wilma,'time','axis','T'
  NCDF_ATTPUT,wilma,'time','calendar','gregorian'
  NCDF_ATTPUT,wilma,'time','valid_min',0.
  NCDF_ATTPUT,wilma,'latitude','long_name','Latitude'
  NCDF_ATTPUT,wilma,'latitude','units','Degrees'
  NCDF_ATTPUT,wilma,'latitude','valid_min',-90.
  NCDF_ATTPUT,wilma,'latitude','valid_max',90.
  NCDF_ATTPUT,wilma,'longitude','long_name','Longitude'
  NCDF_ATTPUT,wilma,'longitude','units','Degrees'
  NCDF_ATTPUT,wilma,'longitude','valid_min',-180.
  NCDF_ATTPUT,wilma,'longitude','valid_max',180.

  CASE loo OF
    0: BEGIN    
      NCDF_ATTPUT,wilma,'q_anoms','long_name','Blended Specific Humidity monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'q_anoms','units','g/kg'
      NCDF_ATTPUT,wilma,'q_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'q_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'q_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'q_anoms','missing_value',mdi
    END
    1: BEGIN    
      NCDF_ATTPUT,wilma,'e_anoms','long_name','Blended Vapour Pressure monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'e_anoms','units','hPa'
      NCDF_ATTPUT,wilma,'e_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'e_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'e_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'e_anoms','missing_value',mdi
    END
    2: BEGIN    
      NCDF_ATTPUT,wilma,'td_anoms','long_name','Blended Dew Point Temperature monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'td_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'td_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'td_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'td_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'td_anoms','missing_value',mdi
    END
    3: BEGIN    
      NCDF_ATTPUT,wilma,'t_anoms','long_name','Blended Air Temperature monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'t_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'t_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'t_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'t_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'t_anoms','missing_value',mdi
    END
    4: BEGIN    
      NCDF_ATTPUT,wilma,'tw_anoms','long_name','Blended Wet Bulb Temperature monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'tw_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'tw_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'tw_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'tw_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'tw_anoms','missing_value',mdi
    END
    5: BEGIN    
      NCDF_ATTPUT,wilma,'rh_anoms','long_name','Blended Relative Humidity monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'rh_anoms','units','%rh'
      NCDF_ATTPUT,wilma,'rh_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'rh_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'rh_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'rh_anoms','missing_value',mdi
    END
    6: BEGIN    
      NCDF_ATTPUT,wilma,'dpd_anoms','long_name','Blended Dew Point Depression monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'dpd_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'dpd_anoms','axis','T'
      valid=WHERE(blendarrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrANM(valid))
        max_t=MAX(blendarrANM(valid))
        NCDF_ATTPUT,wilma,'dpd_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'dpd_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'dpd_anoms','missing_value',mdi
    END

  ENDCASE

  CASE loo OF
    0: BEGIN    
      NCDF_ATTPUT,wilma,'q_abs','long_name','Blended Specific Humidity monthly mean'
      NCDF_ATTPUT,wilma,'q_abs','units','g/kg'
      NCDF_ATTPUT,wilma,'q_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'q_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'q_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'q_abs','missing_value',mdi
    END
    1: BEGIN    
      NCDF_ATTPUT,wilma,'e_abs','long_name','Blended Vapour Pressure monthly mean'
      NCDF_ATTPUT,wilma,'e_abs','units','hPa'
      NCDF_ATTPUT,wilma,'e_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'e_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'e_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'e_abs','missing_value',mdi
    END
    2: BEGIN    
      NCDF_ATTPUT,wilma,'td_abs','long_name','Blended Dew Point Temperature monthly mean'
      NCDF_ATTPUT,wilma,'td_abs','units','deg C'
      NCDF_ATTPUT,wilma,'td_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'td_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'td_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'td_abs','missing_value',mdi
    END
    3: BEGIN    
      NCDF_ATTPUT,wilma,'t_abs','long_name','Blended Air Temperature monthly mean'
      NCDF_ATTPUT,wilma,'t_abs','units','deg C'
      NCDF_ATTPUT,wilma,'t_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'t_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'t_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'t_abs','missing_value',mdi
    END
    4: BEGIN    
      NCDF_ATTPUT,wilma,'tw_abs','long_name','Blended Wet Bulb Temperature monthly mean'
      NCDF_ATTPUT,wilma,'tw_abs','units','deg C'
      NCDF_ATTPUT,wilma,'tw_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'tw_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'tw_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'tw_abs','missing_value',mdi
    END
    5: BEGIN    
      NCDF_ATTPUT,wilma,'rh_abs','long_name','Blended Relative Humidity monthly mean'
      NCDF_ATTPUT,wilma,'rh_abs','units','%rh'
      NCDF_ATTPUT,wilma,'rh_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'rh_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'rh_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'rh_abs','missing_value',mdi
    END
    6: BEGIN    
      NCDF_ATTPUT,wilma,'dpd_abs','long_name','Blended Dew Point Depression monthly mean'
      NCDF_ATTPUT,wilma,'dpd_abs','units','deg C'
      NCDF_ATTPUT,wilma,'dpd_abs','axis','T'
      valid=WHERE(blendarrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(blendarrABS(valid))
        max_t=MAX(blendarrABS(valid))
        NCDF_ATTPUT,wilma,'dpd_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'dpd_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'dpd_abs','missing_value',mdi
    END

  ENDCASE

  current_time=SYSTIME()
  NCDF_ATTPUT,wilma,/GLOBAL,'file_created',STRING(current_time)
  
  NCDF_ATTPUT,wilma,/GLOBAL,'description',"HadISDH monthly mean blended land and marine surface climate monitoring product from 1973 onwards. "+$
                                         "Quality control, homogenisation (land), bias correction (marine), uncertainty estimation, averaging over gridboxes (no smoothing "+$
					 "or interpolation). This combines HadISDH.land."+LandVersion+" and HadISDH.marine."+MarineVersion
  NCDF_ATTPUT,wilma,/GLOBAL,'title',"HadISDH monthly mean blended land and marine surface climate monitoring product from 1973 onwards."
  NCDF_ATTPUT,wilma,/GLOBAL,'institution',"Met Office Hadley Centre (UK), National Climatic Data Centre (USA), Climatic Research Unit (UK), "+$
                                         "National Physical Laboratory (UK), Maynooth University (Ireland), National Oceanography Centre Southampton (NOCS)"
  NCDF_ATTPUT,wilma,/GLOBAL,'history',"Updated "+STRING(current_time)
  NCDF_ATTPUT,wilma,/GLOBAL,'source',"HadISD.2.1.0.2016p (Dunn et al., 2016), www.metoffice.gov.uk/hadobs/hadisd/ and ICOADS.3.0.0 and ICOADS.3.0.1 (Freemand et al., 2016), icoads.noaa.gov"
  NCDF_ATTPUT,wilma,/GLOBAL,'comment'," "
  NCDF_ATTPUT,wilma,/GLOBAL,'reference',"Willett, K. M., Dunn, R. J. H., Thorne, P. W., Bell, S., de Podesta, M., Parker, D. E., Jones, "+$
                                       "P. D., and Williams Jr., C. N.: HadISDH land surface multi-variable humidity and temperature "+$
				       "record for climate monitoring, Clim. Past, 10, 1983-2006, doi:10.5194/cp-10-1983-2014, 2014. and Willett et al. in prep."
  NCDF_ATTPUT,wilma,/GLOBAL,'version',"HadISDH.blend."+BlendVersion
  NCDF_ATTPUT,wilma,/GLOBAL,'Conventions',"CF-1.0"

  NCDF_CONTROL,wilma,/ENDEF

  NCDF_VARPUT, wilma,timesvar, landtims
  NCDF_VARPUT, wilma,latsvar, lats
  NCDF_VARPUT, wilma,lonsvar, lons
  NCDF_VARPUT, wilma,anomvar,blendarrANM
  NCDF_VARPUT, wilma,absvar,blendarrABS

  NCDF_CLOSE,wilma

  ; output gridded Marine product
  CASE loo OF
    0: wilma = NCDF_CREATE(outMarfilq,/clobber)
    1: wilma = NCDF_CREATE(outMarfile,/clobber)
    2: wilma = NCDF_CREATE(outMarfilTd,/clobber)
    3: wilma = NCDF_CREATE(outMarfilT,/clobber)
    4: wilma = NCDF_CREATE(outMarfilTw,/clobber)
    5: wilma = NCDF_CREATE(outMarfilRH,/clobber)
    6: wilma = NCDF_CREATE(outMarfilDPD,/clobber)
  ENDCASE
  
  tid   = NCDF_DIMDEF(wilma,'time',nmons)
  clmid = NCDF_DIMDEF(wilma,'month',12)
  latid = NCDF_DIMDEF(wilma,'latitude',nlats)
  lonid = NCDF_DIMDEF(wilma,'longitude',nlons)
  
  timesvar = NCDF_VARDEF(wilma,'time',[tid],/SHORT)
  latsvar  = NCDF_VARDEF(wilma,'latitude',[latid],/FLOAT)
  lonsvar  = NCDF_VARDEF(wilma,'longitude',[lonid],/FLOAT)
  CASE loo OF
    0:   anomvar = NCDF_VARDEF(wilma,'q_anoms',[lonid,latid,tid],/FLOAT)
    1:   anomvar = NCDF_VARDEF(wilma,'e_anoms',[lonid,latid,tid],/FLOAT)
    2:   anomvar = NCDF_VARDEF(wilma,'td_anoms',[lonid,latid,tid],/FLOAT)
    3:   anomvar = NCDF_VARDEF(wilma,'t_anoms',[lonid,latid,tid],/FLOAT)
    4:   anomvar = NCDF_VARDEF(wilma,'tw_anoms',[lonid,latid,tid],/FLOAT)
    5:   anomvar = NCDF_VARDEF(wilma,'rh_anoms',[lonid,latid,tid],/FLOAT)
    6:   anomvar = NCDF_VARDEF(wilma,'dpd_anoms',[lonid,latid,tid],/FLOAT)
  ENDCASE
  CASE loo OF
    0:   absvar = NCDF_VARDEF(wilma,'q_abs',[lonid,latid,tid],/FLOAT)
    1:   absvar = NCDF_VARDEF(wilma,'e_abs',[lonid,latid,tid],/FLOAT)
    2:   absvar = NCDF_VARDEF(wilma,'td_abs',[lonid,latid,tid],/FLOAT)
    3:   absvar = NCDF_VARDEF(wilma,'t_abs',[lonid,latid,tid],/FLOAT)
    4:   absvar = NCDF_VARDEF(wilma,'tw_abs',[lonid,latid,tid],/FLOAT)
    5:   absvar = NCDF_VARDEF(wilma,'rh_abs',[lonid,latid,tid],/FLOAT)
    6:   absvar = NCDF_VARDEF(wilma,'dpd_abs',[lonid,latid,tid],/FLOAT)
  ENDCASE
  CASE loo OF
    0:   climvar = NCDF_VARDEF(wilma,'q_clims',[lonid,latid,clmid],/FLOAT)
    1:   climvar = NCDF_VARDEF(wilma,'e_clims',[lonid,latid,clmid],/FLOAT)
    2:   climvar = NCDF_VARDEF(wilma,'td_clims',[lonid,latid,clmid],/FLOAT)
    3:   climvar = NCDF_VARDEF(wilma,'t_clims',[lonid,latid,clmid],/FLOAT)
    4:   climvar = NCDF_VARDEF(wilma,'tw_clims',[lonid,latid,clmid],/FLOAT)
    5:   climvar = NCDF_VARDEF(wilma,'rh_clims',[lonid,latid,clmid],/FLOAT)
    6:   climvar = NCDF_VARDEF(wilma,'dpd_clims',[lonid,latid,clmid],/FLOAT)
  ENDCASE
  CASE loo OF
    0:   stdvar = NCDF_VARDEF(wilma,'q_std',[lonid,latid,clmid],/FLOAT)
    1:   stdvar = NCDF_VARDEF(wilma,'e_std',[lonid,latid,clmid],/FLOAT)
    2:   stdvar = NCDF_VARDEF(wilma,'td_std',[lonid,latid,clmid],/FLOAT)
    3:   stdvar = NCDF_VARDEF(wilma,'t_std',[lonid,latid,clmid],/FLOAT)
    4:   stdvar = NCDF_VARDEF(wilma,'tw_std',[lonid,latid,clmid],/FLOAT)
    5:   stdvar = NCDF_VARDEF(wilma,'rh_std',[lonid,latid,clmid],/FLOAT)
    6:   stdvar = NCDF_VARDEF(wilma,'dpd_std',[lonid,latid,clmid],/FLOAT)
  ENDCASE
  
  NCDF_ATTPUT,wilma,'time','long_name','time'
  NCDF_ATTPUT,wilma,'time','units','months beginning Jan 1973'
  NCDF_ATTPUT,wilma,'time','axis','T'
  NCDF_ATTPUT,wilma,'time','calendar','gregorian'
  NCDF_ATTPUT,wilma,'time','valid_min',0.
  NCDF_ATTPUT,wilma,'latitude','long_name','Latitude'
  NCDF_ATTPUT,wilma,'latitude','units','Degrees'
  NCDF_ATTPUT,wilma,'latitude','valid_min',-90.
  NCDF_ATTPUT,wilma,'latitude','valid_max',90.
  NCDF_ATTPUT,wilma,'longitude','long_name','Longitude'
  NCDF_ATTPUT,wilma,'longitude','units','Degrees'
  NCDF_ATTPUT,wilma,'longitude','valid_min',-180.
  NCDF_ATTPUT,wilma,'longitude','valid_max',180.

  CASE loo OF
    0: BEGIN    
      NCDF_ATTPUT,wilma,'q_anoms','long_name','Marine Specific Humidity monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'q_anoms','units','g/kg'
      NCDF_ATTPUT,wilma,'q_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'q_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'q_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'q_anoms','missing_value',mdi
    END
    1: BEGIN    
      NCDF_ATTPUT,wilma,'e_anoms','long_name','Marine Vapour Pressure monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'e_anoms','units','hPa'
      NCDF_ATTPUT,wilma,'e_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'e_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'e_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'e_anoms','missing_value',mdi
    END
    2: BEGIN    
      NCDF_ATTPUT,wilma,'td_anoms','long_name','Marine Dew Point Temperature monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'td_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'td_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'td_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'td_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'td_anoms','missing_value',mdi
    END
    3: BEGIN    
      NCDF_ATTPUT,wilma,'t_anoms','long_name','Marine Air Temperature monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'t_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'t_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'t_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'t_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'t_anoms','missing_value',mdi
    END
    4: BEGIN    
      NCDF_ATTPUT,wilma,'tw_anoms','long_name','Marine Wet Bulb Temperature monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'tw_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'tw_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'tw_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'tw_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'tw_anoms','missing_value',mdi
    END
    5: BEGIN    
      NCDF_ATTPUT,wilma,'rh_anoms','long_name','Marine Relative Humidity monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'rh_anoms','units','%rh'
      NCDF_ATTPUT,wilma,'rh_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'rh_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'rh_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'rh_anoms','missing_value',mdi
    END
    6: BEGIN    
      NCDF_ATTPUT,wilma,'dpd_anoms','long_name','Marine Dew Point Depression monthly mean anomaly ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'dpd_anoms','units','deg C'
      NCDF_ATTPUT,wilma,'dpd_anoms','axis','T'
      valid=WHERE(marinearrANM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrANM(valid))
        max_t=MAX(marinearrANM(valid))
        NCDF_ATTPUT,wilma,'dpd_anoms','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'dpd_anoms','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'dpd_anoms','missing_value',mdi
    END

  ENDCASE

  CASE loo OF
    0: BEGIN    
      NCDF_ATTPUT,wilma,'q_abs','long_name','Marine Specific Humidity monthly mean'
      NCDF_ATTPUT,wilma,'q_abs','units','g/kg'
      NCDF_ATTPUT,wilma,'q_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'q_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'q_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'q_abs','missing_value',mdi
    END
    1: BEGIN    
      NCDF_ATTPUT,wilma,'e_abs','long_name','Marine Vapour Pressure monthly mean'
      NCDF_ATTPUT,wilma,'e_abs','units','hPa'
      NCDF_ATTPUT,wilma,'e_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'e_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'e_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'e_abs','missing_value',mdi
    END
    2: BEGIN    
      NCDF_ATTPUT,wilma,'td_abs','long_name','Marine Dew Point Temperature monthly mean'
      NCDF_ATTPUT,wilma,'td_abs','units','deg C'
      NCDF_ATTPUT,wilma,'td_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'td_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'td_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'td_abs','missing_value',mdi
    END
    3: BEGIN    
      NCDF_ATTPUT,wilma,'t_abs','long_name','Marine Air Temperature monthly mean'
      NCDF_ATTPUT,wilma,'t_abs','units','deg C'
      NCDF_ATTPUT,wilma,'t_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'t_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'t_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'t_abs','missing_value',mdi
    END
    4: BEGIN    
      NCDF_ATTPUT,wilma,'tw_abs','long_name','Marine Wet Bulb Temperature monthly mean'
      NCDF_ATTPUT,wilma,'tw_abs','units','deg C'
      NCDF_ATTPUT,wilma,'tw_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'tw_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'tw_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'tw_abs','missing_value',mdi
    END
    5: BEGIN    
      NCDF_ATTPUT,wilma,'rh_abs','long_name','Marine Relative Humidity monthly mean'
      NCDF_ATTPUT,wilma,'rh_abs','units','%rh'
      NCDF_ATTPUT,wilma,'rh_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'rh_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'rh_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'rh_abs','missing_value',mdi
    END
    6: BEGIN    
      NCDF_ATTPUT,wilma,'dpd_abs','long_name','Marine Dew Point Depression monthly mean'
      NCDF_ATTPUT,wilma,'dpd_abs','units','deg C'
      NCDF_ATTPUT,wilma,'dpd_abs','axis','T'
      valid=WHERE(marinearrABS NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrABS(valid))
        max_t=MAX(marinearrABS(valid))
        NCDF_ATTPUT,wilma,'dpd_abs','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'dpd_abs','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'dpd_abs','missing_value',mdi
    END

  ENDCASE

  CASE loo OF
    0: BEGIN    
      NCDF_ATTPUT,wilma,'q_clims','long_name','Marine Specific Humidity monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'q_clims','units','g/kg'
      NCDF_ATTPUT,wilma,'q_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'q_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'q_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'q_clims','missing_value',mdi
    END
    1: BEGIN    
      NCDF_ATTPUT,wilma,'e_clims','long_name','Marine Vapour Pressure monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'e_clims','units','hPa'
      NCDF_ATTPUT,wilma,'e_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'e_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'e_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'e_clims','missing_value',mdi
    END
    2: BEGIN    
      NCDF_ATTPUT,wilma,'td_clims','long_name','Marine Dew Point Temperature monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'td_clims','units','deg C'
      NCDF_ATTPUT,wilma,'td_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'td_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'td_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'td_clims','missing_value',mdi
    END
    3: BEGIN    
      NCDF_ATTPUT,wilma,'t_clims','long_name','Marine Air Temperature monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'t_clims','units','deg C'
      NCDF_ATTPUT,wilma,'t_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'t_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'t_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'t_clims','missing_value',mdi
    END
    4: BEGIN    
      NCDF_ATTPUT,wilma,'tw_clims','long_name','Marine Wet Bulb Temperature monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'tw_clims','units','deg C'
      NCDF_ATTPUT,wilma,'tw_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'tw_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'tw_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'tw_clims','missing_value',mdi
    END
    5: BEGIN    
      NCDF_ATTPUT,wilma,'rh_clims','long_name','Marine Relative Humidity monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'rh_clims','units','%rh'
      NCDF_ATTPUT,wilma,'rh_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'rh_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'rh_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'rh_clims','missing_value',mdi
    END
    6: BEGIN    
      NCDF_ATTPUT,wilma,'dpd_clims','long_name','Marine Dew Point Depression monthly mean climatology ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'dpd_clims','units','deg C'
      NCDF_ATTPUT,wilma,'dpd_clims','axis','T'
      valid=WHERE(marinearrCLM NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrCLM(valid))
        max_t=MAX(marinearrCLM(valid))
        NCDF_ATTPUT,wilma,'dpd_clims','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'dpd_clims','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'dpd_clims','missing_value',mdi
    END

  ENDCASE

  CASE loo OF
    0: BEGIN    
      NCDF_ATTPUT,wilma,'q_std','long_name','Marine Specific Humidity monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'q_std','units','g/kg'
      NCDF_ATTPUT,wilma,'q_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'q_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'q_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'q_std','missing_value',mdi
    END
    1: BEGIN    
      NCDF_ATTPUT,wilma,'e_std','long_name','Marine Vapour Pressure monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'e_std','units','hPa'
      NCDF_ATTPUT,wilma,'e_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'e_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'e_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'e_std','missing_value',mdi
    END
    2: BEGIN    
      NCDF_ATTPUT,wilma,'td_std','long_name','Marine Dew Point Temperature monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'td_std','units','deg C'
      NCDF_ATTPUT,wilma,'td_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'td_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'td_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'td_std','missing_value',mdi
    END
    3: BEGIN    
      NCDF_ATTPUT,wilma,'t_std','long_name','Marine Air Temperature monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'t_std','units','deg C'
      NCDF_ATTPUT,wilma,'t_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'t_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'t_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'t_std','missing_value',mdi
    END
    4: BEGIN    
      NCDF_ATTPUT,wilma,'tw_std','long_name','Marine Wet Bulb Temperature monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'tw_std','units','deg C'
      NCDF_ATTPUT,wilma,'tw_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'tw_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'tw_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'tw_std','missing_value',mdi
    END
    5: BEGIN    
      NCDF_ATTPUT,wilma,'rh_std','long_name','Marine Relative Humidity monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'rh_std','units','%rh'
      NCDF_ATTPUT,wilma,'rh_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'rh_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'rh_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'rh_std','missing_value',mdi
    END
    6: BEGIN    
      NCDF_ATTPUT,wilma,'dpd_std','long_name','Marine Dew Point Depression monthly mean climatological standard deviation ('+strcompress(clstyr,/remove_all)+'-'+strcompress(cledyr,/remove_all)+')'
      NCDF_ATTPUT,wilma,'dpd_std','units','deg C'
      NCDF_ATTPUT,wilma,'dpd_std','axis','T'
      valid=WHERE(marinearrSDV NE mdi, tc)
      IF tc GE 1 THEN BEGIN
        min_t=MIN(marinearrSDV(valid))
        max_t=MAX(marinearrSDV(valid))
        NCDF_ATTPUT,wilma,'dpd_std','valid_min',min_t(0)
        NCDF_ATTPUT,wilma,'dpd_std','valid_max',max_t(0)
      ENDIF
      NCDF_ATTPUT,wilma,'dpd_std','missing_value',mdi
    END

  ENDCASE

  current_time=SYSTIME()
  NCDF_ATTPUT,wilma,/GLOBAL,'file_created',STRING(current_time)
  NCDF_ATTPUT,wilma,/GLOBAL,'description',"HadISDH monthly mean marine surface climate monitoring product from 1973 onwards. "+$
                                         "Quality control, bias correction, uncertainty estimation, averaging over gridboxes (no smoothing "+$
					 "or interpolation)."
  NCDF_ATTPUT,wilma,/GLOBAL,'title',"HadISDH monthly mean marine surface climate monitoring product from 1973 onwards."
  NCDF_ATTPUT,wilma,/GLOBAL,'institution',"Met Office Hadley Centre (UK), National Oceanography Centre Southampton (NOCS)"
  NCDF_ATTPUT,wilma,/GLOBAL,'history',"Updated "+STRING(current_time)
  NCDF_ATTPUT,wilma,/GLOBAL,'source',"ICOADS.3.0.0 and ICOADS.3.0.1 (Freeman et al., 2016), http://icoads.noaa.gov/"
  NCDF_ATTPUT,wilma,/GLOBAL,'comment'," "
  NCDF_ATTPUT,wilma,/GLOBAL,'reference',"Willett, K. M., Dunn, R. J. H., Parker, D. E., Kennedy, J. F. and Berry, D. I.: "+$
                                       "HadISDH marine surface multi-variable humidity and temperature "+$
				       "record for climate monitoring, in prep."
  NCDF_ATTPUT,wilma,/GLOBAL,'version',"HadISDH.marine."+MarineVersion
  NCDF_ATTPUT,wilma,/GLOBAL,'Conventions',"CF-1.0"

  NCDF_CONTROL,wilma,/ENDEF

  NCDF_VARPUT, wilma,timesvar, landtims
  NCDF_VARPUT, wilma,latsvar, lats
  NCDF_VARPUT, wilma,lonsvar, lons
  NCDF_VARPUT, wilma,anomvar,marinearrANM
  NCDF_VARPUT, wilma,absvar,marinearrABS
  NCDF_VARPUT, wilma,climvar,marinearrCLM
  NCDF_VARPUT, wilma,stdvar,marinearrSDV

  NCDF_CLOSE,wilma


ENDFOR

stop
end
