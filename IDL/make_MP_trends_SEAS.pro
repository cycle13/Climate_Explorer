pro make_MP_trends_SEAS

; to calculate trends - NO missing data tolerance - need a land/sea mask for that

;based on kate_timeseries.pro - WAVE for ICOADS

; written by Kate Willett
; last modified JUL 2012

; could use gridbox standard deviations or uncertainties
;   fit trend to val+unc and val-unc - look at differences
;   bootstrap 100+ timeseries randomly picking out values between val+/-unc, fit trends
;   ENSURE THERE IS A MAXMAX AND MINMIN timeseries?

;---------------------------------------------------
; set up directories and filenames
param='dpd'	;'dpd','td','t','tw','e','q','rh'
param2='DPD'	;'DPD','Td','T','Tw','e','q','RH'
nowmon='JAN'
nowyear='2015'
thenmon='JAN'
thenyear='2015'
homogtype='PHA'	;'PHA','ID','DPD', 'RAW'
version='2.0.1.2014p'

mdi=-1e+30

styr=1973
edyr=2014
sttrd=1973  	; years over which to calculate trends
edtrd=2014  	;
nyrs=(edyr+1)-styr
nmons=nyrs*12
int_mons=indgen(nmons)
latlg=5.
lonlg=5.
stlt=-90+(latlg/2.)
stln=-180+(lonlg/2.)
nlats=180/latlg
nlons=360/lonlg
nbox=LONG(nlats*nlons)

dir='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/'

CASE param OF
  'dpd': BEGIN
    IF (homogtype EQ 'PHA') THEN infile='HadISDH.landDPD.'+version+'_FLATgridPHA5by5_'+thenmon+thenyear
    IF (homogtype EQ 'RAW') THEN infile='HadISDH.landDPD.'+version+'_FLATgridRAW5by5_'+thenmon+thenyear
    units='deg C'
  END
  'td': BEGIN
    IF (homogtype EQ 'PHA') THEN infile='HadISDH.landTd.'+version+'_FLATgridPHA5by5_'+thenmon+thenyear
    IF (homogtype EQ 'DPD') THEN infile='HadISDH.landTd.'+version+'_FLATgridPHADPD5by5_'+thenmon+thenyear
    IF (homogtype EQ 'RAW') THEN infile='HadISDH.landTd.'+version+'_FLATgridRAW5by5_'+thenmon+thenyear
    units='deg C'
  END
  't': BEGIN
    IF (homogtype EQ 'ID') THEN infile='HadISDH.landT.'+version+'_FLATgridIDPHA5by5_'+thenmon+thenyear
    IF (homogtype EQ 'RAW') THEN infile='HadISDH.landT.'+version+'_FLATgridRAW5by5_'+thenmon+thenyear
    units='deg C'
  END
  'tw': BEGIN
    IF (homogtype EQ 'ID') THEN infile='HadISDH.landTw.'+version+'_FLATgridIDPHA5by5_'+thenmon+thenyear
    IF (homogtype EQ 'RAW') THEN infile='HadISDH.landTw.'+version+'_FLATgridRAW5by5_'+thenmon+thenyear
    units='deg C'
  END
  'q': BEGIN
    IF (homogtype EQ 'ID') THEN infile='HadISDH.landq.'+version+'_FLATgridIDPHA5by5_'+thenmon+thenyear
    IF (homogtype EQ 'RAW') THEN infile='HadISDH.landq.'+version+'_FLATgridRAW5by5_'+thenmon+thenyear
    units='g/kg'
  END
  'e': BEGIN
    IF (homogtype EQ 'ID') THEN infile='HadISDH.lande.'+version+'_FLATgridIDPHA5by5_'+thenmon+thenyear
    IF (homogtype EQ 'RAW') THEN infile='HadISDH.lande.'+version+'_FLATgridRAW5by5_'+thenmon+thenyear
    units='hPa'
  END
  'rh': BEGIN
    IF (homogtype EQ 'ID') THEN infile='HadISDH.landRH.'+version+'_FLATgridIDPHA5by5_'+thenmon+thenyear
    IF (homogtype EQ 'RAW') THEN infile='HadISDH.landRH.'+version+'_FLATgridRAW5by5_'+thenmon+thenyear
    units='%rh'
  END
  
ENDCASE
ofile=infile+'_MPtrendsSEAS_'+strcompress(styr,/remove_all)+strcompress(edyr,/remove_all)	;70S-70N

mdi=-1e+30

;styr=1973
;edyr=2014
;sttrd=1973  	; years over which to calculate trends
;edtrd=2014  	;
;nyrs=(edyr+1)-styr
;nmons=nyrs*12
;int_mons=indgen(nmons)
;latlg=5.
;lonlg=5.
;stlt=-90+(latlg/2.)
;stln=-180+(lonlg/2.)
;nlats=180/latlg
;nlons=360/lonlg
;nbox=LONG(nlats*nlons)

trdst=(sttrd-styr)*12	; zero if 1973
trded=(((edtrd+1)-styr)*12)-1	; zero if 1973
;setting trend years and mons which will replace actual
nyrs=(edtrd+1)-sttrd
nmons=nyrs*12

print,trdst,trded,nyrs,nmons


lats=(findgen(nlats)*latlg)+stlt
lons=(findgen(nlons)*lonlg)+stln

djf=[0,1]
mam=[2,3,4]
jja=[5,6,7]
son=[8,9,10]

trendJJA=make_array(nlons,nlats,/float,value=mdi)	;trend
trendJJAU=make_array(nlons,nlats,/float,value=mdi)	;95th pct
trendJJAL=make_array(nlons,nlats,/float,value=mdi)	;5th pct

trendSON=make_array(nlons,nlats,/float,value=mdi)	;trend
trendSONU=make_array(nlons,nlats,/float,value=mdi)	;95th pct
trendSONL=make_array(nlons,nlats,/float,value=mdi)	;5th pct

trendDJF=make_array(nlons,nlats,/float,value=mdi)	;trend
trendDJFU=make_array(nlons,nlats,/float,value=mdi)	;95th pct
trendDJFL=make_array(nlons,nlats,/float,value=mdi)	;5th pct

trendMAM=make_array(nlons,nlats,/float,value=mdi)	;trend
trendMAMU=make_array(nlons,nlats,/float,value=mdi)	;95th pct
trendMAML=make_array(nlons,nlats,/float,value=mdi)	;5th pct

qseas=make_array(nlons,nlats,nyrs,4,/float,value=mdi)	;0=DJF,1=MAM,2=JJA,3=SON
;----------------------------------------------------
; read in files

filee=NCDF_OPEN(dir+'GRIDS/'+infile+'_cf.nc')
timvarid=NCDF_VARID(filee,'time')
longs_varid=NCDF_VARID(filee,'longitude')
lats_varid=NCDF_VARID(filee,'latitude')
CASE param OF
  'q': qvarid=NCDF_VARID(filee,'q_anoms')
  'e': qvarid=NCDF_VARID(filee,'e_anoms')
  'rh': qvarid=NCDF_VARID(filee,'rh_anoms')
  't': qvarid=NCDF_VARID(filee,'t_anoms')
  'td': qvarid=NCDF_VARID(filee,'td_anoms')
  'tw': qvarid=NCDF_VARID(filee,'tw_anoms')
  'dpd': qvarid=NCDF_VARID(filee,'dpd_anoms')
ENDCASE
NCDF_VARGET,filee,timvarid,times
NCDF_VARGET,filee,qvarid,q_values
NCDF_VARGET,filee,lats_varid,lats
NCDF_VARGET,filee,longs_varid,longs
NCDF_CLOSE,filee

q_values=q_values(*,*,trdst:trded)

FOR lnn=0,nlons-1 DO BEGIN
  FOR ltt=0,nlats-1 DO BEGIN
    print,lnn,ltt
    subarr=q_values(lnn,ltt,*)
    subarr=REFORM(subarr,12,nyrs)
    FOR yy=0,nyrs-1 DO BEGIN
      IF (yy GT 0) THEN BEGIN
        gots=WHERE(subarr(djf,yy) NE mdi,count)
        IF (count GT 0) THEN BEGIN
	  tots=TOTAL(subarr(djf(gots),yy)) 
	  IF (subarr(11,yy-1) NE mdi) THEN BEGIN
	    tots=tots+subarr(11,yy-1)
	    count=count+1
	  ENDIF
	  IF (count GE 2) THEN qseas(lnn,ltt,yy,0)=tots/count
	ENDIF
      ENDIF
      gots=WHERE(subarr(mam,yy) NE mdi,count)
      IF (count GE 2) THEN qseas(lnn,ltt,yy,1)=TOTAL(subarr(mam(gots),yy)/count)
      gots=WHERE(subarr(jja,yy) NE mdi,count)
      IF (count GE 2) THEN qseas(lnn,ltt,yy,2)=TOTAL(subarr(jja(gots),yy)/count)
      gots=WHERE(subarr(son,yy) NE mdi,count)
      IF (count GE 2) THEN qseas(lnn,ltt,yy,3)=TOTAL(subarr(son(gots),yy)/count)
    ENDFOR
    
    gots=WHERE(qseas(lnn,ltt,*,0) NE mdi,countq)
    IF (countq GT nyrs*0.7) THEN BEGIN
      trendDJF(lnn,ltt)=(median_pairwise(qseas(lnn,ltt,*,0),mdi,se,lc,uc))*10.	;decadal
      trendDJFU(lnn,ltt)=uc*10.
      trendDJFL(lnn,ltt)=lc*10.
    ENDIF
    gots=WHERE(qseas(lnn,ltt,*,1) NE mdi,countq)
    IF (countq GT nyrs*0.7) THEN BEGIN
      trendMAM(lnn,ltt)=(median_pairwise(qseas(lnn,ltt,*,1),mdi,se,lc,uc))*10.	;decadal
      trendMAMU(lnn,ltt)=uc*10.
      trendMAML(lnn,ltt)=lc*10.
    ENDIF
    gots=WHERE(qseas(lnn,ltt,*,2) NE mdi,countq)
    IF (countq GT nyrs*0.7) THEN BEGIN
      trendJJA(lnn,ltt)=(median_pairwise(qseas(lnn,ltt,*,2),mdi,se,lc,uc))*10.	;decadal
      trendJJAU(lnn,ltt)=uc*10.
      trendJJAL(lnn,ltt)=lc*10.
    ENDIF
    gots=WHERE(qseas(lnn,ltt,*,3) NE mdi,countq)
    IF (countq GT nyrs*0.7) THEN BEGIN
      trendSON(lnn,ltt)=(median_pairwise(qseas(lnn,ltt,*,3),mdi,se,lc,uc))*10.	;decadal
      trendSONU(lnn,ltt)=uc*10.
      trendSONL(lnn,ltt)=lc*10.
    ENDIF
  ENDFOR
ENDFOR
  
;save to file;
;
filename=dir+'TRENDS/'+ofile+'.nc'
file_out=NCDF_CREATE(filename,/clobber)
latid=NCDF_DIMDEF(file_out,'latitude',nlats)
lonid=NCDF_DIMDEF(file_out,'longitude',nlons)
  
latsvar=NCDF_VARDEF(file_out,'latitude',[latid],/FLOAT)
lonsvar=NCDF_VARDEF(file_out,'longitude',[lonid],/FLOAT)
qDJFvar=NCDF_VARDEF(file_out,'DJFMPtrend',[lonid,latid],/FLOAT)
qDJFvarU=NCDF_VARDEF(file_out,'DJFMP95th',[lonid,latid],/FLOAT)
qDJFvarL=NCDF_VARDEF(file_out,'DJFMP5th',[lonid,latid],/FLOAT)
qMAMvar=NCDF_VARDEF(file_out,'MAMMPtrend',[lonid,latid],/FLOAT)
qMAMvarU=NCDF_VARDEF(file_out,'MAMMP95th',[lonid,latid],/FLOAT)
qMAMvarL=NCDF_VARDEF(file_out,'MAMMP5th',[lonid,latid],/FLOAT)
qJJAvar=NCDF_VARDEF(file_out,'JJAMPtrend',[lonid,latid],/FLOAT)
qJJAvarU=NCDF_VARDEF(file_out,'JJAMP95th',[lonid,latid],/FLOAT)
qJJAvarL=NCDF_VARDEF(file_out,'JJAMP5th',[lonid,latid],/FLOAT)
qSONvar=NCDF_VARDEF(file_out,'SONMPtrend',[lonid,latid],/FLOAT)
qSONvarU=NCDF_VARDEF(file_out,'SONMP95th',[lonid,latid],/FLOAT)
qSONvarL=NCDF_VARDEF(file_out,'SONMP5th',[lonid,latid],/FLOAT)

NCDF_ATTPUT,file_out,latsvar,'long_name','Latitude'
NCDF_ATTPUT,file_out,latsvar,'units','Degrees'
NCDF_ATTPUT,file_out,latsvar,'valid_min',-90.
NCDF_ATTPUT,file_out,latsvar,'valid_max',90.
NCDF_ATTPUT,file_out,lonsvar,'long_name','Longitude'
NCDF_ATTPUT,file_out,lonsvar,'units','Degrees'
NCDF_ATTPUT,file_out,lonsvar,'valid_min',-180.
NCDF_ATTPUT,file_out,lonsvar,'valid_max',180.

NCDF_ATTPUT,file_out,qDJFvar,'units',units
NCDF_ATTPUT,file_out,qDJFvar,'long_name','DJF Decadal MP trends'
NCDF_ATTPUT,file_out,qDJFvar,'missing_value',mdi

NCDF_ATTPUT,file_out,qDJFvarU,'units',units
NCDF_ATTPUT,file_out,qDJFvarU,'long_name','DJF Decadal MP 95th pct trends'
NCDF_ATTPUT,file_out,qDJFvarU,'missing_value',mdi

NCDF_ATTPUT,file_out,qDJFvarL,'units',units
NCDF_ATTPUT,file_out,qDJFvarL,'long_name','DJF Decadal MP 5th pct trends'
NCDF_ATTPUT,file_out,qDJFvarL,'missing_value',mdi

NCDF_ATTPUT,file_out,qMAMvar,'units',units
NCDF_ATTPUT,file_out,qMAMvar,'long_name','MAM Decadal MP trends'
NCDF_ATTPUT,file_out,qMAMvar,'missing_value',mdi

NCDF_ATTPUT,file_out,qMAMvarU,'units',units
NCDF_ATTPUT,file_out,qMAMvarU,'long_name','MAM Decadal MP 95th pct trends'
NCDF_ATTPUT,file_out,qMAMvarU,'missing_value',mdi

NCDF_ATTPUT,file_out,qMAMvarL,'units',units
NCDF_ATTPUT,file_out,qMAMvarL,'long_name','MAM Decadal MP 5th pct trends'
NCDF_ATTPUT,file_out,qMAMvarL,'missing_value',mdi

NCDF_ATTPUT,file_out,qJJAvar,'units',units
NCDF_ATTPUT,file_out,qJJAvar,'long_name','JJA Decadal MP trends'
NCDF_ATTPUT,file_out,qJJAvar,'missing_value',mdi

NCDF_ATTPUT,file_out,qJJAvarU,'units',units
NCDF_ATTPUT,file_out,qJJAvarU,'long_name','JJA Decadal MP 95th pct trends'
NCDF_ATTPUT,file_out,qJJAvarU,'missing_value',mdi

NCDF_ATTPUT,file_out,qJJAvarL,'units',units
NCDF_ATTPUT,file_out,qJJAvarL,'long_name','JJA Decadal MP 5th pct trends'
NCDF_ATTPUT,file_out,qJJAvarL,'missing_value',mdi

NCDF_ATTPUT,file_out,qSONvar,'units',units
NCDF_ATTPUT,file_out,qSONvar,'long_name','SON Decadal MP trends'
NCDF_ATTPUT,file_out,qSONvar,'missing_value',mdi

NCDF_ATTPUT,file_out,qSONvarU,'units',units
NCDF_ATTPUT,file_out,qSONvarU,'long_name','SON Decadal MP 95th pct trends'
NCDF_ATTPUT,file_out,qSONvarU,'missing_value',mdi

NCDF_ATTPUT,file_out,qSONvarL,'units',units
NCDF_ATTPUT,file_out,qSONvarL,'long_name','SON Decadal MP 5th pct trends'
NCDF_ATTPUT,file_out,qSONvarL,'missing_value',mdi

;NCDF_ATTPUT,file_out,qbootsvar,'units',units
;NCDF_ATTPUT,file_out,qbootsvar,'long_name','Bootstrap Decadal MP trends'
;NCDF_ATTPUT,file_out,qbootsvar,'missing_value',mdi
;
;NCDF_ATTPUT,file_out,qbootsvarU,'units',units
;NCDF_ATTPUT,file_out,qbootsvarU,'long_name','Bootstrap Decadal MP 95th pct trends'
;NCDF_ATTPUT,file_out,qbootsvarU,'missing_value',mdi
;
;NCDF_ATTPUT,file_out,qbootsvarL,'units',units
;NCDF_ATTPUT,file_out,qbootsvarL,'long_name','Bootstrap Decadal MP 5th pct trends'
;NCDF_ATTPUT,file_out,qbootsvarL,'missing_value',mdi
;
NCDF_CONTROL,file_out,/ENDEF
NCDF_VARPUT,file_out,qDJFvar,trendDJF
NCDF_VARPUT,file_out,qDJFvarU,trendDJFU
NCDF_VARPUT,file_out,qDJFvarL,trendDJFL
NCDF_VARPUT,file_out,qMAMvar,trendMAM
NCDF_VARPUT,file_out,qMAMvarU,trendMAMU
NCDF_VARPUT,file_out,qMAMvarL,trendMAML
NCDF_VARPUT,file_out,qJJAvar,trendJJA
NCDF_VARPUT,file_out,qJJAvarU,trendJJAU
NCDF_VARPUT,file_out,qJJAvarL,trendJJAL
NCDF_VARPUT,file_out,qSONvar,trendSON
NCDF_VARPUT,file_out,qSONvarU,trendSONU
NCDF_VARPUT,file_out,qSONvarL,trendSONL
;NCDF_VARPUT,file_out,qbootsvar,boottrendvals
;NCDF_VARPUT,file_out,qbootsvarU,boottrendvalsU
;NCDF_VARPUT,file_out,qbootsvarL,boottrendvalsL
;NCDF_VARPUT,file_out,bootsvar,int_boots
NCDF_CLOSE,file_out

return

end

