PRO plot_HadISDH_dectrends_JAN2015

param='t'	;'dpd','td','t','tw','e','q','rh'
param2='T'	;'DPD','Td','T','Tw','e','q','RH'
nowmon='MAR'
nowyear='2015'
homogtype='OTHER'	;'PHA','ID','DPD', 'RAW'
version='2.0.1.2014p'
styr='1973'
edyr='2014'

IF (homogtype EQ 'OTHER') $
    THEN indir='/data/local/hadkw/HADCRUH2/UPDATE2014/OTHERDATA/' $
    ELSE indir='/data/local/hadkw/HADCRUH2/UPDATE2014/STATISTICS/TRENDS/'
odir='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/MAPS/'

CASE param OF
  'dpd': BEGIN
    IF (homogtype EQ 'PHA') THEN BEGIN
      infile='HadISDH.landDPD.'+version+'_FLATgridPHA5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landDPD.'+version+'_FLATgridPHA5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF ELSE IF (homogtype EQ 'RAW') THEN BEGIN
      infile='HadISDH.landDPD.'+version+'_FLATgridRAW5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landDPD.'+version+'_FLATgridRAW5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF
    variname='dew point depression'
    letteree='c)'
  END
  'td': BEGIN
    IF (homogtype EQ 'PHA') THEN BEGIN
      infile='HadISDH.landTd.'+version+'_FLATgridPHA5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landTd.'+version+'_FLATgridPHA5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF ELSE IF (homogtype EQ 'DPD') THEN BEGIN
      infile='HadISDH.landTd.'+version+'_FLATgridPHADPD5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landTd.'+version+'_FLATgridPHADPD5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF ELSE IF (homogtype EQ 'RAW') THEN BEGIN
      infile='HadISDH.landTd.'+version+'_FLATgridRAW5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landTd.'+version+'_FLATgridRAW5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF
    variname='dew point temperature'
    letteree='b)'
  END
  't': BEGIN
    IF (homogtype EQ 'ID') THEN BEGIN
      infile='HadISDH.landT.'+version+'_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landT.'+version+'_FLATgridIDPHA5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF ELSE IF (homogtype EQ 'RAW') THEN BEGIN
      infile='HadISDH.landT.'+version+'_FLATgridRAW5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landT.'+version+'_FLATgridRAW5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF ELSE IF (homogtype EQ 'OTHER') THEN BEGIN
      infile='CRUTEM.4.3.0.0.anomalies_MPtrends_19732014'
      outfile='CRUTEM.4.3.0.0.anomalies_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
;      infile='HadCRUT.4.3.0.0.median_MPtrends_19732014'
;      outfile='HadCRUT.4.3.0.0.median_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
;      infile='GHCNM_18802014_MPtrends_19732014'
;      outfile='GHCNM_18802014_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF
    variname='temperature'
    letteree='a)'
  END
  'tw': BEGIN
    IF (homogtype EQ 'ID') THEN BEGIN
      infile='HadISDH.landTw.'+version+'_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landTw.'+version+'_FLATgridIDPHA5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF ELSE IF (homogtype EQ 'RAW') THEN BEGIN
      infile='HadISDH.landTw.'+version+'_FLATgridRAW5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landTw.'+version+'_FLATgridRAW5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF
    variname='wet bulb temperature'
    letteree='a)'
  END
  'q': BEGIN
    IF (homogtype EQ 'ID') THEN BEGIN
      infile='HadISDH.landq.'+version+'_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landq.'+version+'_FLATgridIDPHA5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF ELSE IF (homogtype EQ 'RAW') THEN BEGIN
      infile='HadISDH.landq.'+version+'_FLATgridRAW5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landq.'+version+'_FLATgridRAW5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF ELSE IF (homogtype EQ 'PHA') THEN BEGIN
      infile='HadISDH.landq.'+version+'_FLATgridPHA5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landq.'+version+'_FLATgridPHA5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF
    variname='specific humidity'
    letteree='c)'
  END
  'e': BEGIN
    IF (homogtype EQ 'ID') THEN BEGIN
      infile='HadISDH.lande.'+version+'_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.lande.'+version+'_FLATgridIDPHA5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF ELSE IF (homogtype EQ 'RAW') THEN BEGIN
      infile='HadISDH.lande.'+version+'_FLATgridRAW5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.lande.'+version+'_FLATgridRAW5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF
    variname='vapour pressure'
    letteree='b)'
  END
  'rh': BEGIN
    IF (homogtype EQ 'ID') THEN BEGIN
      infile='HadISDH.landRH.'+version+'_FLATgridIDPHA5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landRH.'+version+'_FLATgridIDPHA5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF ELSE IF (homogtype EQ 'RAW') THEN BEGIN
      infile='HadISDH.landRH.'+version+'_FLATgridRAW5by5_JAN2015_MPtrends_19732014'
      outfile='HadISDH.landRH.'+version+'_FLATgridRAW5by5_'+nowmon+nowyear+'_MPdectrends_19732014.eps'
    ENDIF
    variname='relative humidity'
    letteree='d)'
  END
  
ENDCASE


;--------------------------------------------------------
; variables and arrays
mdi=-1e+30

CASE param OF
  'rh': unitees='%rh'
  'e': unitees='hPa'
  'q': unitees='g kg!E-1!N'
  ELSE: unitees='!Eo!NC'
ENDCASE

latlg=5.
lonlg=5.
stlt=-90+(latlg/2.)
stln=-180+(lonlg/2.)
nlats=180/latlg
nlons=360/lonlg
nbox=LONG(nlats*nlons)

lats=(findgen(nlats)*latlg)+stlt
lons=(findgen(nlons)*lonlg)+stln

;--------------------------------------

filee=NCDF_OPEN(indir+infile+'.nc')
longs_varid=NCDF_VARID(filee,'longitude')
lats_varid=NCDF_VARID(filee,'latitude')

CASE param OF
  'e': BEGIN
    qid=NCDF_VARID(filee,'e_MPtrend') 	; may become uncertainty fields
    q95id=NCDF_VARID(filee,'e_MP95th') 	; may become uncertainty fields
    q5id=NCDF_VARID(filee,'e_MP5th') 	; may become uncertainty fields
  END
  'dpd': BEGIN
    qid=NCDF_VARID(filee,'DPD_MPtrend') 	; may become uncertainty fields
    q95id=NCDF_VARID(filee,'DPD_MP95th') 	; may become uncertainty fields
    q5id=NCDF_VARID(filee,'DPD_MP5th') 	; may become uncertainty fields
  END
  'td': BEGIN
    qid=NCDF_VARID(filee,'Td_MPtrend') 	; may become uncertainty fields
    q95id=NCDF_VARID(filee,'Td_MP95th') 	; may become uncertainty fields
    q5id=NCDF_VARID(filee,'Td_MP5th') 	; may become uncertainty fields
  END
  't': BEGIN
    qid=NCDF_VARID(filee,'T_MPtrend') 	; may become uncertainty fields
    q95id=NCDF_VARID(filee,'T_MP95th') 	; may become uncertainty fields
    q5id=NCDF_VARID(filee,'T_MP5th') 	; may become uncertainty fields
  END
  'tw': BEGIN
    qid=NCDF_VARID(filee,'Tw_MPtrend') 	; may become uncertainty fields
    q95id=NCDF_VARID(filee,'Tw_MP95th') 	; may become uncertainty fields
    q5id=NCDF_VARID(filee,'Tw_MP5th') 	; may become uncertainty fields
  END
  'q': BEGIN
    qid=NCDF_VARID(filee,'q_MPtrend') 	; may become uncertainty fields
    q95id=NCDF_VARID(filee,'q_MP95th') 	; may become uncertainty fields
    q5id=NCDF_VARID(filee,'q_MP5th') 	; may become uncertainty fields
  END
  'rh': BEGIN
    qid=NCDF_VARID(filee,'RH_MPtrend') 	; may become uncertainty fields
    q95id=NCDF_VARID(filee,'RH_MP95th') 	; may become uncertainty fields
    q5id=NCDF_VARID(filee,'RH_MP5th') 	; may become uncertainty fields  
  END
ENDCASE

;qbid=NCDF_VARID(filee,'qhum_MPboottrend') 	; may become uncertainty fields
;qb95id=NCDF_VARID(filee,'qhum_MPboot95th') 	; may become uncertainty fields
;qb5id=NCDF_VARID(filee,'qhum_MPboot5th') 	; may become uncertainty fields
NCDF_VARGET,filee,qid,q_trend
NCDF_VARGET,filee,q95id,q_95th
NCDF_VARGET,filee,q5id,q_5th
;NCDF_VARGET,filee,qbid,q_Btrend
;NCDF_VARGET,filee,qb95id,q_B95th
;NCDF_VARGET,filee,qb5id,q_B5th
NCDF_CLOSE,filee
 
CASE param OF
  'dpd': BEGIN
    kcolsarrST=[255,1,3,4,5,6,8,9,10,11,13,14,15,16,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,10]
  END
  'td': BEGIN
    kcolsarrST=[255,1,4,5,6,7,8,9,10,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,10]
  END
  't': BEGIN
    kcolsarrST=[100,1,3,4,5,6,7,9,10,11,12,13,14,17,18]	; 14 colours + black to ensure no botching
;    kcolsarrST=[255,1,3,4,5,6,7,9,10,11,12,13,14,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,10]
  END
  'tw': BEGIN
    kcolsarrST=[255,1,4,5,6,7,8,9,10,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,10]
  END
  'q': BEGIN
    kcolsarrST=[100,1,4,5,6,7,8,9,10,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
;    kcolsarrST=[255,1,4,5,6,7,8,9,10,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-0.4,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,10]
  END
  'e': BEGIN
    kcolsarrST=[255,1,4,5,6,7,8,9,10,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.,-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0.0,0.05,0.1,0.2,0.3,0.4,0.5,10]
  END
  'rh': BEGIN
;    kcolsarrST=[255,1,4,5,6,7,8,9,10,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    kcolsarrST=[100,1,4,5,6,7,8,9,10,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10,-2,-1.5,-1.,-0.75,-0.5,-0.25,  0.0  ,0.25,0.5,0.75,1.,1.5,2.,10]
  END
ENDCASE  
 
; colour settings - Uni ORegon blue to red
tvlct,0,0,0,0
tvlct,240,240,240,100

ncols=n_elements(kcolsarrST)
colsarrST=kcolsarrST(1:ncols-1)
nlevs=n_elements(levsarrST)-1
labsarrST=strarr(nlevs)
labsarrST(0)=''
labsarrST((nlevs-1))=''
labsarrST(1:nlevs-2)=strcompress(string(levsarrST(2:nlevs-1),format='(f6.2)'),/remove_all)
;------------------------------------------------------------------------------------
;make the plot
set_plot,'PS'
device,filename=odir+outfile,/ENCAPSUL,$
       xsize=26,ysize=20,/landscape,/color,/helvetica,/bold
;       xsize=26,ysize=20,/landscape,/color,/helvetica,/bold

plotsym,0,0.55,/fill
!P.Font=0
!P.Thick=4

CASE param OF
  'dpd': BEGIN
; blue to brown
    tvlct,0,10,40,1
    tvlct,0,40,80,2
    tvlct,0,75,100,3
    tvlct,0,122,153,4
    tvlct,0,169,204,5
    tvlct,50,227,255,6
    tvlct,101,239,255,7
    tvlct,153,248,255,8
    tvlct,204,253,255,9
    tvlct,242,218,205,10
    tvlct,216,175,151,11
    tvlct,204,155,122,12
    tvlct,153,96,53,13
    tvlct,102,47,0,14
    tvlct,51,25,0,15
    tvlct,30,5,0,16
    tvlct,20,1,0,17
    tvlct,10,0,0,18
  END
  't': BEGIN
;blue to red
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
  END
  ELSE: BEGIN
; brown to blue
    tvlct,10,0,0,1
    tvlct,20,1,0,2
    tvlct,30,5,0,3
    tvlct,51,25,0,4
    tvlct,102,47,0,5
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
  END
ENDCASE

  tvlct,230,230,230,100
  
;!P.Position=[0.1,0.15,0.9,0.95]
!P.Position=[0.02,0.01,0.99,0.95]

Map_set,/continents,color=0,/noborder,/noerase,/robinson,/isotropic,/horizon,/grid,glinestyle=1; /grid, glinestyle=2

;; *** WATERMARK
;XYOUTS,0.5,0.5,'DRAFT',/normal,color=100,charsize=16,alignment=0.5,orientation=30


;!P.Position=[0.1,0.15,0.9,0.95]
!P.Position=[0.02,0.01,0.99,0.95]
boxfill,q_trend,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
;!P.Position=[0.1,0.15,0.9,0.95]
!P.Position=[0.02,0.01,0.99,0.95]
Map_set,/continents,color=0,/noborder,/noerase,/robinson,/isotropic,/horizon,/grid,glinestyle=1; ,/grid, glinestyle=2

FOR lnn=0,nlons-1 DO BEGIN
  FOR ltt=0,nlats-1 DO BEGIN
    IF (q_trend(lnn,ltt) NE mdi) THEN BEGIN
      print,q_trend(lnn,ltt),q_95th(lnn,ltt),q_5th(lnn,ltt)
      IF ((q_trend(lnn,ltt) GT 0.) AND (q_5th(lnn,ltt) GT 0.)) OR $
         ((q_trend(lnn,ltt) LT 0.) AND (q_95th(lnn,ltt) LT 0.)) THEN BEGIN
         PLOTS,lons(lnn),lats(ltt),psym=8,color=0,symsize=1.
         print,'SIGNIFICANT!'
      ENDIF
    ENDIF
  ENDFOR
ENDFOR

PLOTS,[-179.9,0],[89.9,89.9],color=0
PLOTS,[0,179.9],[89.9,89.9],color=0
PLOTS,[-179.9,0],[-89.9,-89.9],color=0
PLOTS,[0,179.9],[-89.9,-89.9],color=0

MAKE_KEY,0.05,0.08,0.9,0.03,0.0,-0.03,/NORMAL,COLORS=colsarrST,labels=labsarrST,$
         charsize=1.3,charthick=4,bcolor=0  ;,orientation=1

;XYOUTS,0.5,0.01,styr+' to '+edyr+' trends ('+unitees+' decade!E-1!N)',/normal,color=0,charsize=1.4,alignment=0.5
XYOUTS,0.5,0.01,variname+' trends ('+unitees+' decade!E-1!N)',/normal,color=0,charsize=1.5,alignment=0.5
;XYOUTS,0.05,0.9,letteree,/normal,color=0,charsize=2.5,alignment=0.5

;      MAKE_KEY,0.05,keypos,0.9,0.02,charalign,charmove,/NORMAL,COLORS=colsarr,labels=labsarr,$
;         charsize=keythick,charthick=4,bcolor=0
;      XYOUTS,0.5,0.02,titlees(1),/normal,alignment=0.5,charsize=1.2,color=0


DEVICE,/close


return

END
