PRO plot_HadISDHdectrendsSEAS_FEB2015

version='2.0.1.2014p'
nowmon='JAN'
nowyear='2015'
thenmon='JAN'
thenyear='2015'
homogtype='IDPHA'	;'IDPHA, PHA, PHADPD
param2='T'
param='t'
styr='1973'
edyr='2014'

; files and directories
dir='/data/local/hadkw/HADCRUH2/UPDATE2014/'
filein='STATISTICS/TRENDS/HadISDH.land'+param2+'.'+version+'_FLATgrid'+homogtype+'5by5_'+thenmon+thenyear+'_MPtrendsSEAS_19732014.nc'
fileout='IMAGES/MAPS/HadISDH.land'+param2+'.'+version+'_FLATgrid'+homogtype+'5by5_'+nowmon+nowyear+'_MPtrendsSEAS_19732014.eps'

latlg=5.
lonlg=5.
stlt=-90+(latlg/2.)
stln=-180+(lonlg/2.)
nlats=180/latlg
nlons=360/lonlg
nbox=LONG(nlats*nlons)

lats=(findgen(nlats)*latlg)+stlt
lons=(findgen(nlons)*lonlg)+stln

mdi=-1e+30

;--------------------------------------

filee=NCDF_OPEN(dir+filein)
longs_varid=NCDF_VARID(filee,'longitude')
lats_varid=NCDF_VARID(filee,'latitude')
qDJFid=NCDF_VARID(filee,'DJFMPtrend') 	; may become uncertainty fields
qDJF95id=NCDF_VARID(filee,'DJFMP95th') 	; may become uncertainty fields
qDJF5id=NCDF_VARID(filee,'DJFMP5th') 	; may become uncertainty fields
qMAMid=NCDF_VARID(filee,'MAMMPtrend') 	; may become uncertainty fields
qMAM95id=NCDF_VARID(filee,'MAMMP95th') 	; may become uncertainty fields
qMAM5id=NCDF_VARID(filee,'MAMMP5th') 	; may become uncertainty fields
qJJAid=NCDF_VARID(filee,'JJAMPtrend') 	; may become uncertainty fields
qJJA95id=NCDF_VARID(filee,'JJAMP95th') 	; may become uncertainty fields
qJJA5id=NCDF_VARID(filee,'JJAMP5th') 	; may become uncertainty fields
qSONid=NCDF_VARID(filee,'SONMPtrend') 	; may become uncertainty fields
qSON95id=NCDF_VARID(filee,'SONMP95th') 	; may become uncertainty fields
qSON5id=NCDF_VARID(filee,'SONMP5th') 	; may become uncertainty fields
NCDF_VARGET,filee,qDJFid,qDJF_trend
NCDF_VARGET,filee,qDJF95id,qDJF_95th
NCDF_VARGET,filee,qDJF5id,qDJF_5th
NCDF_VARGET,filee,qMAMid,qMAM_trend
NCDF_VARGET,filee,qMAM95id,qMAM_95th
NCDF_VARGET,filee,qMAM5id,qMAM_5th
NCDF_VARGET,filee,qJJAid,qJJA_trend
NCDF_VARGET,filee,qJJA95id,qJJA_95th
NCDF_VARGET,filee,qJJA5id,qJJA_5th
NCDF_VARGET,filee,qSONid,qSON_trend
NCDF_VARGET,filee,qSON95id,qSON_95th
NCDF_VARGET,filee,qSON5id,qSON_5th
NCDF_CLOSE,filee
    
CASE param OF
  'rh': units='%rh'
  'e': units='hPa'
  'q': units='g kg!E-1!N'
  ELSE: units='!Eo!NC'
ENDCASE

CASE param OF
  'dpd': BEGIN
    kcolsarrST=[100,1,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.,-0.8,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,10]
  END
  'td': BEGIN
    kcolsarrST=[100,1,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.,-0.8,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,10]
  END
  't': BEGIN
    kcolsarrST=[100,1,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.,-0.8,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,10]
  END
  'tw': BEGIN
    kcolsarrST=[100,1,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.,-0.8,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,10]
  END
  'q': BEGIN
    kcolsarrST=[100,1,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05, 0.0 ,0.05,0.1,0.15,0.2,0.25,0.3,0.35,10]
  END
  'e': BEGIN
    kcolsarrST=[100,1,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.,-0.5,-0.4,-0.3,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,10]
  END
  'rh': BEGIN
    kcolsarrST=[100,1,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10,-3,-2,-1.5,-1.,-0.75,-0.5,-0.25,  0.0  ,0.25,0.5,0.75,1.,1.5,2.,3.,10]
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
labsarrST(1:nlevs-2)=strcompress(string(levsarrST(2:nlevs-1),format='(f5.2)'),/remove_all)
;------------------------------------------------------------------------------------
xpos1=[0.05,0.55,0.05,0.55]
xpos2=[0.45,0.95,0.45,0.94]
ypos1=[0.60,0.6,0.15,0.15]
ypos2=[0.95,0.95,0.5,0.5]

;make the plot
set_plot,'PS'
device,filename=dir+fileout,/ENCAPSUL,$
       xsize=26,ysize=16,/landscape,/color,/helvetica,/bold

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

titlees=['DJF','MAM','JJA','SON']

!P.Position=[xpos1(0),ypos1(0),xpos2(0),ypos2(0)]

Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1
XYOUTS,xpos1(0)+((xpos2(0)-xpos1(0))/2.),ypos2(0)+0.01,titlees(0),/NORMAL,alignment=0.5,charsize=1.2,color=0
XYOUTS,xpos1(0)-0.01,ypos2(0)-0.02,'a)',/NORMAL,alignment=0.5,charsize=1.2,color=0


boxfill,qDJF_trend,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1

FOR lnn=0,nlons-1 DO BEGIN
  FOR ltt=0,nlats-1 DO BEGIN
    IF (qDJF_trend(lnn,ltt) NE mdi) THEN BEGIN
      print,qDJF_trend(lnn,ltt),qDJF_95th(lnn,ltt),qDJF_5th(lnn,ltt)
      IF ((qDJF_trend(lnn,ltt) GT 0.) AND (qDJF_5th(lnn,ltt) GT 0.)) OR $
         ((qDJF_trend(lnn,ltt) LT 0.) AND (qDJF_95th(lnn,ltt) LT 0.)) THEN BEGIN
         PLOTS,lons(lnn),lats(ltt),psym=8,color=0,symsize=0.5
         print,'SIGNIFICANT!'
      ENDIF
    ENDIF
  ENDFOR
ENDFOR

PLOTS,[-179.9,0],[89.9,89.9],color=0
PLOTS,[0,179.9],[89.9,89.9],color=0
PLOTS,[-179.9,0],[-89.9,-89.9],color=0
PLOTS,[0,179.9],[-89.9,-89.9],color=0

!P.Position=[xpos1(1),ypos1(1),xpos2(1),ypos2(1)]

Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1
XYOUTS,xpos1(1)+((xpos2(1)-xpos1(1))/2.),ypos2(1)+0.01,titlees(1),/NORMAL,alignment=0.5,charsize=1.2,color=0
XYOUTS,xpos1(1)-0.01,ypos2(1)-0.02,'b)',/NORMAL,alignment=0.5,charsize=1.2,color=0


boxfill,qMAM_trend,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1

FOR lnn=0,nlons-1 DO BEGIN
  FOR ltt=0,nlats-1 DO BEGIN
    IF (qMAM_trend(lnn,ltt) NE mdi) THEN BEGIN
      print,qMAM_trend(lnn,ltt),qMAM_95th(lnn,ltt),qMAM_5th(lnn,ltt)
      IF ((qMAM_trend(lnn,ltt) GT 0.) AND (qMAM_5th(lnn,ltt) GT 0.)) OR $
         ((qMAM_trend(lnn,ltt) LT 0.) AND (qMAM_95th(lnn,ltt) LT 0.)) THEN BEGIN
         PLOTS,lons(lnn),lats(ltt),psym=8,color=0,symsize=0.5
         print,'SIGNIFICANT!'
      ENDIF
    ENDIF
  ENDFOR
ENDFOR

PLOTS,[-179.9,0],[89.9,89.9],color=0
PLOTS,[0,179.9],[89.9,89.9],color=0
PLOTS,[-179.9,0],[-89.9,-89.9],color=0
PLOTS,[0,179.9],[-89.9,-89.9],color=0

!P.Position=[xpos1(2),ypos1(2),xpos2(2),ypos2(2)]

Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1
XYOUTS,xpos1(2)+((xpos2(2)-xpos1(2))/2.),ypos2(2)+0.01,titlees(2),/NORMAL,alignment=0.5,charsize=1.2,color=0
XYOUTS,xpos1(2)-0.01,ypos2(2)-0.02,'c)',/NORMAL,alignment=0.5,charsize=1.2,color=0


boxfill,qJJA_trend,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1

FOR lnn=0,nlons-1 DO BEGIN
  FOR ltt=0,nlats-1 DO BEGIN
    IF (qJJA_trend(lnn,ltt) NE mdi) THEN BEGIN
      print,qJJA_trend(lnn,ltt),qJJA_95th(lnn,ltt),qJJA_5th(lnn,ltt)
      IF ((qJJA_trend(lnn,ltt) GT 0.) AND (qJJA_5th(lnn,ltt) GT 0.)) OR $
         ((qJJA_trend(lnn,ltt) LT 0.) AND (qJJA_95th(lnn,ltt) LT 0.)) THEN BEGIN
         PLOTS,lons(lnn),lats(ltt),psym=8,color=0,symsize=0.5
         print,'SIGNIFICANT!'
      ENDIF
    ENDIF
  ENDFOR
ENDFOR

PLOTS,[-179.9,0],[89.9,89.9],color=0
PLOTS,[0,179.9],[89.9,89.9],color=0
PLOTS,[-179.9,0],[-89.9,-89.9],color=0
PLOTS,[0,179.9],[-89.9,-89.9],color=0

!P.Position=[xpos1(3),ypos1(3),xpos2(3),ypos2(3)]

Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1
XYOUTS,xpos1(3)+((xpos2(3)-xpos1(3))/2.),ypos2(3)+0.01,titlees(3),/NORMAL,alignment=0.5,charsize=1.2,color=0
XYOUTS,xpos1(3)-0.01,ypos2(3)-0.02,'d)',/NORMAL,alignment=0.5,charsize=1.2,color=0


boxfill,qSON_trend,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1

FOR lnn=0,nlons-1 DO BEGIN
  FOR ltt=0,nlats-1 DO BEGIN
    IF (qSON_trend(lnn,ltt) NE mdi) THEN BEGIN
      print,qSON_trend(lnn,ltt),qSON_95th(lnn,ltt),qSON_5th(lnn,ltt)
      IF ((qSON_trend(lnn,ltt) GT 0.) AND (qSON_5th(lnn,ltt) GT 0.)) OR $
         ((qSON_trend(lnn,ltt) LT 0.) AND (qSON_95th(lnn,ltt) LT 0.)) THEN BEGIN
         PLOTS,lons(lnn),lats(ltt),psym=8,color=0,symsize=0.5
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
         charsize=1.2,charthick=4,bcolor=0  ;,orientation=1

XYOUTS,0.5,0.01,styr+' to '+edyr+' trends ('+units+' decade!E-1!N)',/normal,color=0,charsize=1.4,alignment=0.5
DEVICE,/close



return

END
