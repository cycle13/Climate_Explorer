pro plot_seasonalanommaps_FEB2014,year,filein,fileout,plotout,climst,climed,param


CASE param OF
  'q': BEGIN
    titlees='Surface Specific Humidity'
    units='(g kg!E-1!N)'    
    unitees='g/kg'
    kcolsarrST=[100,2,4,5,6,7,9,10,12,14,15,16,17]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.0,-2.0,-1.0,-0.5,-0.25,-0.1,  0.0,  0.1,0.25,0.5,1.0,2.0,10]
  END
  'RH': BEGIN
    titlees='Surface Relative Humidity'
    units='(%rh)'   
    unitees='%rh' 
    kcolsarrST=[100,2,4,5,6,7,9,10,12,14,15,16,17]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-100.0,-10.0,-5.0,-2.5,-1.0,-0.5,  0.0, 0.5,1.0,2.5,5.,10.0,100]
  END
  'e': BEGIN
    titlees='Surface Vapour Pressure'
    units='(hPa)'  
    unitees='hPa'  
    kcolsarrST=[100,2,3,4,5,6,7,9,10,12,13,14,15,16,17]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.0,-2.0,-1.0,-0.75,-0.5,-0.25,-0.1,  0.0,  0.1,0.25,0.5,0.75,1.0,2.0,10]
  END
  'Td': BEGIN
    titlees='Surface Dew Point Temperature'
    units='(!Eo!NC)'
    unitees='deg C'    
    kcolsarrST=[100,2,3,4,5,6,7,9,10,12,13,14,15,16,17]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.0,-2.5,-2.0,-1.5,-1.0,-0.5,-0.25, 0.0, 0.25,0.5,1.0,1.5,2.0,2.5,10]
  END
  'Tw': BEGIN
    titlees='Surface Wet Bulb Temperature'
    units='(!Eo!NC)'
    unitees='deg C'    
;    kcolsarrST=[100,2,3,4,5,6,7,9,10,12,13,14,15,16,17]	; 14 colours + black to ensure no botching
;    levsarrST=[-2e+30,-10.0,-2.0,-1.5,-1.0,-0.75,-0.5,-0.25,  0.0,  0.25,0.5,0.75,1.0,1.5,2.0,10]
    kcolsarrST=[100,2,3,4,5,6,7,9,10,12,13,14,15,16,17]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.0,-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,  0.0,  0.25,0.5,0.75,1.0,1.25,1.5,10]
  END
  'T': BEGIN
    titlees='Surface Temperature'
    units='(!Eo!NC)'    
    unitees='deg C'    
    kcolsarrST=[100,2,3,4,5,6,7,9,10,12,13,14,15,16,17]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.0,-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,  0.0,  0.25,0.5,0.75,1.0,1.25,1.5,10]
  END
  'DPD': BEGIN
    titlees='Surface Dew Point Depression'
    units='(!Eo!NC)'    
    unitees='deg C'    
    kcolsarrST=[100,2,3,4,5,6,7,9,10,12,13,14,15,16,17]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.25,  0.0,  0.25,0.5,1.0,2.0,3.0,4.0,10]
  END

ENDCASE

mdi=-1e30
edyr=2017
styr=1973
nyrs=(edyr+1)-styr
nmons=nyrs*12
nlons=72
nlats=36
maparrDJF=make_array(nlons,nlats,/float,value=mdi)
maparrMAM=make_array(nlons,nlats,/float,value=mdi)
maparrJJA=make_array(nlons,nlats,/float,value=mdi)
maparrSON=make_array(nlons,nlats,/float,value=mdi)

latlg=5.
lonlg=5.
stlt=-90+(latlg/2.)
stln=-180+(lonlg/2.)
nlats=180/latlg
nlons=360/lonlg
nbox=LONG(nlats*nlons)

lats=(findgen(nlats)*latlg)+stlt
lons=(findgen(nlons)*lonlg)+stln

getyrmonst=(year-styr)*12
getyrmoned=(((year+1)-styr)*12)-1

djf=[0,1]
mam=[2,3,4]
jja=[5,6,7]
son=[8,9,10]

inn=NCDF_OPEN(filein)
CASE param OF
  'q': varid=NCDF_VARID(inn,'q_anoms')
  'q': varid=NCDF_VARID(inn,'blendmask_q_anoms')
;  'q': varid=NCDF_VARID(inn,'blend_q_anoms')
  'RH': varid=NCDF_VARID(inn,'rh_anoms')
  'e': varid=NCDF_VARID(inn,'e_anoms')
  'Td': varid=NCDF_VARID(inn,'td_anoms')
  'Tw': varid=NCDF_VARID(inn,'tw_anoms')
  'T': varid=NCDF_VARID(inn,'t_anoms')
  'DPD': varid=NCDF_VARID(inn,'dpd_anoms')
ENDCASE
latid=NCDF_VARID(inn,'latitude')
lonid=NCDF_VARID(inn,'longitude')
NCDF_VARGET,inn,varid,tmp
NCDF_VARGET,inn,latid,lats
NCDF_VARGET,inn,lonid,lons
NCDF_CLOSE,inn

IF (climst NE 1976) THEN BEGIN	; rezero over new climatology period
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

;pull out the 12 months
IF (year NE styr) THEN subarr=tmp(*,*,getyrmonst-1:getyrmoned) ELSE subarr=tmp(*,*,getyrmonst:getyrmoned)


; average over the year for each gridbox

FOR ltt=0,nlats-1 DO BEGIN
  FOR lnn=0,nlons-1 DO BEGIN
    monalls=subarr(lnn,ltt,*)
    IF (year NE styr) THEN BEGIN
      tmp=monalls(mam+1)
      gots=WHERE(tmp NE mdi,count)
      IF (count GE 2) THEN maparrMAM(lnn,ltt)=MEAN(tmp(gots))
      tmp=monalls(jja+1)
      gots=WHERE(tmp NE mdi,count)
      IF (count GE 2) THEN maparrJJA(lnn,ltt)=MEAN(tmp(gots))
      tmp=monalls(son+1)
      gots=WHERE(tmp NE mdi,count)
      IF (count GE 2) THEN maparrSON(lnn,ltt)=MEAN(tmp(gots))
      tmp=monalls([0,djf+1])
      gots=WHERE(tmp NE mdi,count)
      IF (count GE 2) THEN maparrDJF(lnn,ltt)=MEAN(tmp(gots))
    ENDIF ELSE BEGIN
      tmp=monalls(mam)
      gots=WHERE(tmp NE mdi,count)
      IF (count GE 2) THEN maparrMAM(lnn,ltt)=MEAN(tmp(gots))
      tmp=monalls(jja)
      gots=WHERE(tmp NE mdi,count)
      IF (count GE 2) THEN maparrJJA(lnn,ltt)=MEAN(tmp(gots))
      tmp=monalls(son)
      gots=WHERE(tmp NE mdi,count)
      IF (count GE 2) THEN maparrSON(lnn,ltt)=MEAN(tmp(gots))
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
labsarrST(1:nlevs-2)=strcompress(string(levsarrST(2:nlevs-1),format='(f6.2)'),/remove_all)


;------------------------------------------------------------------------------------
xpos1=[0.05,0.55,0.05,0.55]
xpos2=[0.45,0.95,0.45,0.94]
ypos1=[0.60,0.6,0.15,0.15]
ypos2=[0.95,0.95,0.5,0.5]

;make the plot
set_plot,'PS'
device,filename=plotout+'.eps',/ENCAPSUL,$
       xsize=26,ysize=16,/landscape,/color,/helvetica,/bold


!P.Font=0
!P.Thick=4
!X.Thick=4
!Y.Thick=4

IF (param EQ 'q') OR (param EQ 'RH') OR (param EQ 'e') OR (param EQ 'Td') OR (param EQ 'Tw') THEN BEGIN
;  tvlct,10,0,0,1
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
;  tvlct,0,10,40,18
ENDIF ELSE IF (param EQ 'DPD') THEN BEGIN
;  tvlct,10,0,0,18
  tvlct,20,1,0,17
  tvlct,30,5,0,16
  tvlct,51,25,0,15
  tvlct,102,47,0,14
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
;  tvlct,0,10,40,1
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

titlees=['DJF','MAM','JJA','SON']

!P.Position=[xpos1(0),ypos1(0),xpos2(0),ypos2(0)]

Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1
XYOUTS,xpos1(0)+((xpos2(0)-xpos1(0))/2.),ypos2(0)+0.01,titlees(0),/NORMAL,alignment=0.5,charsize=1.2,color=0
XYOUTS,xpos1(0)-0.01,ypos2(0)-0.02,'a)',/NORMAL,alignment=0.5,charsize=1.2,color=0


boxfill,maparrDJF,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1

PLOTS,[-179.9,0],[89.9,89.9],color=0
PLOTS,[0,179.9],[89.9,89.9],color=0
PLOTS,[-179.9,0],[-89.9,-89.9],color=0
PLOTS,[0,179.9],[-89.9,-89.9],color=0

!P.Position=[xpos1(1),ypos1(1),xpos2(1),ypos2(1)]

Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1
XYOUTS,xpos1(1)+((xpos2(1)-xpos1(1))/2.),ypos2(1)+0.01,titlees(1),/NORMAL,alignment=0.5,charsize=1.2,color=0
XYOUTS,xpos1(1)-0.01,ypos2(1)-0.02,'b)',/NORMAL,alignment=0.5,charsize=1.2,color=0


boxfill,maparrMAM,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1

PLOTS,[-179.9,0],[89.9,89.9],color=0
PLOTS,[0,179.9],[89.9,89.9],color=0
PLOTS,[-179.9,0],[-89.9,-89.9],color=0
PLOTS,[0,179.9],[-89.9,-89.9],color=0

!P.Position=[xpos1(2),ypos1(2),xpos2(2),ypos2(2)]

Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1
XYOUTS,xpos1(2)+((xpos2(2)-xpos1(2))/2.),ypos2(2)+0.01,titlees(2),/NORMAL,alignment=0.5,charsize=1.2,color=0
XYOUTS,xpos1(2)-0.01,ypos2(2)-0.02,'c)',/NORMAL,alignment=0.5,charsize=1.2,color=0

boxfill,maparrJJA,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1

PLOTS,[-179.9,0],[89.9,89.9],color=0
PLOTS,[0,179.9],[89.9,89.9],color=0
PLOTS,[-179.9,0],[-89.9,-89.9],color=0
PLOTS,[0,179.9],[-89.9,-89.9],color=0

!P.Position=[xpos1(3),ypos1(3),xpos2(3),ypos2(3)]

Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1
XYOUTS,xpos1(3)+((xpos2(3)-xpos1(3))/2.),ypos2(3)+0.01,titlees(3),/NORMAL,alignment=0.5,charsize=1.2,color=0
XYOUTS,xpos1(3)-0.01,ypos2(3)-0.02,'d)',/NORMAL,alignment=0.5,charsize=1.2,color=0

boxfill,maparrSON,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)
     
; Overplot continents and islands again
Map_set,/continents,color=0,/noborder,/noerase,/robinson,/horizon,/isotropic,/grid,glinestyle=1

PLOTS,[-179.9,0],[89.9,89.9],color=0
PLOTS,[0,179.9],[89.9,89.9],color=0
PLOTS,[-179.9,0],[-89.9,-89.9],color=0
PLOTS,[0,179.9],[-89.9,-89.9],color=0

MAKE_KEY,0.05,0.08,0.9,0.03,0.0,-0.03,/NORMAL,COLORS=colsarrST,labels=labsarrST,$
         charsize=1.2,charthick=4,bcolor=0  ;,orientation=1

XYOUTS,0.5,0.01,'Annual Anomaly relative to '+strcompress(climst,/remove_all)+'-'$
                +strcompress(climed,/remove_all)+' '+units,/normal,color=0,$
		charsize=1.2,alignment=0.5

DEVICE,/close

;save to file;
;
filename=fileout+'.nc'
file_out=NCDF_CREATE(filename,/clobber)
latid=NCDF_DIMDEF(file_out,'latitude',nlats)
lonid=NCDF_DIMDEF(file_out,'longitude',nlons)
  
latsvar=NCDF_VARDEF(file_out,'lat',[latid],/FLOAT)
lonsvar=NCDF_VARDEF(file_out,'lon',[lonid],/FLOAT)
qvarDJF=NCDF_VARDEF(file_out,'anomaliesDJF',[lonid,latid],/FLOAT)
qvarMAM=NCDF_VARDEF(file_out,'anomaliesMAM',[lonid,latid],/FLOAT)
qvarJJA=NCDF_VARDEF(file_out,'anomaliesJJA',[lonid,latid],/FLOAT)
qvarSON=NCDF_VARDEF(file_out,'anomaliesSON',[lonid,latid],/FLOAT)

NCDF_ATTPUT,file_out,latsvar,'long_name','Latitude'
NCDF_ATTPUT,file_out,latsvar,'units','Degrees'
NCDF_ATTPUT,file_out,latsvar,'valid_min',-90.
NCDF_ATTPUT,file_out,latsvar,'valid_max',90.
NCDF_ATTPUT,file_out,lonsvar,'long_name','Longitude'
NCDF_ATTPUT,file_out,lonsvar,'units','Degrees'
NCDF_ATTPUT,file_out,lonsvar,'valid_min',-180.
NCDF_ATTPUT,file_out,lonsvar,'valid_max',180.

NCDF_ATTPUT,file_out,qvarDJF,'units',unitees
NCDF_ATTPUT,file_out,qvarDJF,'long_name','DJF average anomalies'
NCDF_ATTPUT,file_out,qvarDJF,'missing_value',mdi

NCDF_ATTPUT,file_out,qvarMAM,'units',unitees
NCDF_ATTPUT,file_out,qvarMAM,'long_name','MAM average anomalies'
NCDF_ATTPUT,file_out,qvarMAM,'missing_value',mdi

NCDF_ATTPUT,file_out,qvarJJA,'units',unitees
NCDF_ATTPUT,file_out,qvarJJA,'long_name','JJA average anomalies'
NCDF_ATTPUT,file_out,qvarJJA,'missing_value',mdi

NCDF_ATTPUT,file_out,qvarSON,'units',unitees
NCDF_ATTPUT,file_out,qvarSON,'long_name','SON average anomalies'
NCDF_ATTPUT,file_out,qvarSON,'missing_value',mdi

NCDF_CONTROL,file_out,/ENDEF
NCDF_VARPUT,file_out,qvarDJF,maparrDJF
NCDF_VARPUT,file_out,qvarMAM,maparrMAM
NCDF_VARPUT,file_out,qvarJJA,maparrJJA
NCDF_VARPUT,file_out,qvarSON,maparrSON
NCDF_CLOSE,file_out

return

end

