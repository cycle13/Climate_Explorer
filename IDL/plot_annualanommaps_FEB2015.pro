pro plot_annualanommaps_FEB2015,year,filein,fileout,plotout,climst,climed,param


CASE param OF
  'q': BEGIN
    titlees='Surface Specific Humidity'
    units='(g kg!E-1!N)'    
    unitees='g/kg'
    kcolsarrST=[100,2,3,4,5,6,7,9,10,12,13,14,15,16,17]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.0,-2.0,-1.0,-0.75,-0.5,-0.25,-0.1,  0.0,  0.1,0.25,0.5,0.75,1.0,2.0,10]
  END
  'RH': BEGIN
    titlees='Surface Relative Humidity'
    units='(%rh)'   
    unitees='%rh' 
    kcolsarrST=[100,2,3,4,5,6,7,9,10,12,13,14,15,16,17]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-100.0,-10.0,-7.5,-5.0,-2.5,-1.0,-0.5,  0.0, 0.5,1.0,2.5,5.,7.5,10.0,100]
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
  'w': BEGIN
    titlees='Surface Ocean Wind Speed'
    units='(m s!E-1!N)'    
    unitees='m/s'    
    kcolsarrST=[100,2,3,4,5,6,7,9,10,12,13,14,15,16,17]	; 14 colours + black to ensure no botching
    levsarrST=[-2e+30,-10.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.25,  0.0,  0.25,0.5,1.0,2.0,3.0,4.0,10]
  END

ENDCASE

mdi=-1e30
edyr=2014
styr=1973	; 1973, 1971
nyrs=(edyr+1)-styr
nmons=nyrs*12
nlons=72	; 72, 360, 90
nlats=36	; 36, 180, 45
maparr=make_array(nlons,nlats,/float,value=mdi)

latlg=5	;5., 1, 4, 2
lonlg=5	;5., 1, 4, 2
stlt=-90+(latlg/2.)
stln=-180+(lonlg/2.)
nbox=LONG(nlats*nlons)
 
lats=(findgen(nlats)*latlg)+stlt
lons=(findgen(nlons)*lonlg)+stln

getyrmonst=(year-styr)*12
getyrmoned=(((year+1)-styr)*12)-1

inn=NCDF_OPEN(filein)
CASE param OF
;  'q': varid=NCDF_VARID(inn,'anomalies')
;  'q': varid=NCDF_VARID(inn,'q_anoms')
;  'q': varid=NCDF_VARID(inn,'blendmask_q_anoms')
  'q': varid=NCDF_VARID(inn,'blend_q_anoms')
;  'RH': varid=NCDF_VARID(inn,'anomalies')
  'RH': varid=NCDF_VARID(inn,'rh_anoms')
  'e': varid=NCDF_VARID(inn,'e_anoms')
  'Td': varid=NCDF_VARID(inn,'td_anoms')
  'Tw': varid=NCDF_VARID(inn,'tw_anoms')
;  'T': varid=NCDF_VARID(inn,'t_anoms')
  'T': varid=NCDF_VARID(inn,'temperature_anomaly')
  'DPD': varid=NCDF_VARID(inn,'dpd_anoms')
;  'w': varid=NCDF_VARID(inn,'w_anoms')
;  'w': varid=NCDF_VARID(inn,'mask_w_anoms')
;  'w': varid=NCDF_VARID(inn,'Q_w_anoms')
  'w': varid=NCDF_VARID(inn,'sp')
ENDCASE
latid=NCDF_VARID(inn,'latitude')
lonid=NCDF_VARID(inn,'longitude')
NCDF_VARGET,inn,varid,tmp
NCDF_VARGET,inn,latid,lats
NCDF_VARGET,inn,lonid,lons
NCDF_CLOSE,inn

; reduce to desired time frame
;stop
loctims=n_elements(tmp(0,0,*))
print,loctims
tmp=tmp(*,*,loctims-nmons:loctims-1)

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
subarr=tmp(*,*,getyrmonst:getyrmoned)

; average over the year for each gridbox

FOR ltt=0,nlats-1 DO BEGIN
  FOR lnn=0,nlons-1 DO BEGIN
    monalls=subarr(lnn,ltt,*)
    gots=WHERE(monalls NE mdi,count)
    IF (count GE 6) THEN maparr(lnn,ltt)=MEAN(monalls(gots))
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
labsarrST(1:nlevs-2)=string(levsarrST(2:nlevs-1),format='(f6.2)')

;------------------------------------------------------------------------------------
xpos1=0.07
xpos2=0.93
ypos1=0.12
ypos2=0.92
;------------------------------------------------------------------------------------
;make the plot
set_plot,'PS'
device,filename=plotout+'.eps',/ENCAPSUL,$
       xsize=26,ysize=18,/landscape,/color,/helvetica,/bold


!P.Font=0
!P.Thick=4
!X.Thick=4
!Y.Thick=4

tvlct,200,200,200,100
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
ENDIF ELSE IF (param EQ 'w') THEN BEGIN
  tvlct,0,0,0,1
  tvlct,0,0,0,2
  tvlct,0,0,155,3
  tvlct,7,90,255,4
  tvlct,50,118,255,5
  tvlct,89,144,255,6
  tvlct,140,178,255,7
  tvlct,191,212,255,8
  tvlct,229,238,255,9
  tvlct,255,255,153,10
  tvlct,255,255,0,11
  tvlct,255,204,0,12
  tvlct,255,153,0,13
  tvlct,255,102,0,14
  tvlct,255,0,0,15
  tvlct,155,0,0,16
  tvlct,0,0,0,17
  tvlct,0,0,0,18
ENDIF


!P.Position=[xpos1,ypos1,xpos2,ypos2]

Map_set,/continents,color=0,/noborder,/noerase,/robinson,/grid, glinestyle=2

XYOUTS,0.5,ypos2+0.02,titlees+' in '+string(year,format='(a4)'),/NORMAL,alignment=0.5,charsize=2,color=0

!P.Position=[xpos1,ypos1,xpos2,ypos2]
boxfill,maparr,lons,lats,colors=kcolsarrST,levels=levsarrST	    ;,nofill_color=kcolsarrT(0)

Map_set,/continents,color=0,/noborder,/noerase,/robinson,/grid, glinestyle=2

PLOTS,[-179.9,0],[89.9,89.9],color=0
PLOTS,[0,179.9],[89.9,89.9],color=0
PLOTS,[-179.9,0],[-89.9,-89.9],color=0
PLOTS,[0,179.9],[-89.9,-89.9],color=0

MAKE_KEY,0.1,0.07,0.8,0.02,0.0,-0.03,/NORMAL,COLORS=colsarrST,labels=labsarrST,$
         charsize=1.2,charthick=4,bcolor=0
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
qvar=NCDF_VARDEF(file_out,'anomalies',[lonid,latid],/FLOAT)

NCDF_ATTPUT,file_out,latsvar,'long_name','Latitude'
NCDF_ATTPUT,file_out,latsvar,'units','Degrees'
NCDF_ATTPUT,file_out,latsvar,'valid_min',-90.
NCDF_ATTPUT,file_out,latsvar,'valid_max',90.
NCDF_ATTPUT,file_out,lonsvar,'long_name','Longitude'
NCDF_ATTPUT,file_out,lonsvar,'units','Degrees'
NCDF_ATTPUT,file_out,lonsvar,'valid_min',-180.
NCDF_ATTPUT,file_out,lonsvar,'valid_max',180.

NCDF_ATTPUT,file_out,qvar,'units',unitees
NCDF_ATTPUT,file_out,qvar,'long_name','Annual average anomalies'
NCDF_ATTPUT,file_out,qvar,'missing_value',mdi

NCDF_CONTROL,file_out,/ENDEF
NCDF_VARPUT,file_out,qvar,maparr
NCDF_CLOSE,file_out

return

end

