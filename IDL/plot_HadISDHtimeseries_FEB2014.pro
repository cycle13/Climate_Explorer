pro plot_HadISDHtimeseries_FEB2014

; add a comparison of regional time series and trends from the raw and homogenised data.



;-----------------------------------------------------
!Except=2

version='2.0.1.2014p'
nowmon='JAN'
nowyear='2015'
thenmon='JAN'
thenyear='2015'
homogtype='IDPHA'
param='RH'

; files and directories
dir='/data/local/hadkw/HADCRUH2/UPDATE2014/'
filein='STATISTICS/TIMESERIES/HadISDH.land'+param+'.'+version+'_FLATgrid'+homogtype+'5by5_'+thenmon+thenyear+'_areaTS_19732014.nc'
fileout='IMAGES/TIMESERIES/HadISDH.land'+param+'.'+version+'_FLATgrid'+homogtype+'5by5_'+nowmon+nowyear+'_areaTS_19732014.eps'

;--------------------------------------------------------
; variables and arrays
mdi=-1e+30

styr=1973
edyr=2014
nyrs=(edyr+1)-styr
nmons=nyrs*12
int_mons=indgen(nmons)

CASE param OF
  'q': BEGIN
    titlees='Surface Specific Humidity'
    units='g kg!E-1!N' 
    ymin=-1
    ymax=1   
  END
  'RH': BEGIN
    titlees='Surface Relative Humidity'
    units='%rh'   
    ymin=-4.5
    ymax=4.5   
  END
  'e': BEGIN
    titlees='Surface Vapour Pressure'
    units='hPa'  
    ymin=-2
    ymax=2   
  END
  'Td': BEGIN
    titlees='Surface Dew Point Temperature'
    units='!Eo!NC'
    ymin=-2
    ymax=2   
  END
  'Tw': BEGIN
    titlees='Surface Wet Bulb Temperature'
    units='!Eo!NC'
    ymin=-2.5
    ymax=2   
  END
  'T': BEGIN
    titlees='Surface Temperature'
    units='!Eo!NC'    
    ymin=-2
    ymax=2   
  END
  'DPD': BEGIN
    titlees='Surface Dew Point Depression'
    units='!Eo!NC'    
    ymin=-1.5
    ymax=1.5   
  END
ENDCASE


;---------------------------------------------------------
; open area ave TS files
inn=NCDF_OPEN(dir+filein)
CASE param OF
  'q': BEGIN
    varg=NCDF_VARID(inn,'glob_q_anoms')
    varnh=NCDF_VARID(inn,'nhem_q_anoms')
    vart=NCDF_VARID(inn,'trop_q_anoms') 
    varsh=NCDF_VARID(inn,'shem_q_anoms')
  END
  'RH': BEGIN
    varg=NCDF_VARID(inn,'glob_RH_anoms')
    varnh=NCDF_VARID(inn,'nhem_RH_anoms')
    vart=NCDF_VARID(inn,'trop_RH_anoms') 
    varsh=NCDF_VARID(inn,'shem_RH_anoms')
  END
  'e': BEGIN
    varg=NCDF_VARID(inn,'glob_e_anoms')
    varnh=NCDF_VARID(inn,'nhem_e_anoms')
    vart=NCDF_VARID(inn,'trop_e_anoms') 
    varsh=NCDF_VARID(inn,'shem_e_anoms')
  END
  'Td': BEGIN
    varg=NCDF_VARID(inn,'glob_Td_anoms')
    varnh=NCDF_VARID(inn,'nhem_Td_anoms')
    vart=NCDF_VARID(inn,'trop_Td_anoms') 
    varsh=NCDF_VARID(inn,'shem_Td_anoms')
  END
  'Tw': BEGIN
    varg=NCDF_VARID(inn,'glob_Tw_anoms')
    varnh=NCDF_VARID(inn,'nhem_Tw_anoms')
    vart=NCDF_VARID(inn,'trop_Tw_anoms') 
    varsh=NCDF_VARID(inn,'shem_Tw_anoms')
  END
  'T': BEGIN
    varg=NCDF_VARID(inn,'glob_T_anoms')
    varnh=NCDF_VARID(inn,'nhem_T_anoms')
    vart=NCDF_VARID(inn,'trop_T_anoms') 
    varsh=NCDF_VARID(inn,'shem_T_anoms')
  END
  'DPD': BEGIN
    varg=NCDF_VARID(inn,'glob_DPD_anoms')
    varnh=NCDF_VARID(inn,'nhem_DPD_anoms')
    vart=NCDF_VARID(inn,'trop_DPD_anoms') 
    varsh=NCDF_VARID(inn,'shem_DPD_anoms')
  END
ENDCASE
NCDF_VARGET,inn,varg,globval
NCDF_VARGET,inn,varnh,nhemval
NCDF_VARGET,inn,vart,tropval
NCDF_VARGET,inn,varsh,shemval
NCDF_CLOSE,inn


set_plot,'PS'
device,filename=dir+fileout,/color,/ENCAPSUL,xsize=20,ysize=26,/portrait,/helvetica,/bold
!P.Font=0
!P.Thick=4

!P.Position=[0.04,0.52,0.44,0.94]

;-------------------------------------
;the time series plots
tvlct,0,0,0,0		;black			; HadCRUH2


  tvlct,102,47,0,1	; BROWN
  tvlct,0,122,153,2	; BLUE

  tvlct,40,87,255,3	; BLUE
  tvlct,255,61,61,4	; RED

; YEAR***
fulltims=indgen(nmons+1)
zeros=intarr(nmons+1)

yyrange=ymax-ymin

allnames=['1975','1980','1985','1990','1995','2000','2005','2010']

x1pos=[0.1,0.1,0.1,0.1]
x2pos=[0.95,0.95,0.95,0.95]
y1pos=[0.74,0.51,0.28,0.05]
y2pos=[0.97,0.74,0.51,0.28]

IF (param EQ 'T') THEN BEGIN
  collow=3
  colhigh=4
ENDIF ELSE IF (param EQ 'DPD') THEN BEGIN
  collow=2
  colhigh=1
ENDIF ELSE BEGIN
  collow=1
  colhigh=2
ENDELSE


!P.Position=[x1pos(0),y1pos(0),x2pos(0),y2pos(0)]
plot,fulltims,zeros,min_value=-100,ystyle=1,xstyle=5,yrange=[ymin,ymax],ythick=4,charsize=0.9,$
     ytitle='anomalies ('+units+')',charthick=1,color=0,/noerase,thick=4,/nodata	
FOR i=0,nmons-1 DO BEGIN
  IF (globval(i) NE mdi) AND (globval(i) LT 0.) THEN PLOTS,[i,i],[0.,globval(i)],color=collow 
  IF (globval(i) NE mdi) AND (globval(i) GE 0.) THEN PLOTS,[i,i],[0.,globval(i)],color=colhigh 
ENDFOR

oplot,fulltims(0:nmons),zeros,min_value=-100,color=0,thick=1
;oplot,fulltims(0:nmons-1),globval,min_value=-100,color=0,thick=4
PLOTS,[0,nmons],[ymax,ymax],color=0	; YEAR***
PLOTS,[0,nmons],[ymin,ymin],color=0	; YEAR***
FOR tt=1,nyrs-1 DO BEGIN
  IF (((tt+3.)/5.)-FLOOR(((tt+3.)/5.)) GT 0.) THEN BEGIN
    PLOTS,[tt*12,tt*12],[ymin,ymin*0.97],color=0,thick=4 
    PLOTS,[tt*12,tt*12],[ymax,ymax*0.97],color=0,thick=4 
  ENDIF ELSE BEGIN
    PLOTS,[tt*12,tt*12],[ymin,ymin*0.94],color=0,thick=4
    PLOTS,[tt*12,tt*12],[ymax,ymax*0.94],color=0,thick=4
  ENDELSE
ENDFOR

XYOUTS,0.5,y2pos(0)-0.03,'GLOBE (70S-70N)',/normal,alignment=0.5,color=0,charsize=1.
trend=median_pairwise(globval,mdi,se,lc,uc)
XYOUTS,0.94,y1pos(0)+0.01,string(trend*120.,format='(f5.2)')+' '+units+' dec!E-1!N ('+string(lc*120.,format='(f5.2)')+$
                 ' to '+string(uc*120.,format='(f5.2)')+')',/normal,color=0,charsize=1.2,alignment=1


!P.Position=[x1pos(1),y1pos(1),x2pos(1),y2pos(1)]
plot,fulltims,zeros,min_value=-100,ystyle=1,xstyle=5,yrange=[ymin,ymax],ythick=4,charsize=0.9,$
     ytitle='anomalies ('+units+')',charthick=1,color=0,/noerase,thick=4,/nodata	
FOR i=0,nmons-1 DO BEGIN
  IF (nhemval(i) NE mdi) AND (nhemval(i) LT 0.) THEN PLOTS,[i,i],[0.,nhemval(i)],color=collow 
  IF (nhemval(i) NE mdi) AND (nhemval(i) GE 0.) THEN PLOTS,[i,i],[0.,nhemval(i)],color=colhigh 
ENDFOR
oplot,fulltims(0:nmons),zeros,min_value=-100,color=0,thick=1
;oplot,fulltims(0:nmons-1),nhemval,min_value=-100,color=0,thick=4
PLOTS,[0,nmons],[ymax,ymax],color=0	; YEAR***
PLOTS,[0,nmons],[ymin,ymin],color=0	; YEAR***
FOR tt=1,nyrs-1 DO BEGIN
  IF (((tt+3.)/5.)-FLOOR(((tt+3.)/5.)) GT 0.) THEN BEGIN
    PLOTS,[tt*12,tt*12],[ymin,ymin*0.97],color=0,thick=4 
    PLOTS,[tt*12,tt*12],[ymax,ymax*0.97],color=0,thick=4 
  ENDIF ELSE BEGIN
    PLOTS,[tt*12,tt*12],[ymin,ymin*0.94],color=0,thick=4
    PLOTS,[tt*12,tt*12],[ymax,ymax*0.94],color=0,thick=4
  ENDELSE
ENDFOR

XYOUTS,0.5,y2pos(1)-0.03,'N HEMISPHERE (20N-70N)',/normal,alignment=0.5,color=0,charsize=1.
trend=median_pairwise(nhemval,mdi,se,lc,uc)
XYOUTS,0.94,y1pos(1)+0.01,string(trend*120.,format='(f5.2)')+' '+units+' dec!E-1!N ('+string(lc*120.,format='(f5.2)')+$
                 ' to '+string(uc*120.,format='(f5.2)')+')',/normal,color=0,charsize=1.2,alignment=1


!P.Position=[x1pos(2),y1pos(2),x2pos(2),y2pos(2)]
plot,fulltims,zeros,min_value=-100,ystyle=1,xstyle=5,yrange=[ymin,ymax],ythick=4,charsize=0.9,$
     ytitle='anomalies ('+units+')',charthick=1,color=0,/noerase,thick=4,/nodata	
FOR i=0,nmons-1 DO BEGIN
  IF (tropval(i) NE mdi) AND (tropval(i) LT 0.) THEN PLOTS,[i,i],[0.,tropval(i)],color=collow 
  IF (tropval(i) NE mdi) AND (tropval(i) GE 0.) THEN PLOTS,[i,i],[0.,tropval(i)],color=colhigh 
ENDFOR
oplot,fulltims(0:nmons),zeros,min_value=-100,color=0,thick=1
;oplot,fulltims(0:nmons-1),tropval,min_value=-100,color=0,thick=4
PLOTS,[0,nmons],[ymax,ymax],color=0	; YEAR***
PLOTS,[0,nmons],[ymin,ymin],color=0	; YEAR***
FOR tt=1,nyrs-1 DO BEGIN
  IF (((tt+3.)/5.)-FLOOR(((tt+3.)/5.)) GT 0.) THEN BEGIN
    PLOTS,[tt*12,tt*12],[ymin,ymin*0.97],color=0,thick=4 
    PLOTS,[tt*12,tt*12],[ymax,ymax*0.97],color=0,thick=4 
  ENDIF ELSE BEGIN
    PLOTS,[tt*12,tt*12],[ymin,ymin*0.94],color=0,thick=4
    PLOTS,[tt*12,tt*12],[ymax,ymax*0.94],color=0,thick=4
  ENDELSE
ENDFOR

XYOUTS,0.5,y2pos(2)-0.03,'TROPICS (20S-20N)',/normal,alignment=0.5,color=0,charsize=1.
trend=median_pairwise(tropval,mdi,se,lc,uc)
XYOUTS,0.94,y1pos(2)+0.01,string(trend*120.,format='(f5.2)')+' '+units+' dec!E-1!N ('+string(lc*120.,format='(f5.2)')+$
                 ' to '+string(uc*120.,format='(f5.2)')+')',/normal,color=0,charsize=1.2,alignment=1


!P.Position=[x1pos(3),y1pos(3),x2pos(3),y2pos(3)]
plot,fulltims,zeros,min_value=-100,ystyle=1,xstyle=5,yrange=[ymin,ymax],ythick=4,charsize=0.9,$
     ytitle='anomalies ('+units+')',charthick=1,color=0,/noerase,thick=4,/nodata	
FOR i=0,nmons-1 DO BEGIN
  IF (shemval(i) NE mdi) AND (shemval(i) LT 0.) THEN PLOTS,[i,i],[0.,shemval(i)],color=collow 
  IF (shemval(i) NE mdi) AND (shemval(i) GE 0.) THEN PLOTS,[i,i],[0.,shemval(i)],color=colhigh 
ENDFOR
oplot,fulltims(0:nmons),zeros,min_value=-100,color=0,thick=1
;oplot,fulltims(0:nmons-1),shemval,min_value=-100,color=0,thick=4
PLOTS,[0,nmons],[ymax,ymax],color=0	; YEAR***
PLOTS,[0,nmons],[ymin,ymin],color=0	; YEAR***
yrcount=0
FOR tt=1,nyrs-1 DO BEGIN
  IF (((tt+3.)/5.)-FLOOR(((tt+3.)/5.)) GT 0.) THEN BEGIN
    PLOTS,[tt*12,tt*12],[ymin,ymin*0.97],color=0,thick=4 
    PLOTS,[tt*12,tt*12],[ymax,ymax*0.97],color=0,thick=4 
  ENDIF ELSE BEGIN
    PLOTS,[tt*12,tt*12],[ymin,ymin*0.94],color=0,thick=4
    PLOTS,[tt*12,tt*12],[ymax,ymax*0.94],color=0,thick=4
    XYOUTS,tt*12,ymin*1.15,allnames(yrcount),alignment=0.5,color=0,charsize=1.
    yrcount=yrcount+1
  ENDELSE
ENDFOR
XYOUTS,0.5,y2pos(3)-0.03,'S HEMISPHERE (70S-20S)',/normal,alignment=0.5,color=0,charsize=1.
trend=median_pairwise(shemval,mdi,se,lc,uc)
XYOUTS,0.94,y1pos(3)+0.01,string(trend*120.,format='(f5.2)')+' '+units+' dec!E-1!N ('+string(lc*120.,format='(f5.2)')+$
                 ' to '+string(uc*120.,format='(f5.2)')+')',/normal,color=0,charsize=1.2,alignment=1

XYOUTS,0.5,0.01,'YEAR',/normal,alignment=0.5,charsize=1.1,color=0


device,/close

end
