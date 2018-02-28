pro plot_GIRLSMATHSGRAPH_SEP2013

;---------------------------------------------------
; set up directories and filenames


odir='/home/h04/hadkw/Desktop/images/'
ofileGLOB='GirlsIntoMaths_GLOBAL_graphMAR2015.eps'	
ofileUK='GirlsIntoMaths_UK_graphMAR2015.eps'	
ofileSW='GirlsIntoMaths_SW_graphMAR2015.eps'	

fyr=1960
lyr=2020
nyears=(lyr-fyr)+1

mdi=-1.e+30


set_plot,'ps'
device,filename=odir+ofileGLOB,/ENCAPSUL,xsize=26,ysize=20,/portrait,/color,/helvetica,/bold
!P.thick=4
!P.font=0
zeros=intarr(nyears+1)
yarr=indgen(nyears)+fyr

tvlct,0,0,0,0
ymax=16.
ymin=13.
yrange=ymax-ymin

ytitle='Annual Temperature (degrees C)'
xtitle='Years'
titlee='Global Temperatures (land & ocean)'
 
 
!P.Position=[0.1,0.1,0.93,0.85] 
plot,zeros,xrange=[0,nyears],yrange=[ymin,ymax],xstyle=5,ystyle=1,ytickname=REPLICATE(' ',7),xtitle=xtitle,$
     color=0,charsize=1.5,charthick=4,ythick=4,xthick=4
axis,yaxis=1,charsize=1.1,charthick=3,ythick=4,yticks=30,yticklayout=1,ytickformat='(f4.1)',yrange=[ymin,ymax]   
axis,yaxis=0,charsize=1.1,charthick=3,ythick=4,yticks=30,yticklayout=1,ytickformat='(f4.1)',yrange=[ymin,ymax]   

XYOUTS,0.5,0.93,titlee,/normal,color=0,charsize=2,charthick=4,alignment=0.5
XYOUTS,0.5,0.01,xtitle,/normal,color=0,charsize=2,charthick=4,alignment=0.5
XYOUTS,0.04,0.515,ytitle,/normal,color=0,charsize=2,charthick=4,orientation=90,alignment=0.5

plots,[0,nyears],[ymin,ymin],color=0
plots,[0,nyears],[ymax,ymax],color=0

FOR loo=1,nyears DO BEGIN
  PLOTS,[loo,loo],[ymin,ymax],color=0,thick=1
  PLOTS,[loo,loo],[ymin,ymin+(0.02*yrange)],color=0
  PLOTS,[loo,loo],[ymax,ymax-(0.02*yrange)],color=0
  XYOUTS,loo-0.1,ymin-(0.06*yrange),strcompress(yarr(loo-1),/remove_all),color=0,charsize=0.8,orientation=90  
  XYOUTS,loo-0.1,ymax+(0.01*yrange),strcompress(yarr(loo-1),/remove_all),color=0,charsize=0.8,orientation=90  
ENDFOR
  
  
FOR loo=ymin,ymax,0.1 DO BEGIN
  print,loo
  PLOTS,[0,nyears],[loo,loo],color=0,thick=1
ENDFOR  
device,/close


set_plot,'ps'
device,filename=odir+ofileUK,/ENCAPSUL,xsize=26,ysize=20,/portrait,/color,/helvetica,/bold
!P.thick=4
!P.font=0
zeros=intarr(nyears+1)
yarr=indgen(nyears)+fyr

tvlct,0,0,0,0
ymax=10.
ymin=7.
yrange=ymax-ymin

ytitle='Annual Temperature (degrees C)'
xtitle='Years'
titlee='UK Temperatures'
  
!P.Position=[0.1,0.1,0.93,0.85] 
plot,zeros,xrange=[0,nyears],yrange=[ymin,ymax],xstyle=5,ystyle=1,ytickname=REPLICATE(' ',7),xtitle=xtitle,$
     color=0,charsize=1.5,charthick=4,ythick=4,xthick=4
axis,yaxis=1,charsize=1.1,charthick=3,ythick=4,yticks=30,yticklayout=1,ytickformat='(f4.1)',yrange=[ymin,ymax]   
axis,yaxis=0,charsize=1.1,charthick=3,ythick=4,yticks=30,yticklayout=1,ytickformat='(f4.1)',yrange=[ymin,ymax]   

XYOUTS,0.5,0.93,titlee,/normal,color=0,charsize=2,charthick=4,alignment=0.5
XYOUTS,0.5,0.01,xtitle,/normal,color=0,charsize=2,charthick=4,alignment=0.5
XYOUTS,0.04,0.515,ytitle,/normal,color=0,charsize=2,charthick=4,orientation=90,alignment=0.5

plots,[0,nyears],[ymin,ymin],color=0
plots,[0,nyears],[ymax,ymax],color=0

FOR loo=1,nyears DO BEGIN
  PLOTS,[loo,loo],[ymin,ymax],color=0,thick=1
  PLOTS,[loo,loo],[ymin,ymin+(0.02*yrange)],color=0
  PLOTS,[loo,loo],[ymax,ymax-(0.02*yrange)],color=0
  XYOUTS,loo-0.1,ymin-(0.06*yrange),strcompress(yarr(loo-1),/remove_all),color=0,charsize=0.8,orientation=90  
  XYOUTS,loo-0.1,ymax+(0.01*yrange),strcompress(yarr(loo-1),/remove_all),color=0,charsize=0.8,orientation=90  
ENDFOR
  
  
FOR loo=ymin,ymax,0.1 DO BEGIN
  print,loo
  PLOTS,[0,nyears],[loo,loo],color=0,thick=1
ENDFOR  
device,/close

set_plot,'ps'
device,filename=odir+ofileSW,/ENCAPSUL,xsize=26,ysize=20,/portrait,/color,/helvetica,/bold
!P.thick=4
!P.font=0
zeros=intarr(nyears+1)
yarr=indgen(nyears)+fyr

tvlct,0,0,0,0
ymax=11.
ymin=8.
yrange=ymax-ymin

ytitle='Annual Temperature (degrees C)'
xtitle='Years'
titlee='SW England and S Wales Temperatures'
 
 
!P.Position=[0.1,0.1,0.93,0.85] 
plot,zeros,xrange=[0,nyears],yrange=[ymin,ymax],xstyle=5,ystyle=1,ytickname=REPLICATE(' ',7),xtitle=xtitle,$
     color=0,charsize=1.5,charthick=4,ythick=4,xthick=4
axis,yaxis=1,charsize=1.1,charthick=3,ythick=4,yticks=30,yticklayout=1,ytickformat='(f4.1)',yrange=[ymin,ymax]   
axis,yaxis=0,charsize=1.1,charthick=3,ythick=4,yticks=30,yticklayout=1,ytickformat='(f4.1)',yrange=[ymin,ymax]   

XYOUTS,0.5,0.93,titlee,/normal,color=0,charsize=2,charthick=4,alignment=0.5
XYOUTS,0.5,0.01,xtitle,/normal,color=0,charsize=2,charthick=4,alignment=0.5
XYOUTS,0.04,0.515,ytitle,/normal,color=0,charsize=2,charthick=4,orientation=90,alignment=0.5

plots,[0,nyears],[ymin,ymin],color=0
plots,[0,nyears],[ymax,ymax],color=0

FOR loo=1,nyears DO BEGIN
  PLOTS,[loo,loo],[ymin,ymax],color=0,thick=1
  PLOTS,[loo,loo],[ymin,ymin+(0.02*yrange)],color=0
  PLOTS,[loo,loo],[ymax,ymax-(0.02*yrange)],color=0
  XYOUTS,loo-0.1,ymin-(0.06*yrange),strcompress(yarr(loo-1),/remove_all),color=0,charsize=0.8,orientation=90  
  XYOUTS,loo-0.1,ymax+(0.01*yrange),strcompress(yarr(loo-1),/remove_all),color=0,charsize=0.8,orientation=90  
ENDFOR
  
  
FOR loo=ymin,ymax,0.1 DO BEGIN
  print,loo
  PLOTS,[0,nyears],[loo,loo],color=0,thick=1
ENDFOR  

device,/close

stop
end

