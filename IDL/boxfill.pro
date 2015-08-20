pro boxfill,data,lon,lat,colors=colors,levels=levels,nofill_color=nofill_color

; give gridbox centres

if not keyword_set(nofill_color) then nofill_color=255
;print,nofill_color
;*mdi=-1.e+30
lond=(lon(1)-lon(0))/2.
latd=abs(lat(1)-lat(0))/2.
for i=0,n_elements(lon)-1 do begin  ;changed from -2 
  for j=0,n_elements(lat)-1 do begin
; j was j=1,n_elements but this looked a bit odd - 0 seems to work better - last loop covers j=0 but with half boxes
; perhaps this is to cope with global data

	s=min(where(levels gt data(i,j)))-1
print,s,data(i,j)
;*	if (data(i,j) NE mdi) THEN BEGIN
	  if colors(s(0)) ne nofill_color then $
	  polyfill,[lon(i)-lond,lon(i)-lond,lon(i)+lond,lon(i)+lond,lon(i)-lond],$
		[lat(j)-latd,lat(j)+latd,lat(j)+latd,lat(j)-latd,lat(j)-latd],$
		color=colors(s)
;	print,lon(i)-lond,lon(i)-lond,lon(i)+lond,lon(i)+lond,lon(i)-lond,$
;		lat(j)-latd,lat(j)+latd,lat(j)+latd,lat(j)-latd,lat(j)-latd
;*	endif
  endfor
  j=0
  s=min(where(levels gt data(i,j)))-1
;*  if(data(i,j) NE mdi) THEN BEGIN
    if colors(s(0)) ne nofill_color then $
    polyfill,[lon(i)-lond,lon(i)-lond,lon(i)+lond,lon(i)+lond,lon(i)-lond],$
    [lat(j)-latd,lat(j),lat(j),lat(j)-latd,lat(j)-latd],$
    color=colors(s)
;*  endif
endfor
end
