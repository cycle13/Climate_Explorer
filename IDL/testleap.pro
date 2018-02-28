function testleap,year

; function to test if a year is a leap year
; returns 0.0 if it is a leap year
; returns a non-zero number if it is not a leap year

; first test - is it divisible by 4?
leapoo=(year/4.)-round(year/4.)	

; second test - if it is divisible by 100. then is it also divisible by 400?

FOR i=0,n_elements(leapoo)-1 DO BEGIN
  IF (((year(i)/100.) - round(year(i)/100.)) EQ 0.) THEN BEGIN
    leapoo(i)=leapoo(i)+((year(i)/400.) - round(year(i)/400.))
  END
ENDFOR  

return,leapoo

end
