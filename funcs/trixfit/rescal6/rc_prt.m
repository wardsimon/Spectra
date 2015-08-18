function rc_prt
figure(findobj('name','Rescal: Projections'))
print temp.ps -dps
pause(1)
!lpr -Pps temp.ps
