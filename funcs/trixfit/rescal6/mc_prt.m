function rc_prt
figure(findobj('Tag','Rescal: Simulation Plot'))
print temp.ps -dps
pause(1)
!print/queue=clps4_av6 temp.ps
