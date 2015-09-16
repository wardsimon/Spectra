% leasqrexamp : leasqr and simplex fit example/test

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: leasqr example

fprintf(1,'\nThis is a fitting example with Leasqr V3.b.');
fprintf(1,'\nCommands to be typed follow...\n');

fprintf(1,'\nCreating time axis, initial parameters and non noised data.');
fprintf(1,'\nFitting function ''leasqrfunc'' is : y=p(1)*exp(-p(2)*x) with p= [ 1 ; 0.1 ].');
fprintf(1,'\nYou should use your own function file ''myfunc.m''');
t = [1:100];
p=[1; 0.1];
data=leasqrfunc(t,p);
fprintf(1,'\n  t = [1:100];\n  p=[1; 0.1];\n  data=leasqrfunc(t,p);');

fprintf(1,'\nAdding Random Noise.\n data=data+0.05*rand(1,100);\n');
data=data+0.05*rand(1,100);

fprintf(1,'\nLaunching Least Square Marquardt-Levenberg Fit.');
fprintf(1,'\n  leasqr(t,data,[.8;.05],''leasqrfunc'')\n');

[a,b,c,d,e]=leasqr(t,data,[.8;.05],'leasqrfunc',[],[],[],[],[],[],'verbose');

fprintf(1,'\nLaunching Simplex Fit.');
simplex(t,data,[.8;.05],'leasqrfunc',[],[],[],[],'normal');
