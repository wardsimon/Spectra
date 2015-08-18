function rc_projs(R0,NP)
%
% MATLAB function to plot the projections of the resolution ellipse
% of a triple axis
%
% DFM 10.11.95

figure(findobj('Tag','Rescal: Projections'));

A=NP;
const=1.17741;

%----- Remove the vertical component from the matrix.

B=[A(1,1:2),A(1,4);
   A(2,1:2),A(2,4);
   A(4,1:2),A(4,4)];

%----- Work out projections for different cuts through the ellipse

%----- S is the rotation matrix that diagonalises the projected ellipse

%----- 1. Qx, Qy plane

[R0P,MP]=rc_int(3,R0,B);

theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
S=[];
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta);
axes('Position',[0.09 0.1 0.21 0.21])
line(x,y)
fill(x,y,'r')
xlabel('Qx (1/Angs)')
ylabel('Qy (1/Angs)')

%---------------- Add slice through Qx,Qy plane ----------------------

MP=A(1:2,1:2);

theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
S=[];
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta);
hold on
line(x,y)
fill(x,y,'g')
hold off

%----- 2. Qx, W plane

[R0P,MP]=rc_int(2,R0,B);

theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
S=[];
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta);
axes('Position',[0.41 0.1 0.21 0.21])
line(x,y)
fill(x,y,'r')
xlabel('Qx (1/Angs)')
ylabel('Energy (meV)')
%sxlabel('\14\times Q_x (�^{-1})')
%sylabel('\14\times Energy (meV)')

%---------------- Add slice through Qx,W plane ----------------------

MP=[A(1,1) A(1,4);A(4,1) A(4,4)];

theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
S=[];
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta);
hold on
line(x,y)
fill(x,y,'g')
hold off

%----- 3. Qy, W plane

[R0P,MP]=rc_int(1,R0,B);

theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
S=[];
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta);
axes('Position',[0.73 0.1 0.21 0.21])
line(x,y)
fill(x,y,'r')
xlabel('Qy (1/Angs)')
ylabel('Energy (meV)')
%sxlabel('\14\times Q_y (�^{-1})')
%sylabel('\14\times Energy (meV)')

%---------------- Add slice through Qy,W plane ----------------------

MP=[A(2,2) A(2,4);A(4,2) A(4,4)];

theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
S=[];
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta);
hold on
line(x,y)
fill(x,y,'g')
hold off





