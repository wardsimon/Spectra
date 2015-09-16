function rc_test(file)
%
% MATLAB function to calculate the resolution function
% of a triple axis
%
% D.F.M.

%----- f converts from energy units into k^2
%      f=0.4826 for meV and f=1.996854 for THz

f=0.4826;


%----- Read in parameter file

if nargin ~=0

   fid=fopen(file,'r');
   data=[];
   [p]=fscanf(fid,'%f');
   fclose(fid);

else

   hrc_paras=get(findobj('Tag','hrc_paras'),'Userdata');

   for i=1:length(hrc_paras)
       p(i)=str2num(get(hrc_paras(i),'String')); 
   end

end

%----- Calculate Q vector

[Q2c,Qmag]=rc_re2rc([p(19) p(20) p(21)], [p(22) p(23) p(24)], [p(31) p(32) p(33)]);

%----- Calculate resolution matrix in Q frame

[R0,NP,vi,vf]=rc_cnmat(f,Qmag,p,1);

%----- Test of Normalisation constant

[R0P,NPP]=rc_int(1,R0,NP);
[R0P,NPP]=rc_int(1,R0P,NPP);
[R0P,NPP]=rc_int(1,R0P,NPP);


disp(' Resolution Matrix in frame Qx, Qy, Qz')

NP

%----- Calculate Bragg widths

disp(' Bragg widths')

[bragw]=rc_bragg(NP)

%----- Diagonalise resolution matrix in Q frame

[a,b]=eig(NP);

a'*NP*a;

inv(a)'*b*inv(a);

%----- Now work out transformations

Q2c;

A1=[p(25) p(26) p(27)]';
A2=[p(28) p(29) p(30)]';

V1=Q2c*A1;
V2=Q2c*A2;

%----- Form unit vectors V1, V2, V3 in scattering plane

V3=cross(V1,V2);
V2=cross(V3,V1);
V3=V3/sqrt(sum(V3.*V3));
V2=V2/sqrt(sum(V2.*V2));
V1=V1/sqrt(sum(V1.*V1));

U=[V1';V2';V3'];

%----- S transformation matrix from (h,k,l) to V1,V2,V3

S=U*Q2c;
S_save=S;

%----- Work out angle of Q wrt to V1, V2

TT=S*[p(31) p(32) p(33)]';
cos_theta=TT(1)/sqrt(sum(TT.*TT));
sin_theta=TT(2)/sqrt(sum(TT.*TT));

%----- Rotation matrix from Q to V1,V2,V3

R=[cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1];

T=zeros(4,4);
T(4,4)=1;
T(1:3,1:3)=R*S;

disp(' Resolution matrix in frame V1, V2, V3')

T'*NP*T

% ===================== Phonon Widths =================================

%----- Transform scan and dispersion vectors into A-1 and meV vectors
%      in the V1,V2,V3 coordinate frame. 

DQ_cart=(T*[p(35) p(36) p(37) p(38)]')';

scan=DQ_cart/sqrt(sum(DQ_cart.*DQ_cart));
plane=(T(1:3,1:3)*[p(39) p(40) p(41)]')';
plane=[p(42)*plane/sqrt(sum(plane.*plane)) -1];
plane=plane/sqrt(sum(plane.*plane));

% ----- Calculate the phonon width.

[phoni,phonw]=rc_phon(1,NP,plane);
phonw=abs((phonw/(sum(scan.*plane))))

% =====================================================================

%----- Diagonalised matrix in RLU

[a,b]=eig(T'*NP*T);

%----- Projections of the ellipse in Qx, Qy, Qz, W space

A=NP;
const=1.17741;

%----- Remove the vertical component from the matrix.

B=[A(1,1:2),A(1,4);A(2,1:2),A(2,4);A(4,1:2),A(4,4)];

%----- Work out projections for different cuts through the ellipse

rc_projs=figure('Name','Rescal: Projections',...
		'NumberTitle','off',...
		'Visible','on');

%----- Set WYSIWYG characteristics

unis = get(rc_projs,'units');
ppos = get(rc_projs,'paperposition');
set(rc_projs,'units',get(rc_projs,'paperunits'));
pos = get(rc_projs,'position');
pos(3:4) = ppos(3:4);
set(rc_projs,'position',pos);
set(rc_projs,'units',unis);

%----- Write parameters to figure

rc_lab(rc_projs,bragw,phonw,R0,vi,vf)

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
plot(x,y)
fill(x,y,'r')
sxlabel('\14\times Q_x (Å^{-1})')
sylabel('\14\times Q_y (Å^{-1})')

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
plot(x,y)
fill(x,y,'r')
sxlabel('\14\times Q_x (Å^{-1})')
sylabel('\14\times Energy (meV)')


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
plot(x,y)
fill(x,y,'r')
sxlabel('\14\times Q_y (Å^{-1})')
sylabel('\14\times Energy (meV)')


