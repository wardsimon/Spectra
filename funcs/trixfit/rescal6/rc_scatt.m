function rc_scatt
%
% MATLAB function to calculate the resolution function
% of a triple axis
%
% DAT 21.11.97 

% Updating to include :
% 		Typical scan 15 pts.
% 		Show scan direction
%		Give reciprocal space details. 

%----- f converts from energy units into k^2
%      f=0.4826 for meV and f=1.996854 for THz

f=0.4826;

%----- Get method from Parameter window

method=get(findobj('tag','hrc_rescal_method'),'Userdata');

%----- Get parameters from Parameter and Instrumentation window

p=rc_savp('respars');

%----- Calculate Q vector

[Q2c,Qmag]=rc_re2rc([p(19) p(20) p(21)], [p(22) p(23) p(24)], [p(31) p(32) p(33)]);

%----- Calculate resolution matrix in Q frame, and check that scattering triangle closes

[R0,NP,vi,vf,Error]=feval(method,f,Qmag,p,0);

if Error ~= 0; disp('Scattering triangle does not close'), return; end

%----- Test of Normalisation constant

[R0P,NPP]=rc_int(1,R0,NP);
[R0P,NPP]=rc_int(1,R0P,NPP);
[R0P,NPP]=rc_int(1,R0P,NPP);

%----- Calculate Bragg widths

disp(' Bragg widths')

[bragw]=rc_bragg(NP)

%----- Diagonalise resolution matrix in Q frame

[a,b]=eig(NP);
inv(a)'*b*inv(a);

%----- Now work out transformations

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

if~(TT(3)<eps && TT(3)>-eps) 
  disp('Q not in scattering plane!!!!')
end

%----- Rotation matrix from V1,V2,V3 to Q

R=[cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1];

%----- Now make up the T {transform matrix} that 
%----- takes us from Y{hkl}->Y{qx,qy,qz} the Cooper Nathans
%----- coordinates. 


T=zeros(4,4);
T(4,4)=1;
T(1:3,1:3)=R*S;

VR=zeros(4,4);
VR(4,4)=1;
VR(1:3,1:3)=R;

VR'*NP*VR;
VNP=VR'*NP*VR;

%----- Phonon Widths

%----- Transform scan and dispersion vectors into A-1 and meV vectors
%----- in the qx,qy,qz,en coordinate frame (Cooper Nathans frame). 

DQ_cart=(T*[p(35) p(36) p(37) p(38)]')';

dscan=sqrt(sum(DQ_cart.*DQ_cart));
scan=DQ_cart/sqrt(sum(DQ_cart.*DQ_cart));
plane=(T(1:3,1:3)*[p(39) p(40) p(41)]')';
plane=[p(42)*plane/sqrt(sum(plane.*plane)) -1];
plane=plane/sqrt(sum(plane.*plane));

% ----- Calculate the phonon width.

[phoni,phonw]=rc_phon(R0,NP,plane);
phonw=abs((phonw/(sum(scan.*plane))));
disp('Phonon width in (meV,1/Angs mixture)');
disp(phonw)
disp('and intensity ');
disp(phoni)


%----- This is to become a new routine ----------------------------
%----- Used Currently for Development -----------------------------
%------------------------------------------------------------------
%----- Reciprocal space plots -------------------------------------
%------------------------------------------------------------------


if ~isempty(findobj('Tag','Rescal: Reciprocal'))
   close(findobj('Tag','Rescal: Reciprocal'))
end

h=figure('Name','Rescal: Reciprocal',...
                'Tag','Rescal: Reciprocal',...
                'NumberTitle','off',...
                'Menubar','None',...
                'Visible','on','Position',[500 150 500 500]);
figure(h);
whitebg(h); % use white background as it looks better

h1=uimenu(h,'Label','Print');
uimenu(h1,'Label','Print Figure','Callback','print');

%----- Set WYSIWYG characteristics

unis = get(h,'units');
%set(h,'units',get(h,'paperunits'));
%pos = [450 100 700 700];
%set(h,'position',pos);
set(h,'units',unis);

%----------------------------------
dm=p(1);            % monochromator d-spacing in Angs.
da=p(2);            % analyser d-spacing in Angs.
sm=p(6);            % scattering sense of monochromator (left=+1,right=-1)
ss=p(7);            % scattering sense of sample (left=+1,right=-1)
sa=p(8);            % scattering sense of analyser (left=+1,right=-1)
kfix=p(9);          % fixed momentum component in ang-1.
fx=p(10);           % fx=1 for fixed incident and 2 for scattered wavevector.
w=p(34);            % energy transfer.
q0=Qmag; 

% In addition the parameters f, energy pre-multiplier-f*w
% where f=0.48 for meV to ang-2 - and q0 which is the wavevector
% transfer in ang-1, are passed over.

f=0.48;

% Calculate ki and kf

ki=sqrt(kfix^2+(fx-1)*f*w);  % kinematical equations.
kf=sqrt(kfix^2-(2-fx)*f*w);

% Test if scattering triangle is closed
cos_2theta=(ki^2+kf^2-q0^2)/(2*ki*kf);

thetaa=asin(pi/(da*kf));      % theta angles for analyser
thetam=asin(pi/(dm*ki));      % and monochromator.
  


% Now construct a scattering diagram

% work out the scattering vector Q
S=S_save;
Qhkl=[p(31) p(32) p(33)]';
Qv=S*Qhkl;
Qmag=sqrt(Qv'*Qv);

% now calculate angles of the scattering triangle
phi=atan2(Qv(2),Qv(1));
twotheta=acos(cos_2theta);
cosbeta=(kf^2+q0^2-ki^2)/(2*kf*q0)
beta=-ss*acos(cosbeta);

% construct the meeting point Spt of Ki and Kf vectors
Spt_v=[kf*cos(phi+beta), kf*sin(phi+beta), 0];

% plot up the three scattering points

xpts=[0 Qv(1) Spt_v(1) 0]
ypts=[0 Qv(2) Spt_v(2) 0]

% set frame size
prop=2;
Dx=prop*(max(xpts)-min(xpts));
Dy=prop*(max(ypts)-min(ypts));
Cenx=(max(xpts)+min(xpts))/2;
Ceny=(max(ypts)+min(ypts))/2;
dim=max([Dx Dy]);
Dx=dim;Dy=dim;
frame=[Cenx-dim/2 Cenx+dim/2 Ceny-dim/2 Ceny+dim/2];


figure(h);
set(gca,'visible','off');
axes('Position',[0.65 0.60 0.3 0.3])
plot(xpts,ypts,'w*');
axis(frame);

title('Scattering Triangle');
xlabel('V1 (1/�)')
ylabel('V2 (1/�)')
grid on;

% --- Calculate and superimpose RL vectors. ---------
% number of points along x.
% The user defined vectors A1 and A2 are
A1V=S*A1;
A2V=S*A2;
% number of RLPs along x and y:
nx=floor(Dx/abs(A1V(1)))+1;
ny=floor(Dy/abs(A2V(2)))+1;
% number of reciproval lattice points within the frame:
Nrlp=nx*ny;

% Now need to calculate the Irlp(nx,ny) and Jrlp(nx,ny) 
% indices of the Nrlp RLPs.

% Min and max Jrlp can be calculated directly as only
% A2 has a V2 component. 
xmin=frame(1);
xmax=frame(2);
ymin=frame(3);
ymax=frame(4);

minJrlp=min([fix(ymin/A2V(2)) fix(ymax/A2V(2))]);
indJ=minJrlp+[0:ny-1];   % the set of J indices
Jrlp=ones(nx,1)*indJ;

% now construct Irlp
minI0rlp=min([fix(xmin/A1V(1)) fix(xmax/A1V(1)) ]);
indI=minI0rlp+[0:nx-1];
I0rlp=indI'*ones(1,ny);
% have to calculate the shift for each line due to J
Ishift=-ceil(Jrlp*A2V(1)/abs(A1V(1)));
Irlp=I0rlp+Ishift;

% calculate the X and Y components (V1,V2 components)
% of the RLP's

XRrlp=Irlp*A1V(1)+Jrlp*A2V(1);
YRrlp=Irlp*A1V(2)+Jrlp*A2V(2);

% Now we have found all the indices that are contained
% we can plot them onto the figure

hold on;
hRLP=plot(XRrlp,YRrlp,'ko');
set(hRLP,'linewidth',2)
hold off;

% Add text to the RLPs
% Need to make up a text array for each point
%Xtext=reshape(XRrlp',nx*ny,1)+.1;
%Ytext=reshape(YRrlp',nx*ny,1)+.1;
%
%Stext=[];
%for i=1:nx
%  for j=1:ny
%    indrlp=Irlp(i,j)*A1+Jrlp(i,j)*A2;
%%    string=['(' num2str(Irlp(i,j)) ',' num2str(Jrlp(i,j)) ')' ];
%    string=['(' num2str(indrlp(1)) ',' num2str(indrlp(2)) ','...
%    num2str(indrlp(3)) ')' ];
%    Stext=str2mat(Stext,string);
%  end
%end

htext=text(-0.1*dim,-.1*dim,'(0,0,0)');
htext1=text(0.1*dim+Qv(1),-.1*dim+Qv(2),'Q pos');

% and include the reciprocal lattice vectors
% in the bottom left hand corner:
xoffset=frame(1)+.2*dim;
yoffset=frame(3)+.2*dim;
%X=[0 A1V(1)]+xoffset;Y=[0 A1V(2)]+yoffset;
%hA1=gf_arrow(X,Y,[0 1 0],1,'A1');
%X=[0 A2V(1)]+xoffset;Y=[0 A2V(2)]+yoffset;
%hA2=gf_arrow(X,Y,[0 1 0],1,'A2');


% Now draw the scattering vectors out:

hQ=gf_arrow(xpts(1:2),ypts(1:2),[0 0 0],3,'Q');
hKi=gf_arrow(xpts([3 2]),ypts([3 2]),[1 0 0],2,'Ki');
hKf=gf_arrow(xpts(3:4),ypts(3:4),[1 0 0],2,'Kf');
X=[0 A1V(1)]+xoffset;Y=[0 A1V(2)]+yoffset;
%hA1=gf_arrow(X,Y,[0 1 1],2,'A1',.1);
X=[0 A2V(1)]+xoffset;Y=[0 A2V(2)]+yoffset;
%hA2=gf_arrow(X,Y,[0 1 1],2,'A2',-.2);

% ============ work out all the phonon dispersion lines and widths ====

%----- Transform scan and dispersion vectors into A-1 and meV vectors
%----- in the V1,V2,V3,en coordinate frame. 

DQ_V=[(S*[p(35) p(36) p(37)]')' p(38)];

scan_V=DQ_cart/sqrt(sum(DQ_cart.*DQ_cart));
G_V=[(S*[p(39) p(40) p(41)]')' p(42)];
plane_V=(S*[p(39) p(40) p(41)]')';
plane_V=[p(42)*plane_V/sqrt(sum(plane_V.*plane_V)) -1];
plane_V=plane_V/sqrt(sum(plane_V.*plane_V));

% ----- Calculate the phonon width.
disp(plane_V)
xx=[-plane_V(2) plane_V(2)]+Qv(1);
yy=[plane_V(1) -plane_V(1)]+Qv(2);
%hdisp=line(xx,yy);
%set(hdisp,'linewidth',2,'color',[0 0 0],'linestyle','--')


% ======================================================================
% Lets show the resolution ellipsoid now

%----- S is the rotation matrix that diagonalises the projected ellipse


A=VNP;
const=1.17741;

%----- Remove the vertical component from the matrix.

B=[A(1,1:2),A(1,4);A(2,1:2),A(2,4);A(4,1:2),A(4,4)];


%----- 1. Plot on the scattering plane graph

[R0P,MP]=rc_int(3,R0,B);

theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
S=[];
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta);
x=5*x+Qv(1);
y=5*y+Qv(2);
%axes('Position',[0.09 0.1 0.21 0.21])
hold on
line(x,y);
fill(x,y,'r')
hold off
%xlabel('V1 (1/Angs)')
%ylabel('V2 (1/Angs)')

%------

%----- 1. Qx, Qy plane

[R0P,MP]=rc_int(3,R0,B);

theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
S=[];
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta);
axes('Position',[0.15 0.6 0.3 0.3])
line(x,y)
fill(x,y,'r')
axis([-.1 .1 -.1 .1])
xlabel('V1 (1/Angs)')
ylabel('V2 (1/Angs)')
title('Q Orientation')

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
grid on
hold off
 
%--------------- add dispersion plane -------------------------------- 

axis(axis)
hdisp=line(xx-Qv(1),yy-Qv(2));
set(hdisp,'linewidth',2,'color',[0 0 0],'linestyle','--')

% add an arrow indicating scan direction
hscan=gf_arrow([0 5*DQ_V(1)],[0 5*DQ_V(2)],[0 0 0],2,'DQE');

%----- 2. Qx, W plane

% rotate ellipsoid into dispersion frame.
theta=atan2(G_V(1),G_V(2));
RotGtoV=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]
BG=RotGtoV'*B*RotGtoV;
[R0P,MP]=rc_int(2,R0,BG);
DQG_scan=RotGtoV'*[DQ_V(1);DQ_V(2);DQ_V(4)];

theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
S=[];
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta);
axes('Position',[0.15 0.15 0.3 0.3])
line(x,y)
fill(x,y,'r')
xlabel('[GH,GK,GL] (1/Angs)')
ylabel('Energy (meV)')
%sxlabel('\14\times Q_x (�^{-1})')
%sylabel('\14\times Energy (meV)')


axis(axis)
hdisp=line([-1 1],[-G_V(4) G_V(4)]);
set(hdisp,'linewidth',2,'color',[0 0 0],'linestyle','-')
grid on
htext=text(-.08,-1.5,['slope = ' num2str(p(42)) ' meV/Angs'])
title('Q,E Focussing');

gf_arrow([0 2*DQG_scan(1)],[0 2*DQG_scan(3)],[0 0 0],2,'DQE');

phi
twotheta
beta

% --------------------------------------------------------------------------

% Now calculate scan!

disp(DQ_cart)
disp(dscan)
disp(phonw)
axes('Position',[0.65 0.15 0.3 0.3]);
xp=[-10:1:10]*dscan; % 21 points with defined spacing
xx=[-phonw:phonw/50:phonw]*1.75;
sigma=phonw/(sqrt(2)*2*log(2));
yp=phoni*exp(-1/2*(xp/sigma).^2)*1e9;
yy=phoni*exp(-1/2*(xx/sigma).^2)*1e9;
hp=plot(xp,yp,'ko',xx,yy,'r-');
axis([min(xx) max(xx) 0 1.5*phoni*1e9]);
title('Phonon Scan');
label=['DQE=[' num2str(p(35)) ',' ...
num2str(p(36)) ',' num2str(p(37)) ',' num2str(p(38)) ...
'] (meV Angs^-^1)' ];
xlabel(label);
ylabel('Normalized Intensity');
set(hp,'linewidth',2)
tt1=text(min(xx)+(max(xx)-min(xx))*.05,...
1.35*max(yy),['FWHM ' num2str(phonw) 'meVA^-^1']);

return
