function rc_spec
%
% MATLAB function to calculate the resolution function
% of a triple axis
%
% DAT 21.11.97 


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

if(TT(3)~=0) 
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

scan=DQ_cart/sqrt(sum(DQ_cart.*DQ_cart));
plane=(T(1:3,1:3)*[p(39) p(40) p(41)]')';
plane=[p(42)*plane/sqrt(sum(plane.*plane)) -1];
plane=plane/sqrt(sum(plane.*plane));

% ----- Calculate the phonon width.

[phoni,phonw]=rc_phon(1,NP,plane);
phonw=abs((phonw/(sum(scan.*plane))));
disp('Phonon width in (meV,1/Angs mixture)');
disp(phonw)


%----- This is to become a new routine ----------------------------
%----- Used Currently for Development -----------------------------
%------------------------------------------------------------------
%----- Reciprocal space plots -------------------------------------
%------------------------------------------------------------------


if ~isempty(findobj('Tag','Rescal: Spectrometer'))
   close(findobj('Tag','Rescal: Spectrometer'))
end

h=figure('Name','Rescal: Spectrometer',...
                'Tag','Rescal: Spectrometer',...
                'NumberTitle','off',...
                'Menubar','None',...
                'Visible','on','Position',[300 50 500 500]);
figure(h);
whitebg(h);
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
xlabel('V1 (1/Å)')
ylabel('V2 (1/Å)')
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
                         
                         
% --------------------------------------------------------------------------
disp('************ spectrometer configuration ************************')

disp(['ki=' num2str(ki) ' Angs-1  kf=' num2str(kf) ...
' Angs-1  q=' num2str(Qmag) ' Angs-1']);

disp(['The scattering angle phi (degrees) = ' num2str(ss*180*twotheta/pi)]);

disp(['The monochromator angle 2thetam (degrees) = ' ...
 num2str(sm*180*2*thetam/pi)]);
disp(['The analyser angle 2thetaa (degrees) = ' ...
 num2str(sa*180*2*thetaa/pi)]);
disp('The scattering sense is such that left(+ve) right(-ve)')
disp(['The angle psi between Kf and V1 is (degrees) ' ...
 num2str(180*(pi-phi-beta)/pi)]);
 disp('***************************************************************')
 
% Now draw up the spectrometer configuration

L=[1 1 1 1];
% _____________________Extra Parameters________________________________________
offset=42;

L0=p(20+offset);     % distance between source and monochromator (cm).
L1=p(21+offset);     % distance between monochromator and sample (cm).
L2=p(22+offset);     % distance between sample and analyser (cm).
L3=p(23+offset);     % distance between analyser and detector (cm).
L=[L0 L1 L2 L3]; 

alf0=p(11);     % horizontal pre-monochromator collimation.
alf1=p(12);     % horizontal pre-sample collimation.
alf2=p(13);     % horizontal post-sample collimation.
alf3=p(14);     % horizontal post-analyser collimation.   

                                                                   
DD=[L(1) 0; L(2)*sin(sm*2*thetam) L(2)*cos(sm*2*thetam)];
DD=[DD ; L(3)*sin(sm*2*thetam+ss*twotheta) ...
 L(3)*cos(sm*2*thetam+ss*twotheta)];
DD=[DD ; L(4)*sin(sm*2*thetam+ss*twotheta+sa*2*thetaa) ...
L(3)*cos(sm*2*thetam+ss*twotheta+sa*2*thetaa)];

angles=[0 ; sm*2*thetam ; sm*2*thetam+ss*twotheta ; ...
sm*2*thetam+ss*twotheta+sa*2*thetaa ]

DD= [ L'.*cos(angles) L'.*sin(angles)];

spec=[0 0 ;DD(1,:);DD(1,:)+DD(2,:);...
DD(1,:)+DD(2,:)+DD(3,:); DD(1,:)+DD(2,:)+DD(3,:)+DD(4,:)];
axes('Position',[0.65 0.15 0.3 0.3])
hp=plot(spec(:,1),spec(:,2),'ro');
set(hp,'markersize',15);
hSM=gf_arrow(spec(1:2,1),spec(1:2,2),[0 0 0],3,[num2str(alf0) '`']);
hMS=gf_arrow(spec(2:3,1),spec(2:3,2),[0 0 0],3,[num2str(alf1) '`']);
hSA=gf_arrow(spec(3:4,1),spec(3:4,2),[0 0 0],3,[num2str(alf2) '`']);
hAD=gf_arrow(spec(4:5,1),spec(4:5,2),[0 0 0],3,[num2str(alf3) '`']);
set(gca,'visible','off');
xmax=max(spec(:,1));
xmin=min(spec(:,1));
ymax=max(spec(:,2));
ymin=min(spec(:,2));
cen=[(xmax(1)+xmin(1))/2 (ymax(1)+ymin(1))/2];
dax=1.05*max([(xmax(1)-xmin(1)) (ymax(1)-ymin(1))]);
axis([cen(1)-dax(1)/2 cen(1)+dax(1)/2 cen(2)-dax(1)/2 cen(2)+dax(1)/2 ])
title('Spectrometer Configuration');


cen
dax

% -------------------- PLOT UP THE FLUX ----------------------------

T=300; % moderator temperature.
conv=1/11.8; % K->meV
kbT=conv*T; % characteristic temperature in meV.
E=[0:kbT/50:5*kbT];
flux=maxwell(E,T);
E=[E 5*kbT];
flux=[flux 0];

axes('Position',[0.15 0.6 0.3 0.3])
line(E,flux)
fill(E,flux,'r')
xlabel('Energy (meV)')
ylabel('Normalised Reactor Flux')
title('Flux Distribution')
axis([0 5*kbT 0 1.5])
Ei=1/f*ki^2;
fEi=maxwell(Ei,T);
hl=line([Ei Ei],[0 fEi]);
set(hl,'linewidth',2,'color',[0 0 0]);
text(kbT,1.2,['Ei= ' num2str(Ei) 'meV']);
str=[sprintf('Flux is %2.0f',fEi*100) '% of max.'];
text(kbT,1.4,str);


ha=axes('Position',[0.15 0.15 0.3 0.3])
axis([0 1 0 1])
set(ha,'visible','off');

% --------------------------------------------------------------------------
str=str2mat(['Spectrometer Angles:']);
str=str2mat(str,['ki=' num2str(ki) ' Å-1  kf=' num2str(kf) ' Å-1']);
str=str2mat(str,['q=' num2str(Qmag) ' Å-1']);
str=str2mat(str,['phi (degrees) = ' num2str(ss*180*twotheta/pi)]);
str=str2mat(str,['2thetam (degrees) = ' num2str(sm*180*2*thetam/pi)]);
str=str2mat(str,['2thetaa (degrees) = ' num2str(sa*180*2*thetaa/pi)]);

num=6;
xpos=zeros(1,num);
ypos=[0.55 0.4 0.3 0.2 0.1 0];
for i=1:num
  text(xpos(i),ypos(i),str(i,:));
end
