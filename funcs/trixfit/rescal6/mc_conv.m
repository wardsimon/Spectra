function mc_conv
%
% MATLAB routine to calculate and plot the convolution
% of the resolution ellipsoid with a dispersion.
% D.A.T. 22.10.95
% D.F.M. 14.11.95
%
% Called by: mc_win
% Calls  to: rc_re2rec, rc_cnmat, mc_monte, mc_sqw
%
% Units: f converts from energy units into k^2
%        f=0.4826 for meV and f=1.996854 for THz
%
%

f=0.4826;
mon_flag=0;      % 1 for monitor, 0 for time
method=get(findobj('tag','hrc_rescal_method'),'Userdata');

%----- Read rescal parameters from window

p=rc_savp('respars')';

%----- Read scan parameters from window

pscan=rc_savp('spars');

%----- Number of points in scan (NSC), and number of Monte Carlo steps (NMC)

NSC=pscan(9);
NMC=pscan(10);

%===================   Initialise everything   =====================

%----- Calculate Q2c matrix

[Q2c,Qmag]=rc_re2rc([p(19) p(20) p(21)], [p(22) p(23) p(24)], [p(31) p(32) p(33)]);

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

S=U*Q2c;     % This is used to bring the CN matrix into a defined frame.

%============ End Initialisation ===================================

intensity=[];
x_axis=[];

%----- Calculate step size

dscan=(pscan(5:8)-pscan(1:4))/(NSC-1);

%----- Find scan variable; choose first non-zero element of dscan

scan_var=find(dscan);

%----- Loop over number of points in scan

for j=1:NSC

%----- Calculate Q,E vector

   p(31:34)=pscan(1:4)+(j-1)*dscan(1:4);

%----- Q vector in cartesian coordinates

   Qcart=Q2c*p(31:33)';
   Qmag=sqrt(sum(Qcart.*Qcart));
 
%----- Calculate resolution matrix in Q frame. NP is resolution matrix in Qx, Qy & Qz

   [R0,NP,vi,vf,Error]=feval(method,f,Qmag,p,mon_flag);
   R0_corrected=R0/(sqrt(det(NP))/(2*pi)^2); 	% corrected resolution
                                                % volume as Monte Carlo
						% integral is over a normalised
						% ellipsoid.
%----- Work out angle of Q wrt to V1, V2

   TT=S*[p(31) p(32) p(33)]';
   cos_theta=TT(1)/sqrt(sum(TT.*TT));
   sin_theta=TT(2)/sqrt(sum(TT.*TT));

%----- Rotation matrix from Q to V1,V2,V3

   R=[cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1];

   T=zeros(4,4);
   T(4,4)=1;
   T(1:3,1:3)=R*S;

%----- Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV]) 

   M=T'*NP*T; 

%----- Perform Monte-Carlo integration

   if Error == 0

%   [X]=mc_monte(M,NMC);
%   s=mc_sqw(pscan(11:22),X(1,:)+p(31),X(2,:)+p(32),X(3,:)+p(33),X(4,:)+p(34));
%   x_axis=[x_axis p(31+(scan_var(1)-1))];
%   intensity=[intensity R0_corrected*sum(s)/NMC];

%----- Alternative Monte-Carlo code for mex file 

[V,E]=eig(M);
fwhm=zeros(1,4);
fwhm(1)=1/sqrt(E(1,1));
fwhm(2)=1/sqrt(E(2,2));
fwhm(3)=1/sqrt(E(3,3));
fwhm(4)=1/sqrt(E(4,4));
b_mat=reshape(inv((V'))',1,16);
x_axis=[x_axis p(31+(scan_var(1)-1))];
[smcmex]=mc_tasfit(NMC,p(31:34),pscan(11:22),fwhm,b_mat);
intensity=[intensity R0_corrected*smcmex];

%----- End of alternative routine.
   
      set(findobj('Tag','MC Message Window'),'String',[' Point: ' num2str(j)]);
      drawnow

   else

      set(findobj('Tag','MC Message Window'),'String',...
         ['Scattering triangle will not close for point: ' num2str(j)]);      
      drawnow

   end	

end

% set(findobj('Tag','MC Message Window'),'String',[' Done(' num2str(toc) 's)']);

%----- Create figure window and plot results

if ~isempty(findobj('Tag','MC Plot'))
   close(findobj('Tag','MC Plot'))
end

hmc_plot=figure('Name','Rescal: Simulation Plot',...
                'Tag','Rescal: Simulation Plot',...
                'Menubar','None',...
                'Numbertitle','off');

hmc_plot_print=uimenu(hmc_plot,'Label','Print');
uimenu(hmc_plot_print,'Label','Print Figure','Callback','mc_prt');

plot(x_axis,intensity,'yo')
hold on
norm=max(intensity);

%----- Make labels

[n,m]=size(scan_var);

if m ==1					% Single variable scan

   x_matlabs=str2mat('H','K','L','Energy');
   x_lab=x_matlabs(scan_var(1),:);

elseif m==2 & scan_var(2) ~= 4  % Two Q components

   x_matlabs=str2mat(' ',' ','H (HK scan)','H (HL Scan)',' K (KL Scan)');
   x_lab=x_matlabs(scan_var(1)+scan_var(2),:);  

elseif m==2 & scan_var(2) == 4 % One Q and Energy

   x_matlabs=str2mat('H (Mixed Q+E Scan)','K (Mixed Q+E scan)','L (Mixed Q+E scan)');
   x_lab=x_matlabs(scan_var(1),:);

elseif m>=3

   x_lab='This is a complicated Q+E Scan';
 
end

%----- Label axes

xlabel(x_lab)
ylabel('Intensity')
%output=[x_axis' 100*intensity' sqrt(intensity)'];
%save test.out output -ascii 

%----- Plot unconvoluted scattering function

intensity=[];
x_axis=[];

%----- Calculate step size

dscan=(pscan(5:8)-pscan(1:4))/(NSC-1);

for j=1:NSC

%----- Calculate Q,E vector

   p(31:34)=pscan(1:4)+(j-1)*dscan(1:4);

%----- Q vector in cartesian coordinates

   Qcart=Q2c*p(31:33)';
   Qmag=sqrt(sum(Qcart.*Qcart));
 
%----- Calculate resolution matrix in Q frame. NP is resolution matrix in Qx, Qy & Qz

   [R0,NP,vi,vf,Error]=feval(method,f,Qmag,p,mon_flag);


%----- Calculate unconvoluted scattering function

   if Error == 0

      s=mc_sqw(pscan(11:22),p(31),p(32),p(33),p(34));
      intensity=[intensity s];
      x_axis=[x_axis p(31+(scan_var(1)-1))];   
      set(findobj('Tag','MC Message Window'),'String',[' Point: ' num2str(j)]);
      drawnow

   else
      set(findobj('Tag','MC Message Window'),'String',...
         ['Scattering triangle will not close for point: ' num2str(j)]);      
      drawnow
   end	

end

% set(findobj('Tag','MC Message Window'),'String',[' Done(' num2str(toc) 's)']);

%----- Overlay unconvoluted scattering function

%kill=find(isnan(intensity));
kill=find(isnan(intensity) | ~isfinite(intensity));

intensity(kill)=[];
x_axis(kill)=[];
norm=norm/max(intensity);
plot(x_axis,norm*intensity,'r--')

%----- Add legend

legend('Covoluted S(Q,w)','Un-Convoluted S(Q,w)')

%----- Save Convoluted data to file

output=[x_axis' norm*intensity'];

