function [y, name, pnames, pin]=trix2(x,p, flag)
% trix2     : 4D convoluted fit function for TAS v2
% MFIT function [y, name, pnames, pin]=trix2(x,p, flag)
% for triple-axis data. 
% This function uses a quadratic phonon dispersion curve
% and a Damped Harmonic Oscillator. 
% (Fak and Dorner, Physica B 234-236 (1997) 1107)
% The whole stuff is convoluted with the 4D resolution function.
%
% Calls to: rescal, trixsqw, rc_re2rc, rc_cnmat, rc_savp, mc_monte
%	
% Des McMorrow 25.1.96
% E.Farhi 09.97

global GLOBAL_LENX GLOBAL_PRES GLOBAL_PSCAN GLOBAL_METHOD
global GLOBAL_B_MAT GLOBAL_FWHM GLOBAL_R0

if nargin==2;

%----- Get method from Parameter window
   method = findobj('tag','hrc_rescal_method');
   if ~isempty(method)
	method=get(method,'Userdata');
   else
	return
   end

%----- Pre-allocate matrices in for loops

   hmf_data=findobj('tag','mf_DataWindow');
   data=get(hmf_data,'userdata');
   index=data(:,4);
   xdata=data(:,1); % get limits and whole data
   xstart=xdata(1);
   xend=xdata(length(xdata));

   lx = length(x);

   y=zeros(size(x));

%----- Read rescal parameters from window

   pres=rc_savp('respars')'; % Rescal:Parameters:sections Spectro + Lattice in row

%----- Read scan parameters from window trixpar window

   pscan=rc_savp('tpars')'; % Rescal:Parameters:section Scan with Mfit in row

   dscan=(pscan(5:8)-pscan(1:4));

%----- Check to see if the scan parameters or spectrometer parameters
%----- have changed - if not use the previously calculated values of
%----- the resolution ellipsoid, otherwise recalculate these quantities. 
compared = 0;
if ~isempty(GLOBAL_LENX)
	compared = compared + sum(lx==GLOBAL_LENX);
end
if ~isempty(GLOBAL_PRES)
	compared = compared + sum(pres==GLOBAL_PRES);
end
if ~isempty(GLOBAL_PSCAN)
	compared = compared + sum(pscan==GLOBAL_PSCAN);
end
if ~isempty(GLOBAL_METHOD)
	compared = compared + sum(method==GLOBAL_METHOD);
end

totlen=(length(lx)+length(pres)+length(pscan)+length(method));

if (compared~=totlen)

%----- Units: f converts from energy units into k^2 (0.4826 for meV, 1.996854 for THz)
%      At the moment, only works for meV!

   f=0.4826;
   mon_flag=p(16);               %  mon_flag: 1=monitor, 0=time

   GLOBAL_LENX=lx;
   GLOBAL_PRES=pres;
   GLOBAL_PSCAN=pscan;
   GLOBAL_METHOD=method;

   disp('Rescal : Compute 4D-TAS Resolution function')

%----- Calculate Q2c matrix with Rescal

   Q2c = rc_re2rc( pres(19:21), pres(22:24), pres(31:33) );

%----- Now work out transformations

   A1=pres(25:27)'; % axis 'A' in Lattice in cols
   A2=pres(28:30)'; % axis 'B' 

   V1=Q2c*A1; % in cols
   V2=Q2c*A2;

%----- Form unit vectors V1, V2, V3 in scattering plane

   V3=cross(V1,V2);
   V2=cross(V3,V1);
   V3=V3/norm(V3);
   V2=V2/norm(V2);
   V1=V1/norm(V1); % in cols

   U=[V1 V2 V3]';

%----- S transformation matrix from (h,k,l) to V1,V2,V3

   S=U*Q2c;     % This is used to bring the CN matrix into a defined frame.

%----- Calculate step size

%   dscan=(pscan(5:8)-pscan(1:4)); % hkle step

   GLOBAL_B_MAT = zeros(lx,16);
   GLOBAL_FWHM = zeros(lx,4);
   GLOBAL_R0 = zeros(lx,1);
                                                
   for j=1:lx

      	pres(31:34)=pscan(1:4)+(x(j)-xstart)/(xend-xstart)*dscan(1:4); % current hkle 

%----- Q vector in cartesian coordinates

      	Qcart=Q2c*pres(31:33)';
      	Qmag=norm(Qcart);
 
%----- Calculate resolution matrix in Q frame. NP is resolution matrix in Qx, Qy & Qz

      	[R0,NP,vi,vf,Error]=feval(method,f,Qmag,pres,mon_flag);
      	R0_corrected=R0/(sqrt(det(NP))/(2*pi)^2); % corrected resolution
						% volume as Monte Carlo
						% integral is over a normalised
						% ellipsoid.

%----- Work out angle of Q wrt to V1, V2

      	TT=S*pres(31:33)';
      	cos_theta=TT(1)/norm(TT);
      	sin_theta=TT(2)/norm(TT);

%----- Rotation matrix from Q to V1,V2,V3

      	R=[cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1];

      	T=zeros(4);
      	T(4,4)=1;
   	T(1:3,1:3)=R*S;

%----- Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV]) 

      	M=T'*NP*T; 

	[V,E]=eig(M);
	fwhm=1./sqrt(diag(E))';
	b_mat=reshape(inv((V'))',1,16);

	GLOBAL_B_MAT(j,:)=b_mat; % lx*16 axes
	GLOBAL_FWHM(j,:)=fwhm; % lx*4 widths
	GLOBAL_R0(j)=R0_corrected; % lx*1 MC volume

    end % for

    fvcpx = sum(abs(imag(GLOBAL_FWHM')));
    if any(fvcpx)
	disp('Warning : Physical instrumental limits reached.');
	fvrel = find(~fvcpx);
	for j=find(fvcpx)
		[dummy,k]=min(abs(fvrel-j));	% closest element real ellipsoid.
		GLOBAL_FWHM(j,:) = GLOBAL_FWHM(k,:);
		GLOBAL_R0(j) = GLOBAL_R0(k);
		GLOBAL_B_MAT(j,:) = GLOBAL_B_MAT(k,:);
	end
	fprintf(1,'          Fixing %i ellipsoids on scan.\n',length(find(fvcpx)));
	GLOBAL_FWHM = real(GLOBAL_FWHM);
	GLOBAL_R0 = real(GLOBAL_R0);
	GLOBAL_B_MAT = real(GLOBAL_B_MAT);
    end

end % if totlen

%----- Do scan  ...
            
    NMC=p(15);
    V0 = [ p(1) p(2) p(3) ]; % valley main axis (V0)
    V0 = V0/norm(V0);
    Aaxe = pres(25:27); % lattice 'A' axis
    Baxe = pres(28:30); % lattice 'B' axis
    Caxe = cross(Aaxe,Baxe); % lattice 'C' axis
    V1   = cross(Caxe,V0);
    V1   = V1/norm(V1);
    V2   = cross(V0,V1);

    partemp = [ V0 p(4) p(5) p(6) p(7) p(8) p(9) p(10) p(11) p(12) p(13) V1 V2 ];

    for j=1:lx
	% position in scan HKLE
	pres(31:34) = pscan(1:4)+(x(j)-xstart)/(xend-xstart)*dscan(1:4);

	% call MC routine
	smcmex = mcint2(NMC,pres(31:34),...
			partemp,GLOBAL_FWHM(j,1:4),GLOBAL_B_MAT(j,1:16));
	y(j)=GLOBAL_R0(j)*smcmex;   % correct for MC volume
    end

   y = p(11)*y + p(14); % add background

else

%----- Mfit parameter window

	y=[];
	name='TRIX 2';
	mf_msg('Click on Guess button to have some more info...');
	pnames=str2mat('DC:Xdir','DC:Ydir','DC:Zdir','DC:Lx','DC:Qy','DC:Qz',...
			'Ph:Xo','Ph:Yo','Ph:Zo','Ph:Energy',...
			'DHO:Amplitude','DHO:Damping',...
			'DHO:Temp','Background','NMC','Mon Flag');

	if flag==1, pin=[ 1 0 0 20 0 0 0 0.05 0 1 10000 0.05 10 1 100 1 ]; else pin = p; end
	if flag==2
		mf_msg('Sorry, this function does not support guess.')
		disp('TRIX v2 4D-TAS function')
		disp('Dispersion curve (DC) is :');
		disp('Linear on axis (XYZ)dir, quadratic on two other base vectors');
		disp('Negative quadratic curvatures for Qy and Qz indicate inverted concavity.')
		disp('The XYZdir is usually the scan direction along phonon propagation vector.')
		disp('For q = (Xo,Yo,Zo) E = Ph:Energy');
		pres=rc_savp('respars')';
		V0 = [ p(1) p(2) p(3) ]; % valley main axis (V0)
		V0 = V0/norm(V0);
 		Aaxe = pres(25:27); % lattice 'A' axis
  		Baxe = pres(28:30); % lattice 'B' axis
   		Caxe = cross(Aaxe,Baxe); % lattice 'C' axis
    		V1   = cross(Caxe,V0);
    		V1   = V1/norm(V1);
    		V2   = cross(V0,V1);
		disp('Valley base is (in rows):')
		disp(V0)
		disp(V1)
		disp(V2)
		pin = p;
	end

%----- TRIX parameter window for resolution and scan parameters

	rescal('trixmode');


end
