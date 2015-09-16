function [y, name, pnames, pin]=rescon(x,p, flag)

% rescon    :  3D resolution convolution of Lorentzian+Lorentzian squared
%        scattering function with triangular out-of-plane and Gaussian
%        in-plane resolution function. See m-file comments for parameter syntax
%        MZ 30.5.95

% Note on speed: gains of up to ~50% are possible with small modifications
% at a cost of clarity; for instance, it is only necessary to evaluate the
% gauss basepoints and resolution function if the res parameters change. However,
% large gains (factors of 3) would mean major recoding.

if nargin==2;
   nres=3;                 % Number of resn widths in integration
   
%====== Assign model parameters from p vector============================

%------ Resolution and integration parameters----------------------------
   k=p(1);                 % Triangular vertical resolution half-width half max  
   yresw=p(2);             % In-plane resolution width (gaussian sd, perp. to scan dir)   
   zresw=p(3);             % In-plane resolution width (gaussian sd, par. to scan dir)   
   phi=p(4);               % Angle between scan dir and zres
   ngss=p(5);                % Order of gauss rule

%------ Lorentzian scattering function parameters------------------------
	lamp=p(6);              % Lorentzian amplitude
	lwid=p(7);              % Lorentzian width
   lyoff=0;                % Lorentzian centre: y...
   lzoff=p(8);             % ...and z

%------ Lorentzian scattering function parameters------------------------
	l2amp=p(9);             % Lorentzian squared amplitude
	l2wid=p(10);            % Lorentzian squared width
   l2yoff=0;               % Lorentzian squared centre
	l2zoff=p(8)+p(11);           % ...and z

%------ Other model parameters ------------------------------------------
	bg=p(12);               % Background
	bgslope=p(13);          % Background slope
   temp=p(13);             % Temperature
	qz=x(:);                % Rename dep. variable from x to qz
	qy=zeros(size(qz));     % Keep qy for future extension: =0 here
   y=zeros(size(qz));      % preassign output vector
   yl=y;
   yl2=y;
   
%======= Calculate Gauss basepoints and weights ==============
   [bpy bpz wyz]=grule2d(ngss, ngss);                    % Calculate Gauss bp's and weights
   bpy=bpy*yresw*nres;                                   % Scale to cover rectangular region 
   bpz=bpz*zresw*nres;                                   % nres res widths      

   resfun=nres^2/pi/2*wyz...                             % norm. res. func * weights
          .*exp(-0.5*((bpy/yresw).^2+(bpz/zresw).^2));
   CosPhi=cos(phi*pi/180);                               % rotation of res. ellipse      
   SinPhi=sin(phi*pi/180);

%======= Do convolution ==================================================
   for i=1:length(qz);

%------- Lorentzian convolution ------------------------------------------
      qpy=qy(i)-lyoff;                                      % Transform integration
      qpz=qz(i)-lzoff;                                      % variables
      qppy= qpy*CosPhi+qpz*SinPhi;
      qppz=-qpy*SinPhi+qpz*CosPhi;

      t=k./sqrt(lwid^2+(bpy-qppy).^2+(bpz-qppz).^2);        % t=k/lambda
      lorfun=2*t.*atan(t)-log(1+t.^2);
      yl(i)=sum(sum(resfun.*lorfun));   

%------- Lorentzian squared convolution ----------------------------------
      qpy=qy(i)-l2yoff;                                      % Transform integration
      qpz=qz(i)-l2zoff;                                      % variables
      qppy= qpy*CosPhi+qpz*SinPhi;
      qppz=-qpy*SinPhi+qpz*CosPhi;

      t=k./sqrt(l2wid^2+(bpy-qppy).^2+(bpz-qppz).^2);        % t=k/delta
      lor2fun=t.^3.*atan(t);
      yl2(i)=sum(sum(resfun.*lor2fun));   

   end
   y=yl*lamp*(lwid/k)^2 + yl2*l2amp*(l2wid/k)^4 + bg + bgslope*(qz-mean(qz));

else

%======== Either initialization or guess requested ===========================
	y=[];
	name='Lor+Lor^2 with res. conv.';
	pnames=str2mat('k','yresw','zresw','Phi','Gauss.order','A','gamma','zlor');
	pnames=str2mat(pnames,'B','alpha','zlor2.offset','Background','Bg.slope',...
	'Temperature');
	if flag==1, pin=[1 1 1 0 20 0 1 0 0 1 0 0 0 0]; else pin = p; end

%-------- Interactively get starting parameters ------------------------------
	if flag==2

		mf_msg('Click on background');
		[x bg]=ginput(1);

		mf_msg('Click on lorentzian peak');
		[lcen lamp]=ginput(1);
		mf_msg('Click on lorentzian width');
		[lwidth y]=ginput(1);
		lwidth=abs(lwidth-lcen);
		lamp=lamp-bg;

		mf_msg('Click on lorentzian squared peak');
		[l2cen l2amp]=ginput(1);
		mf_msg('Click on lorentzian squared width');
		[l2width y]=ginput(1);
		l2width=abs(l2width-l2cen);
		l2amp=l2amp-bg;

		pin=[0.069 0.00075 0.00044 0 20 lamp lwidth lcen ...
			  l2amp l2width lcen-l2cen bg 0 0];
	end
end
