function [y, name, pnames, pin]=refl(x,p, flag)
% refl      : reflectivity function
%  MFIT function [y, name, pnames, pin]=refl(x,p, flag)
%	Reflectivity data, N layers, using the full Dynamic theory
%	and the equations of Parratt.  Note that:-
% 	i)   x is Scattering vector, given in reciprocal Angstroms
%	ii)  the electron density is given in e/Angstroms^3
%	iii) The linear absorption coefficient in reciprocal Angstroms
%		(should be of order 1E-6)
%	AW  12.4.95
% 
%  Incoporated resolution in the form of a Gaussian
%     	AW  28.11.96

global nint res flag_aw
c = 8*pi;
ro = 2.82e-5;
lam = 1.54051;

if nargin==2;
%
%	Initialise everything
%

   F=[];
   R=[];
   
   ngss=p(length(p));
   res=p(length(p)-1);

   if ngss==0
      x=x(:);
   else
      Q=x(:);		% Save this for later
      x=[x(1)-3*res/2.355:min(x(2)-x(1),res/ngss):x(length(x))+6*res/2.355]';
   end;

   f=zeros(length(x),p(1)+1);
   a=zeros(length(x),p(1));
   rms=zeros(length(x),p(1));

   f(:,1) = x;
   a(:,1)=ones(length(x),1);

%
%	Set up f(n) (Fresnel coefficients for the various interfaces)
%
   for ii = 2:1:p(1)
     f(:,ii) = sqrt(x.^2-2*c*ro*p(4*ii-6)-c*j*p(4*ii-5)/lam);
     rms(:,ii-1) = exp(-(x*p(4*ii-4)).^2); % The Roughness coefficient
     if p(1) > 1
       a(:,ii) = exp(-j*f(:,ii)*p(4*ii-3)); % Calculate a(n) (Amplitude factor for half perpendicular depth)
     end;  
   end;

   f(:,p(1)+1) = sqrt(x.^2-2*c*ro*p(length(p)-7)-c*j*p(length(p)-6)/lam);
   rms(:,p(1)) = exp(-(x*p(length(p)-5)).^2);

%
%	Calculate reflectivity via recursion relation
%
   R = 0;
   for i = p(1)+1:-1:2
     F = rms(:,i-1).*(f(:,i-1)-f(:,i))./(f(:,i-1)+f(:,i));
%	  F = (f(:,i-1)-f(:,i))./(f(:,i-1)+f(:,i));
     R = a(:,i-1).*(R+F)./(R.*F+1);
   end;
%
%	Calculate geometric factors
%
   L = p(length(p)-3);
   sigma = 0.274;
   tm = 0.697/(sqrt(2)*sigma);
   lim = L*lam*x/(8*pi*sqrt(2)*sigma);
   geom = erf(lim)/erf(tm);
%   geom = L*x*lam/(4*pi*0.488);
   temp = find(geom > 1);
   geom(temp) = ones(size(temp));

%
%	And finally the function...
%
   func = p(length(p)-2)+p(4*p(1)+1)*geom.*abs(R).^2;

%	Now for resolution
   if ngss == 0
      y=func;
   else
      [bp,wf]=grule(ngss);
      bp=bp*3*res/2.355;
      for i = 1:length(Q)
         vbp = Q(i)+bp;
         I = interp1(x,func,vbp);
         y(i) = sum(I.*wf');
      end;
      y=y';
   end;

else
   y=[];
   name='Reflectivity';
   if size(flag_aw,2) < 4
      flag_aw = [flag_aw zeros(1,4)];
   end;
   if flag_aw(1:4) ~= 'refl'
      nint=input('Number of interfaces:');
      res=input('Enter resolution width (‰-1): ');
      flag_aw = str2mat('refl');
   end;
   pnames = str2mat('N-interfaces');
   for i = 1:nint-1
     pnames=str2mat(pnames,...
 	    sprintf('e-density_%d',i),...
	    sprintf('Absorption_%d',i),...
	    sprintf('Roughness_%d',i),...
	    sprintf('Thickness_%d',i) );
   end;
   pnames=str2mat(pnames,...
          sprintf('e-density_s'),...
  	  sprintf('absorption_s'),...
	  sprintf('Roughness_s'),...
	  sprintf('Incident Intensity'),...
	  sprintf('Sample width(mm)'),...
	  sprintf('Background'),...
          sprintf('Res width'),...
          sprintf('Res order') );
   pin = [nint;zeros(4*nint+2,1);res;0];
end



