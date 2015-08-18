function [y, name, pnames, pin]=voigt(x,p, flag)
% voigt     : Voigt
% function [y, {name, pnames, pin}]=voigt(x,p, {flag})
%
% MFIT Voigt fitting function
% p = [ Amplitude Centre Gauss_width Lorz_width Background ]

% Author:  MZ <mzinkin@sghms.ac.uk> adapted from DFM
% Description:  Voigt
nv = (length(p)-3)/4;

if nargin==2;
    y=zeros(size(x));
    for i=1:nv
        off=(i-1)*4;
        N = 16;
        b = -sqrt(log(2))/p(3+off);
        a = b*p(4+off);
        b = b*2*1i;
        z = a + b*(x-p(2+off));
        
        M=2*N; M2=2*M; k=[-M+1:1:M-1]';
        L=sqrt(N/sqrt(2));
        tt=(L*tan(k*pi/M2)).^2;
        f=[0; exp(-tt).*(L^2+tt)];
        a=real(fft(fftshift(f)))/M2;
        a=flipud(a(2:N+1));
        l=L-z;
        Z=(L+z)./l;
        pp=polyval(a,Z);
        y=y+p(1+off)*real(2*pp ./l.^2+(1/sqrt(pi))*ones(size(z)) ./l);
    end
	y=y+(p(5+off)+p(6+off)*x + p(7+off)*x.^2);
else
   name=' N Voigt';
   pin = p;
   y=zeros(size(x));
   if flag==2
   	ngss=0;

      	mf_msg('Click on background');    
	[x bg]=ginput(1);

	but=1;
	pin=[];
		
	while but==1
		mf_msg(sprintf('Click on peak %d (right button to end)',ngss+1));
		[cen amp but]=ginput(1);
		if but==1
			mf_msg(sprintf('Click on width %d (right button to end)',ngss+1));
			[width y but]=ginput(1);
			width=abs(width-cen);
			amp=amp-bg;
			pin=[pin amp cen width];
			ngss=ngss+1;
		end
	end
	pin=[ pin bg];
   end

   y=[];
   pnames=str2mat('Int_1','Centre_1','FWHM_1');
   if ngss>0 & any(p)
	name=sprintf('%d Gaussians',ngss);
   else
	name='n gauss  : clik on Guess button to set n.';
	if flag==1, pin=[ 0 0 1 0 ]; else pin = p; end
   end

   for i=2:ngss
	pnames=str2mat(pnames,...
               sprintf('Int_%d',i),...
               sprintf('Centre_%d',i),...
               sprintf('FWHM_%d ',i));
   end;
   pnames=str2mat(pnames,'Background');
	
end
