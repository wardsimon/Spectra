function [y, name, pnames, pin]=pow1(x, p, flag)
% pow1      : Power law y=0 x<x0
% function [y, {name, pnames, pin}]=pow1(x,p, {flag})
%
% MFIT Power law fitting function y=0 x<x0
% p = [ Amp Xoffset Exponent Background ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description:  Power law y=0 x<x0

if nargin==2
	y= p(1)*(x>p(2)).*abs(x-p(2)).^p(3) + p(4);
else
	% ---------- Return initializtion data -----------------------
	y=[];
	name='Power law (y=0 x<x0)';
	pnames=str2mat('Amplitude','X.offset','Exponent','Background');
	if flag==1, pin=[0 0 0 0]; else pin = p; end
   if flag==2
        %-------------- Get starting guesses 
        % Get starting guess for x offset
        mf_msg('Click on x offset estimate (y=0 x<x0)');	    
        [xoff dummy]=ginput(1);

        % Get starting guess for background 
        mf_msg('Click on background estimate');	    
        [dummy bg]=ginput(1);

        % Get two points on the 'curve'
        mf_msg('Click on two points on curve');	    
        [xp yp]=ginput(2);			   

        % Work out starting guess for the exponent
        power=log((yp(1)-bg)/(yp(2)-bg))/log((xp(1)-xoff)/(xp(2)-xoff));

        % Work out starting guess for the amplitude
        amp=abs((yp(1)-bg)/((xp(1)-xoff)^power));
    
        pin=[amp xoff power bg];
    end
end

end
