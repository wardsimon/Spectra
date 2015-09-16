function [y, name, pnames, pin]=cusp(x, p, flag)
% cusp      : Power law cusp.
% [y, {name, pnames, pin}]=cusp(x,p, {flag})
%
%  MFIT Power law cusp fitting function.
%  p(1) = A1 (x<xc)
%  p(2) = A2 (x>=xc)
%  p(3) = xc
%  p(4) = pow1 (x<xc)
%  p(5) = pow2/pow1 (exp2: x>xc)
%  p(6) = bg

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  Power law cusp.

%	MZ 12.12.94 rev EF 27.06.97
%
if nargin==2
	y= p(1)*(x< p(3)).*abs(x-p(3)).^p(4) + ...
	   p(2)*(x>=p(3)).*abs(x-p(3)).^p(5) + p(6);
else
	% ---------- Return initialization data -----------------------
	y=[];
	name='Power law cusp';
	pnames=str2mat('Amplitude_1',...
						'Amplitude_2',...
						'X.offset',...
						'Exponent_1',...
						'Exponent_2',...
						'Background');
	if flag==1, pin=[1 1 0 1 1 1]; else pin = p; end
   if flag==2
        %-------------- Get starting guesses 
        % Get starting guess for x offset
        mf_msg('Click on x offset (threshold) estimate');	    
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
    
        pin=[amp amp xoff power power bg];
    end
end

