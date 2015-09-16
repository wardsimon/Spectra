function w = window(N,wt)
% window : Window functions (for spctrm)
%Syntax: w = window(N,wt)
%
%  generate a window function
%
%  N = length of desired window
%  wt = window type desired
%  'squa','rect' = rectangular        
%  'bart','tria' = triangular (Bartlett)            
%         'blac' = Blackman           'welc'  = Welch
%         'hann' = Hanning            'hamm'  = Hamming
%  w = row vector containing samples of the desired window

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  window functions

nn = N-1;
pn = 2*pi*(0:nn)/nn;

wt = wt + 0;
wt = wt(1:4);

if (strcmp(wt,'rect') | strcmp(wt,'squa'))
                        w = ones(1,N);
elseif (strcmp(wt,'tria') | strcmp(wt,'bart'))
                        m = nn/2;
                        w = (0:m)/m;
                        w = [w w(ceil(m):-1:1)];
elseif strcmp(wt,'hann')
                        w = 0.5*(1 - cos(pn));
elseif strcmp(wt,'hamm')
                        w = .54 - .46*cos(pn);
elseif strcmp(wt,'blac')
                        w = .42 -.5*cos(pn) + .08*cos(2*pn);
elseif strcmp(wt,'welc')
			m = nn/2;
			w = 1 - (((0:nn)-m)/m).^2 ;
else
                        disp('Incorrect Window type requested')
end
