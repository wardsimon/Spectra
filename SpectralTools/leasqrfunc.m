function y = leasqrfunc(x,p)
% leasqrfunc : this is a function of 'x' and 'p' parameters for leasqrexamp.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: leasqr example fit function

% example of function.
y=p(1)*exp(-p(2)*x);

% If you want to use some global variables, you must declare them in the body of that function.
% ex: 'global c;' with 'c' beeing global variable in Octave/MatLab.

%x=x(26:284);
%y=-p(2).*imag(1./((1+(i.*x*p(3))).^1)./x);
%y=convlv(y,c,2);
%y = y + p(1);

%end
