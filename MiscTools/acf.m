function [ak, lags] = acf(x, m, w)
% acf : Autocorrelation function
%Syntax: [ak, lags] = acf(x, m [,w]) compute autocorrelation function
%---        at m lags via Rader's method based on the FFT.
%       ==> works for complex-valued signals
%
%   Usage:   [ak, lags] = acf(x, m [,w])
%
%    x  : input signal
%    m  : number of lags (also = FFT length, if even)
%    w  : Hann window applied to acf, if there are 3 args
%    ak : autocorrelation function
%  lags : vector of lags [-m/2:m/2] for the output ACF
%

%---------------------------------------------------------------
% copyright 1994, by C.S. Burrus, J.H. McClellan, A.V. Oppenheim,
% T.W. Parks, R.W. Schafer, & H.W. Schussler.  For use with the book
% "Computer-Based Exercises for Signal Processing Using MATLAB"
% (Prentice-Hall, 1994).
%---------------------------------------------------------------

% Aufruf testen:
if nargin<2 | nargin>3 | nargout~=2,
   error('  Function call should be: [ak,la] = acf(x,m [,w]);');
   return
end;

x  = x(:).';
m = m - rem(m,2);   %<-- make m even, because m = FFT length
m2 = m/2;
Lx = length(x);
mu = (-1).^(0:m-1);
ak = zeros(1,m);
z  = zeros(1,m);
n1 = 1;
while( n1 <= Lx )
   n2 = min( n1+m2-1, Lx );
   zi = fft(x(n1:n2), m);
   ak = ak + zi.*conj( zi + mu.*z );
   z  = zi;
   n1 = n1 + m2;
end;
ak = ifft(ak)/Lx;
if( ~any(imag(x)) ),
   ak = real(ak);  end
ak = [conj(ak(m2+1:-1:2)) ak(1:m2+1)];
if nargin==3,
   ak = ak.*window(m+1,'hann')'; end;
lags = -m2:m2;
