function  [P, c] = spctrm(xs, N, M, wintype,last)
% spctrm : Power Spectrum Estimation by Welch's method
%Syntax: [P, c] = spctrm(xs, N, M, wintype,mode='silent')   
%-----
%   Usage:      [P,c] = spctrm(x {N,M,wintype})
%
%       P : power spectrum by Welch's method
%       c : correlation function = inverse of Welch power spectrum
%       x : input signal vector
%       N : FFT length
%       M : window length
%   wintype = window type, 'hamming' by default (see window.m)

%---------------------------------------------------------------
% copyright 1994, by C.S. Burrus, J.H. McClellan, A.V. Oppenheim,
% T.W. Parks, R.W. Schafer, & H.W. Schussler.  For use with the book
% 'Computer-Based Exercises for Signal Processing Using MATLAB'
% (Prentice-Hall, 1994).
%---------------------------------------------------------------

if ~exist('last')
	last = 'silent';
end
if isempty(last)
	last = 'silent';
end

if strcmp(last,'silent')
	tmp = 0;
else
	tmp = 1;
end

x = xs(:);
if ~exist('M') 
	M = [];
end
if isempty(M)
	M =length(x);
end
if ~exist('N') 
	N = [];
end
if isempty(N)
	N = 2^(ceil(log (length(x) -1) / log(2)));
end
if ~exist('wintype') 
	wintype = [];
end
if isempty(wintype)
    wintype = 'hamm';
end

if (tmp>0)
	fprintf(1,'Power Spectrum computation with %i data points, "%s" window (%i points).\n', N, wintype, M);
end
w = window(M,wintype);
K = fix((length(x)-M/2)/(M/2));
KLU = K*sum(w.^2);
P = zeros(1,N);
for i=1:K
    ibeg = (M/2)*(i-1) + 1;
    X = fft(x(ibeg:ibeg+M-1).*w, N);
    P = P + abs(X).^2;
end
P = P/KLU;
c = ifft(P);
c = real(c(1:M/2));
X = P(1:N/2+1);
if( tmp > 1 )
  axis([0 1 0 max(P)])
  title(['Power Spectrum: Window length=',...
      num2str(M), '  Number sections=', num2str(K)])
  xlabel('FREQUENCY  (omega/pi)')
  plot((0:N/2)/(N/2),X)
  pause
  title(['Log Power Spectrum: Window length=',...
      num2str(M), '  Number sections=', num2str(K)])
  xlabel('FREQUENCY  (omega/pi)')
  plot((0:N/2)/(N/2),10*log10(X))
  pause
  title(['Correlation Function: Window length=',...
      num2str(M), '  Number sections=', num2str(K)])
  xlabel('LAG INDEX')
  plot(0:M/2-1, c)
end

if (size(xs,1)==1)
	P = P';
	c = c';
end
