function yt=sft(y)
% yt = sft(y) Slow Fourier transform function
%
% Input
%    y - Time series
% Output
%    yt - Discrete Fourier transform
N = length(y);   % Length of the time series
twopiN = -2*pi*sqrt(-1)/N;
for k=0:N-1
 temp = exp(twopiN*(0:N-1)*k);
 yt(k+1) = sum(y .* temp);
end
%return;
