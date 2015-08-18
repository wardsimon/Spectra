function [mag,ph] = fftplot(Y,T)
% fftplot : FFT plot(yfft,{period})
%Syntax: [mag,ph] = fftplot(Y,T) compute and/or plot the magnitude
%                          in dB and phase in degrees
%                          of the FFT data in Y
%
%  T = sampling period, if this parameter is not given it is assumed
%      to be equal to one.
%
%  if called with both output arguments, i.e. / [m,p] = fftplot(Y)
%  the magnitude and phase are returned in m and p respectively
%  and no plot is generated
%
%  if called with one output argument, i.e.  / m = fftplot(Y)
%  the magnitude of Y in dB is returned and a plot of the magnitude
%  of Y is generated in the range 0 s w s pi
%

% Author:  I don't know
% Description: plots some fft data

lab = ['DIGITAL FREQUENCY / pi'
       '   FREQUENCY HERTZ    '];
ly = log(Y);
mag = 20 * real(ly) / log(10);
if nargout == 2
                ph = imag(ly);
elseif nargout <= 1
                n = max(size(Y))/2;
                f = (0:n-1)/n;
                if nargin == 2, f = 0.5*f/T; end
                plot(f,mag(1:n)),title('MAGNITUDE'),ylabel('dB'),
                xlabel(lab(nargin,:))
end
