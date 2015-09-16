% Signal processing functions
%
% simplex : Simplex routine (see fmins matlab function)
% derivative : Derivative approximation
% leasqr : Non Linear Least Square multivariable fit.
% leasqrexamp : leasqr and simplex fit example/test
% leasqrfunc : this is a function of 'x' and 'p' parameters for leasqrexamp.
% axisrescale : Rescales an axis by Lagrange interpolation
% convlv : Convolution (by FFTs), handles side effects.
% smooth : Data smoothing by Savitzky-Golay method
% dfdp : Partials by finite differencies.
% correl : Correlation
% u=isinft(z) inverse sine transform, assume that size(z)=[nx,ny]=[2^k-1,ny]
% u=sinft(z) sine transform, assume that size(z)=[nx,ny]=[2^k-1,ny]
% yt = sft(y) Slow Fourier transform function
% window : Window functions (for spctrm)
% spctrm : Power Spectrum Estimation by Welch's method
% findpeaks : Quick search of peaks in a signal.
% removepeak : Removes some peaks in a signal
% trapez : Trapezoidal integration
% interpsp : Lagrange 1D interpolation/smoothing.
