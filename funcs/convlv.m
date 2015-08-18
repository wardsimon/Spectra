function [y,D,R,lrm] = convlv (ds, r, op, N, last)
% convlv : Convolution (by FFTs), handles side effects.
%Syntax: [y,D,R] = convlv (data, respns, {op=2, N=[], mode='silent'})
%
% Convolves/deconvolves 'data' by 'respns' or autoconvolution.
% 'op' is 1 or for convolution, -1 for deconvolution.
% if op=2, side effects are reduced by signal extension.
% If n is supplied by user, a n-points FFT is used, else using powers of 2.
% FFT's are also returned.
% Note : 'respns' signal should better be centered.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: convolution (by FFTs), handles side effects.

% uses : fft, ifft,vect2row
% From : Numerical Recipes in C. p538
% Part of 'Spectral tools'. E.Farhi. 07/96 rev 09/97

if (nargin < 1)
  disp('y = convlv ( data, respns, {op=2, N=[], mode=''silent''} ) or convlv(data)');
  return
end

if (nargin == 1)
  r=ds;
end

if nargin < 3, op = []; end
if nargin < 4, N = []; end
if nargin < 5, last = []; end

if isempty(op)
  op = 2;
end

if isempty(last)
  last = 'silent';
end

if strcmp(last,'silent')
  tmp = 0;
else
  tmp = 1;
end

if length(ds) < length(r)
  tmp2 = ds;
  ds = r;
  r = tmp2;
end


  d = ds(:); d=d';
  r = r(:); r=r';


ld = length(d);
lr = length(r);
lrm = ceil(lr/2);

if min(size(ds) ~= 1) | min(size(r) ~= 1)
  disp(size(ds))
  disp(size(r))
  y = [];
  D = []; R=[];
  lrm = 0;
  error('Should be used with vectors.');
end

if (op == 2)  % extending data ends with data ends values
  d = [ d(1)*ones(1,lrm) d d(ld)*ones(1,lrm) ];
  ld = ld + 2*lrm;
  lrm = 2*lrm;
  if tmp
    disp('Extending signal for convolution.');
  end
end

% Use FFT with the smallest power of 2 which is >= length (x) +
% length (b) - 1 as number of points ...

if isempty(N)
    N = 2^(ceil (log (ld + lr - 1) / log(2)));
end

    D = fft(d, N);
    R = fft(r, N);
    if ( (op == 1) | (op == 2))
  y = ifft(D.*R);   % convolves
  if tmp
    disp('Convolution.');
  end

    elseif (op == -1)
  if tmp
    disp('Deconvolution.');
  end

  tmp2 = find(R == 0);
  if (~isempty(tmp2))
    warning('Deconvolving : zeros found in response-FFT.');
  end
  y = ifft(D./R);
    else
  error('No significance for operation parameter.');
    end

  y = y(lrm:(length(ds)+lrm-1));

% Final cleanups:  if both x and b are real respectively integer, y
% should also be

  if (~ (any (imag (d)) | any (imag (r))))
    y = real (y);
  end
  if (~ (any (d - round (d)) | any (r - round (r))))
    y = round (y);
  end

if (size(ds,2) == 1)
  y=y';
end






