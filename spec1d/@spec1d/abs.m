function sout=abs(varargin)
%
% function r=abs(s1)
%
% @SPEC1D/abs function to give the absolute value of spectrum s1.
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

s1 = [varargin{:}];

for n = 1:length(s1)
    x = s1(n).x(:);
    y = s1(n).y(:);
    e = s1(n).e(:);
    yfit = s1(n).yfit(:);

    yabs = abs(y);

    if ~isempty(yfit)
        yfit = abs(yfit);
    end

    r = s1(n).copy;
    r.y = yabs;
    r.e = e;
    r.yfit = yfit;

    sout(n) = feval(class(r),r);
    
end