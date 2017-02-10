function s_out = abs(varargin)
% Function to give the absolute value of spectrum.
%
% s_out = abs(s_in)
%
% Input:
%
% s_in  A single or array of spec1d objects
%
% This function returns a spec1d object of which the absolute value of the 
% signal and y-fit values are returned (if available) 
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

    s_out(n) = feval(class(r),r);    
end