function rvec = times(s1,s2)
%
% function r = times(s1,s2)
%
% @SPEC1D/times function to give the product of spectrum s1 with a value or spectrum s2.
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

if isa(s1,'spec1d')
    s(1) = s1;
    if isa(s2,'spec1d')
        s(2) = s2;
        fac = [];
    else
        fac = s2;
    end
else
    fac = s1;
    s(1) = s2;
end

if isempty(fac)
    % Multiply spec1d objects
    if length(s(1).x)~=length(s(2).x)
        warning('Objects are not the same length. Using interpolation for second object')
        s(2) = interpolate(s(2),s(1));
    end
    
    rvec = s(1);
    rvec.y = rvec.y .* s(2).y;
    rvec.e = sqrt((s(1).e./s(1).y).^2 + (s(2).e./s(2).y).^2).*abs(rvec.y);
    rvec.yfit = rvec.yfit .* s(2).yfit;
else
    rvec   = s(1);
    rvec.y = fac*rvec.y;
    rvec.e = abs(fac)*rvec.e;
    rvec.yfit = fac*rvec.yfit;
end
    