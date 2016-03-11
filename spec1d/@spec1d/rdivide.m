function s_out = rdivide(s1,s2)
%
% function r = rdivide(s1,s2)
%
% @SPEC1D/times function to give the fraction of spectrum s1 with a value or spectrum s2.
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

if isa(s1,'spec1d')
    s = s1;
    if isa(s2,'spec1d')
        fac = [];
        if length(s2)>1 && length(s)==1
            temp = s2;
            s2 = s;
            s = temp;
            clear temp;
        elseif length(s2)>1 && length(s)>=1
            error('spec1d:divide:InputLengthError','Divide can not be used for two spectrum both with dimensions greater than 1')
        end
    else
        fac = s2;
        s2 = [];
        r = 0;
    end
else
    fac = s1;
    s = s2;
    s2 = [];
    r = 1;
end


if isempty(fac)
    for i = 1:length(s)
        % Divide spec1d objects
        if length(s(i).x)~=length(s2.x)
            warning('Objects are not the same length. Using interpolation for second object')
            s2 = interpolate(s2,s(i),'method','builtin');
        end
        
        r = s(i);
        r.y = r.y ./ s2.y;
        r.e = sqrt((s(i).e./s(i).y).^2 + (s2.e./s2.y).^2).*r.y;
        r.yfit = r.yfit ./ s2.yfit;
        s_out(i) = feval(class(r),r);
    end
else
    for i = 1:length(s)
        if r
            r   = s(i);
            r.y = fac./r.y;
            r.e = fac*r.e./s(i).y.^2;
            r.yfit = fac./r.yfit;
        else
            r   = s(i);
            r.y = r.y/fac;
            r.e = r.e/fac;
            r.yfit = r.yfit/fac;
        end
        s_out(i) = feval(class(r),r);
    end
end
