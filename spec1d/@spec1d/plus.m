function s_out = plus(s1,s2)
%
% function r = plus(s1,s2)
%
% @SPEC1D/plus function to give the value of spectrum s1 plus a value or spectrum s2.
%
% 1. r = s1+s2     , where s1 and s2 are spectra
% 2. r = s1+a      , where a is to be added to the y axis
% 3. r = s1+[xc,yc], where xc (yc) is added the x (y) axis
%
% If spectra are not the same length, s2 is interpolated to
% the x values of s1.
%
% Plus should not be used to combine counts taken with different
% monitors - use combine(s1,s2)
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

if isa(s1,'spec1d')
    s = s1;
    if isa(s2,'spec1d')
        val = [];
        if length(s2)>1 && length(s)==1
            temp = s2;
            s2 = s;
            s = temp;
            clear temp;
        elseif length(s2)>1 && length(s)>=1
            error('spec1d:times:InputLengthError','Times can not be used for two spectrum both with dimensions greater than 1')
        end
    else
        val = s2;
        s2 = [];
    end
else
    val = s1;
    s = s2;
    s2 = [];
end

for i = 1:length(s)
    if isempty(val)
        if length(s(i).x)~=length(s2.x)
            warning('spec1d:plus:UnevenSpectra','Number of points in s1 and s2 are not equal. Using Interpolation')
            s2 = interpolate(s2,s(i),'method','builtin');
        end
        
        r = s(i);
        r.y = r.y + s2.y;
        r.e = sqrt(r.e.^2 + s2.e.^2);
        if ~all([isempty(r.yfit) isempty(s2.yfit)])
            if isempty(r.yfit) || isempty(s2.yfit)
                warning('spec1d:plus:EmptyYFit','The y-fit of spectrum s1 or s2 is empty. No summation has occured on the fit channel')
                if ~isempty(s2.yfit)
                    r.yfit = s2.yfit;
                end
            else
                r.yfit = r.yfit + s2.yfit;
            end
        end
    else
        r   = s(i);
        if length(val)==1
            r.y = val + r.y;
            if ~isempty(r.yfit)
                r.yfit = val + r.yfit;
            end
        elseif length(val)==2
            r.x = r.x + val(1);
            r.y = r.y + val(2);
            if ~isempty(r.yfit)
                r.yfit = r.yfit + val(2);
            end
        else
            error('spec1d:plus:InvalidNumberArray','Length of constant is invalid. See documentation')
        end
    end
    s_out(i) = feval(class(r),r);
end
