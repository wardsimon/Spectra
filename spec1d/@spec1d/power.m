function s_out = power(s,p)

    for i = 1:length(s)
        r = s(i);
        r.y = r.y.^p;
        if any(~isreal(r.y))
            error('spec1d:power:imaginaryY','Some or all y-values in the spec1d object are imaginary. spec1d does not support imaginary numbers.')
        end
        r.e = p*r.y .* (s(i).e./s(i).y);
        if any(isnan(r.e)) || any(isinf(r.e)) || any(~isreal(r.e))
            warning('spec1d:power:errorError','Some or all e-values in the spec1d object are imaginary/NaN/Inf. These points have been removed.')
            ind = isnan(r.e) | isinf(r.e) | ~isreal(r.e);
            r.x(ind) = [];
            r.y(ind) = [];
            r.e(ind) = [];
            if ~ isempty(r.yfit)
                r.yfit(ind)  = [];
            end
        end
        if ~ isempty(r.yfit)
            r.yfit  = r.yfit.^p;
        end
        s_out(i) = feval(class(r),r);
    end
    
        