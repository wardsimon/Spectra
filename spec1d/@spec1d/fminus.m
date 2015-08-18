function s_out  = fminus(varargin)
    
    k = 1;
    for i = 1:length(varargin)
        temp = varargin{i};
        if isa(temp,'spec1d')
            for j= 1:length(temp)
                s_in(k) = temp(j);
                k = k+1;
            end
        elseif isa(temp,'struct')
            for j = 1:length(temp)
                p(j) = temp(j);
            end
        end
    end
    
    for i = 1:length(s_in)
        x = getfield(s_in(i),'x'); x = x(:);
        y = getfield(s_in(i),'y'); y = y(:);
        e = getfield(s_in(i),'e'); e = e(:);
        if length(p) == 1
            y_calc = feval(p.function,x,p.pvals);
            s_in(i) = setfield(s_in(i),'yfit',y_calc);
            [e1, e2] = confidence(s_in(i),p,0.05);
        else
            y_calc = feval(p(i).function,x,p(i).pvals);
            s_in(i) = setfield(s_in(i),'yfit',y_calc);
            [e1, e2] = confidence(s_in(i),p(i),0.05);
        end
        
        e1 = e1(:); e2 = e2(:);
        
        % We have differenct upper and lower bounds. Have to try something!
        e = mean([sqrt(e.^2 + (y-e1).^2) sqrt(e.^2 + (y-e2).^2)],2);
        
        s_out(i) = setfield(s_in(i),'y',y-y_calc);
        s_out(i) = setfield(s_out(i),'e',e);
        s_out(i) = setfield(s_out(i),'yfit',[]);
    end
end