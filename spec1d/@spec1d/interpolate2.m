function s_out = interpolate(s_in, xnew, varargin)
    % function s=interpolate(s, xnew)
    %
    % SPED1D/Interpolates a spectrum to new x-values
    %
    
    if isa(xnew,'spec1d')
        xnew = xnew.x;
    end
    method = [];
    if nargin == 3
        if ischar(varargin{1})
            method = varargin{1};
        end
    end
    
    for i = 1:length(s_in)
        x_temp = s_in(i).x;
        y_temp = s_in(i).y;
        e_temp = s_in(i).e;        
        try
            if isempty(method)
                s(i).y = interp1(x_temp,y_temp,xnew);
            else
                s(i).y = interp1(x_temp,y_temp,xnew,method);
            end
            s(i).x = xnew(:);
            s(i).y = s(i).y(:);
            
            for j=1:length(xnew)
                n1=find(x_temp<=xnew(j));
                n2=find(x_temp> xnew(j));
                if isempty(n1)
                    n1=n2(1);
                    n2=n2(2);
                elseif isempty(n2)
                    n2=n1(end);
                    n1=n1(end-1);
                else
                    n1=n1(end);
                    n2=n2(1);
                end
                s(i).e(j)=sqrt((e_temp(n1)*(x_temp(n2)-xnew(j))/(x_temp(n2)-x_temp(n1)))^2+...
                    (e_temp(n2)*(x_temp(n1)-xnew(j))/(x_temp(n1)-x_temp(n2)))^2);
            end
            s_out(i)=spec1d(s(i));
        catch
            warning('spec1d:interpolate','Can not interpolate s(%i)',i)
        end
    end