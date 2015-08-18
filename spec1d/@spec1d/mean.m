function varargout = mean(varargin)
    k = 1;
    for i = 1:length(varargin)
        temp = varargin{i};
        for j = 1:length(temp)
            if isa(temp(j),'spec1d')
                s(k) = temp(j);
                k = k+1;
            end
        end
    end
    m = zeros(size(s));
    e = zeros(size(s));
    
    for i=1:length(s)
        m(i) = mean(s(i).y);
        if all(abs(diff(s(i).e))==abs((s(i).e(2)-s(i).e(1))))
            % The standard error
        e(i) = std(s(i).e)/sqrt(length(s(i).e));
        else
            % General error
            e(i) = sqrt(sum(s(i).e.^2));
        end
    end
    varargout{1} = m;
    if nargout == 2
        varargout{2} = e;
    end
end