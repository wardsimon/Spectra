function varargout = weigted_mean(varargin)
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
        m(i) = sum(s(i).y./(s(i).e.^2))./sum(1./(s(i).e.^2));
        e(i) = sqrt(1./sum(1./(s(i).e.^2)));
    end
    
    varargout{1} = m;
    if nargout == 2
        varargout{2} = e;
    end
end