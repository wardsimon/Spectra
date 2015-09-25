function varargout = min(varargin)
    k=1;
    for i=1:length(varargin)
        temp=varargin{i};
        for j=1:length(temp)
            if isa(temp(j),'spec1d')
                s(k)=temp(j);
                k=k+1;
            end
        end
    end
    
    for i=1:length(s)
        [y_max, ind]=min(s.y);
        if nargout==1
            varargout{1}(i)=y_max;
        elseif nargout<=2
            varargout{1}(i)=y_max;
            varargout{2}(i)=ind;
        elseif nargout<=3
            varargout{1}(i)=y_max;
            varargout{2}(i)=ind;
            varargout{3}(i)=s(i).x(ind);
        end
    end
end