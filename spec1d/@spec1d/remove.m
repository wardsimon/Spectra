function varargout = remove(varargin)
    k=1;
    for i=1:length(varargin)
        temp=varargin{i};
        for j=1:length(temp)
            if isa(temp(j),'spec1d')
                s(k)=temp(j);
                k=k+1;
            elseif isa(temp(j),'numeric') | isa(temp(j),'logical')
                rm(j)=temp(j);
            end
        end
    end
    for i=1:length(s)
       varargout{1}(i)=spec1d(s(i).x(~rm),s(i).y(~rm),s(i).e(~rm)); 
    end
end
            
            
        
    