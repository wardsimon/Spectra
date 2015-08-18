function v = isspec(varargin)
    
    k=1;
    for i=1:nargin
        temp=varargin{i};
        for j=1:length(temp)
            v(k)=isa(temp(j),'spec1d');
            k=k+1;
        end
    end