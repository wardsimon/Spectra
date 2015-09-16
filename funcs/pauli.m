function m = pauli(varargin)
    t.x=[0 1 ; 1 0];
    t.y=[0 -1j; 1j 0];
    t.z=[1 0; 0 -1];
    m=eye(2);
    for i=1:length(varargin)
        if iscell(varargin{i})
            for j=1:length(varargin{i})
                m=m*t.(varargin{i}{j});
            end
        else
            m=m*t.(varargin{i});
        end
    end
end
    