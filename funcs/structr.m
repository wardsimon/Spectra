function out = structr(varargin)
    narginchk(2,3)
    
    if nargin == 3
        s=varargin{1};
        ind=varargin{2};
        r=varargin{3};
    elseif nargin==2
        s=varargin{1};
        r=varargin{2};
        ind =1:length(s);
    end
    
    ind = ind(:)';
    if ~ischar(r)
        error('Field must be a string')
    end
    j=1;
    for i=ind
        if i>length(s)
            error('Index %i is greater than the length of the structure.',i)
        end
        if iscell(s);
            t=s{i};
        else
            t=s(i);
        end
        if isfield(t,r)
            out(j,:)=t.(r);
        else
            out(j,:)=NaN;
            warning('Index %i could not be extracted.',i)
        end
        j=j+1;
    end
    if iscell(out)
        out=cell2mat(out);
    end
end