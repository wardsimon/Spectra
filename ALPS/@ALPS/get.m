function varargout=get(varargin)

r=1;    m=1;
f='';
for i=1:nargin
    temp=varargin{i};
    if isa(temp,'ALPS')
        for j=1:length(temp)
            in(r)=temp(j);
            r=r+1;
        end
    elseif isa(temp,'char')  
        f{m}=temp;
        m=m+1;
    end
end
r=r-1;

k=1;
for i=1:r
    fn=fieldnames(in(i));
    for j=1:length(fn)
        for m=1:length(f)
            if strcmp(f{m},fn(j))
                varargout{k}=eval(sprintf('in(%i).%s',i,fn{j}));
                k=k+1;
            end
        end
    end
end
        