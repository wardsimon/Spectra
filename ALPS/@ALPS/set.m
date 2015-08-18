function varargout=set(varargin)

r=1;    m=1;    f='';
skip=0;
for i=1:nargin
    temp=varargin{i};
    if ~skip
        if isa(temp,'ALPS')
            for j=1:length(temp)
                in(r)=temp(j);
                r=r+1;
            end
        elseif isa(temp,'char')
            f{m}=temp;
            m=m+1;
            skip=1;
        end
    else
       in_c{m-1}=temp;
       skip=0;
    end
end
r=r-1;

if length(in_c)~=(m-1)
    error('Missing number of results to arguemnts')
end

k=1;
for i=1:r
    fn=fieldnames(in(i));
    out(k)=in(k);
    for j=1:length(fn)
        for m=1:length(f)
            if strcmp(f{m},fn(j))
                st=sprintf('out(%i).%s=in_c{%i};',k,fn{j},k);
                eval(st);
                varargout{k}=out(k);
                k=k+1;
            end
        end
    end
end
        