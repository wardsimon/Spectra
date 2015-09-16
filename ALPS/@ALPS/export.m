function varargout = export(varargin)

r=1;    m=1;    n=1;
sstring='';
ffn=false;
for i=1:nargin
    temp=varargin{i};
    % Find filenames
    if isa(temp,'char') && ~ffn
        filename{n}=temp;
        n=n+1;
    end
    
    % Find data
    if isa(temp,'ALPS')
        ffn=true;
        for j=1:length(temp)
            in(r)=temp(j);
            r=r+1;
        end
    end
    
    % Find save options
    if isa(temp,'char') && ffn
        if strcmp(temp(1),'-')
            sstring=sprintf('%s , ''%s''',sstring,temp);
        else
            s(m)=temp;
            m=m+1;
        end
    end
end
r=r-1;  n=n-1;

if n~=r
    error('You need equal number of filenames to data files!');
end


if isempty(sstring)
    sstring=',''-ascii'', ''-double'', ''-tabs''';
end

if m>3
    warning('More than 2 outputs selected. Errors will occour if they are not the same size!')
end

for i=1:r
    for j=1:length(s)
        v(:,j)=get(in(i),s(j));
    end
    try
        eval(sprintf('save(''%s'',''v''%s)',filename{i},sstring));
        varargout{i}=0;
    catch
        varargout{i}=1;
    end
end
