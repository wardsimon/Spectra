function varargout=divide(varargin)

k=1;
l=1;
for i=1:length(varargin)
    temp=varargin{i};
    if isa(temp,'spec1d')
        for j=1:length(temp)
            s(k)=temp(j);
            k=k+1;
        end
    elseif isa(temp,'double')
        dv(l:l+length(temp)-1)=temp;
        l=l+length(temp);
    end
end

if length(dv)==(k-1) || length(dv)==1
    if length(dv)==1;
        dv=dv*ones(1,(k-1));
    end
    for i=1:(k-1)
        varargout{1}(i)=(1/dv(i))*s(i);
    end
end