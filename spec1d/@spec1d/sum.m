function s=sum(varargin)
k=1;
for i=1:length(varargin)
    for j=1:length(varargin{i})
        if isa(varargin{i}(j),'spec1d')
            s1(k)=varargin{i}(j);
            k=k+1;
        end
    end
end

for i=1:(k-1)
    r.x(i)=i;
    r.y(i)=sum(s1(i).y);
    r.e(i)=sum(s1(i).e);
    if ~isempty(s1(i).yfit)
        r.yfit(i)=sum(s1(i).yfit);
    end
end
s=spec1d(r);