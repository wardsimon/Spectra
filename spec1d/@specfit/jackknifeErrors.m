function [s, obj2] = jackknifeErrors(obj,s,varargin)
% Perform a n-1 jackknife on the dataset

if ~exist('s','var')
    try
        s = evalin('caller','s');
    catch
        error();
    end
else
    if ~isa(s,'spec1d')
        temp = cell(1,length(varargin)+1);
        temp{1} = s;
        temp{2:end} = varargin{:};
        varargin = temp;
        try
            s = evalin('caller','s');
        catch
            error();
        end
    end
end

thetai = zeros(length(s.x),length(obj.pvals));
for i = 1:1:length(obj.pvals);
    if obj.notfixed(i)
        nf = zeros(size(obj.pvals));
        nf(i) = 1;
    else
        continue
    end
    if isempty(varargin)
        [~, b] = arrayfun(@(j)fits(s.set('x',s.x([1:(j-1),(j+1):length(s.x)]),'y',s.y([1:(j-1),(j+1):length(s.x)]),'e',s.e([1:(j-1),(j+1):length(s.x)])),obj.func,obj.pvals,nf),1:length(s.x),'UniformOutput',0);
    else
        [~, b] = arrayfun(@(j)fits(s.set('x',s.x([1:(j-1),(j+1):length(s.x)]),'y',s.y([1:(j-1),(j+1):length(s.x)]),'e',s.e([1:(j-1),(j+1):length(s.x)])),obj.func,obj.pvals,nf,varargin{:}),1:length(s.x),'UniformOutput',0);
    end
    thetai(:,i) = cellfun(@(x)x.pvals(i),b);
end

% Standard blocked error estimate.
se = sqrt(((length(s.x)-1)/length(s.x))*sum(bsxfun(@minus,thetai,mean(thetai)).^2));
% Bias in parameters.
bias = (length(s.x)-1)*(mean(thetai)-obj.pvals');

obj2 = obj.copy;
obj2.pvals = obj2.pvals - bias';
obj2.evals = se(:);

s = s.set('yfit',feval(obj2.func,s.x,obj2.pvals));
r = corrcoef([s.y(:),s.yfit(:)]);
obj2.rsq = r(1,2).^2;
v = length(s.y)-sum(logical(obj.notfixed));
obj2.chisq = sum(((s.y-s.yfit)./s.e).^2 )/v;
s = s.set('fitdata',obj2);
end

