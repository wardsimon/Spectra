function s = bootstrapErrors(obj,s,varargin)

if any(strcmp('Statistics and Machine Learning Toolbox',...
        arrayfun(@(x)x.Name,ver,'UniformOutput',0))) 
    s = bootstrapErrorsST(obj,s,varargin{:});
    return
end

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

if isempty(obj)
    error()
end


p = inputParser;
p.KeepUnmatched = true;
p.addParamValue('sigma',0.9544,@isnumeric)
p.addParamValue('plot',0,@isnumeric)
p.addParamValue('iter',200,@isnumeric)
p.addParamValue('biasCor',1,@isnumeric)
p.addParamValue('method','p',@ischar)
p.parse(varargin{:});
r = p.Results;
n = r.iter;
conf = r.sigma; % 2 sigma.
if conf >= 1
    conf = 1 - (1-erf(conf/sqrt(2)));
end
biasCorrect = r.biasCor;
plot_bootstrap = r.plot;
method = r.method;


if ~isempty(p.Unmatched)
    varargin = {};
else
    f = fieldnames(p.Unmatched);
    varargin = cell(1,2*length(f));
    for i = 1:2:(2*length(f))
        varargin{i} = f{i};
        varargin{i+1} = p.Unmatched.(f{i});
    end
end

res = s.y-feval(obj.func,s.x,obj.pvals);
sd = std(res);

T = NaN(n,length(obj.pvals));
TT = zeros(n,length(obj.pvals));
fprintf('Bootstrapping fit with %i itterations. This may take some time...\n',n)

% Loop through all available parameters
for i = 1:1:length(obj.pvals);
    if ~obj.notfixed(i)
        continue
    end
    nf = zeros(size(obj.pvals));
    nf(i) = 1;
    % Create a synthetic dataset.
    switch lower(method(1))
        case 'e' % empirical
            ind = ceil(length(s.x) * rand(length(s.x),n));
            rsynth = zeros(size(ind));
        case 'p' % parametric
            rsynth = sd*randn(length(s.x),n);
            ind = ceil(length(s.x) * rand(length(s.x),n));
        case 'b' % Baysean
            varargin{end +1} = 'criteria';
            varargin{end +1} = @bayesianMinimiser;
            ind = ceil(length(s.x) * rand(length(s.x),n));
            %             ind = repmat((1:length(s.x))',1,n);
            rsynth = zeros(size(ind));
        otherwise
            error
    end
    r = arrayfun(@(i)s.set('x',s.x(ind(:,i)),'y',s.y(ind(:,i))+rsynth(:,i),'e',s.e(ind(:,i))),1:n,'UniformOutput',0);
    r = [r{:}];
    % Fit the synthetic datasets
    if isempty(varargin)
        [~,b] = arrayfun(@(x)fits(x,obj.func,obj.pvals,nf),r,'UniformOutput',0);
    else
        [~,b] = arrayfun(@(x)fits(x,obj.func,obj.pvals,nf,varargin{:}),r,'UniformOutput',0);
    end
    % Get Synthetic fit values.
    T(:,i) = cellfun(@(x)x.pvals(i),b);
    TT(:,i) = sort(T(:,i));
end

if biasCorrect % Semi-Jackknife approach to remove cross correlation.
    fprintf('Calculating the bias-corrected and accelerated parameters z0 and a. Please wait...\n')
    z0 = sqrt(2)*erfinv(2*(sum(bsxfun(@(x,y) x<y,TT,obj.pvals'))/n)-1); % bias correction.
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
    a = sum((bsxfun(@minus,mean(thetai),thetai)).^3)./(6*(sum( (bsxfun(@minus,mean(thetai),thetai)).^2).^(3/2)));
else
    z0 = 0;
    a = 0;
end
% Calculate confidence ranges
zLo = sqrt(2)*erfinv(2*conf-1);
zHi = sqrt(2)*erfinv(2*(1-conf)-1);
zClo = z0 + (z0+zLo)./(1-a.*(z0+zLo));

[m,nn] = size(zClo);
z = (zClo - 0*ones(m,nn))/sqrt(1);
bcaLo = 0.5 + erf(z/sqrt(2))/2;

zChi = z0 + (z0+zHi)./(1-a.*(z0+zHi));
z = (zChi - 0*ones(m,nn))/sqrt(1);
bcaHi = 0.5 + erf(z/sqrt(2))/2;

% Set values and plot if needed.
obj2 = obj.copy;
for i = 1:length(obj2.pvals)
    if obj2.notfixed(i)
        obj2.pvals(i) = median(TT(:,i));
        indi = round([n*bcaHi(i) n*bcaLo(i)]);
        indi(indi==0) = 1;
        obj2.evals(i) = max(abs(TT(indi,i)-obj2.pvals(i)));
        if plot_bootstrap
            figure;
            histogram(TT(:,i))
            hold on
            plot([obj2.pvals(i) obj2.pvals(i)],get(gca,'YLim'),'r-')
            plot(TT([indi(1) indi(1)],i),get(gca,'YLim'),'b-')
            plot(TT([indi(2) indi(2)],i),get(gca,'YLim'),'b-')
            try
                title(sprintf('Bootstrap of parameter %s',obj2.pnames(i,:)))
                xlabel(obj2.pnames(i,:))
            catch
                try
                    title(sprintf('Bootstrap of parameter %s',obj2.pnames{i}))
                    xlabel(obj2.pnames{i})
                catch
                    xlabel(obj2.pnames)
                end
            end
            ylabel(sprintf('Occourance/%i',n))
        end
    end
end
r = corrcoef([s.y(:),s.yfit(:)]);
v = length(s.y)-sum(logical(obj.notfixed));
s = s.set('yfit',feval(obj2.func,s.x,obj2.pvals));
%,'rsq',r(1,2).^2,'chisq',sum(((s.y-s.yfit)./s.e).^2 )/v);
s = s.setfitdata(obj2);

end

function r = drchrnd(a,n)
% take n samples from a dirichlet distribution with length of a.
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);
end

function f = bayesianMinimiser(Signal, Error, Model)
f = least_square(Signal, Error./drchrnd(ones(size(Signal(:)))',1)', Model);
end