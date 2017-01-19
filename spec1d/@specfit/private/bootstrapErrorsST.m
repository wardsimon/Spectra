function s = bootstrapErrorsST(obj,s,varargin)

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
p.addParamValue('sigma',1,@isnumeric)
p.addParamValue('plot',0,@isnumeric)
p.addParamValue('iter',200,@isnumeric)
p.addParamValue('biasCor',0,@isnumeric)
p.addParamValue('method','p',@ischar)
p.addParamValue('fast',0,@isnumeric)
p.addParamValue('parallel',0,@isnumeric)

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
fast = r.fast;

if ~isempty(varargin) && fast 
    fitfunc = @fastfit;
    getp = @(ff) ff;
else
    fitfunc = @fits;
    getp = @(ff) ff.fitdata.pvals';
end

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

switch lower(method(1))
    case 'e' % empirical
        if isempty(varargin)
            fitit = @(xp) feval(fitfunc,subindex(s,xp),obj.func,obj.pvals,obj.notfixed);
        else
            fitit = @(xp) feval(fitfunc,subindex(s,xp),obj.func,obj.pvals,obj.notfixed,varargin{:});
        end
        ind = 1:length(s.x);
    case 'p' % parametric - For normally distributed residuals
        foo = fitdist(res, 'normal');
        if isempty(varargin)
            fitit = @(xp) feval(fitfunc,setfield(s,'y',bsxfun(@minus,get(s,'y'),random(foo, [length(xp) 1]))),obj.func,obj.pvals,obj.notfixed);
        else
            fitit = @(xp) feval(fitfunc,setfield(s,'y',bsxfun(@minus,get(s,'y'),random(foo, [length(xp) 1]))),obj.func,obj.pvals,obj.notfixed,varargin{:});
        end
        ind = res;
    case 'b' % Baysean
        varargin2 = varargin;
        varargin2{end +1} = 'criteria';
        varargin2{end +1} = @bayesianMinimiser;
        fitit = @(xp) fits(subindex(s,ceil(length(s.x) * rand(length(xp),1))),obj.func,obj.pvals,obj.notfixed,varargin2{:});
        ind = zeros(size(s.x));
    otherwise
        error('The bootstrapping method is not recognised.')
end

pstat = bootstrp(n,@(in) getp(fitit(in)),ind,'Options',statset('UseParallel',r.parallel));

if biasCorrect % Semi-Jackknife approach to remove cross correlation.
    pstat = sort(pstat);
    if isempty(varargin)
        fitit = @(xp) feval(fitfunc,subindex(s,xp),obj.func,obj.pvals,obj.notfixed);
    else
        fitit = @(xp) feval(fitfunc,subindex(s,xp),obj.func,obj.pvals,obj.notfixed,varargin{:});
    end
    ind = 1:length(s.x);
    [ci,thetai] = bootci(n,{@(in) getp(fitit(in)),ind},'alpha',conf,'Options',statset('UseParallel',r.parallel));
    jak_e = max(abs(bsxfun(@minus,ci,mean(thetai))));
end

for i = 1:size(pstat,2)
    pd(i) = fitdist(pstat(:,i),'normal');
    if plot_bootstrap
        figure
        histogram(pstat(:,i),'Normalization','pdf','DisplayName',sprintf('Histogram p(%i)',i))
        hold on
        if biasCorrect
            histogram(thetai(:,i),'Normalization','pdf','DisplayName',sprintf('Histogram jk p(%i)',i))
        end
        ax = gca;
        x_values = linspace(ax.XLim(1),ax.XLim(2),200);
        y_values = pdf(pd(i),x_values);
        plot(x_values,y_values,'LineWidth',2,'DisplayName','Bootstrap')
        if biasCorrect
            plot(x_values,pdf('Normal',x_values,mean(thetai(:,i)),jak_e(i)),'LineWidth',2,'DisplayName','BCa max')
        end
        plot(x_values,pdf('Normal',x_values,obj.pvals(i),obj.evals(i)),'LineWidth',2,'DisplayName','Gaussian')
    end
end

obj2 = specfit(arrayfun(@(x) x.mu,pd)',arrayfun(@(x) x.sigma,pd)',obj.func, obj.notfixed,s);
s = set(s,'yfit',feval(obj.func,s.x,arrayfun(@(x) x.mu,pd)));
s = s.setfitdata(obj2);

end

function r = drchrnd(a,n)
% take n samples from a dirichlet distribution with length of a.
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);
end

function y = normpdf(x,mu,sigma)
 y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
end

function r = normrand()
% Avoid    a+(b-a)*rand   in case   a-b > realmax
a2 = a/2;
b2 = b/2;
mu = a2+b2;
sig = b2-a2;

r = mu + sig .* (2*rand(sizeOut,'like',mu)-1);

% Fill in elements corresponding to illegal parameter values
if ~isscalar(a) || ~isscalar(b)
    r(a > b) = NaN;
elseif a > b
    r(:) = NaN;
end
end

function f = bayesianMinimiser(Signal, Error, Model)
f = least_square(Signal, Error./drchrnd(ones(size(Signal(:)))',1)', Model);
end