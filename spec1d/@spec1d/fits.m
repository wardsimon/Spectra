function [sout,fitdata,optional]=fits(s1,func,pin,notfixed,varargin)
%% function [sout,fitdata]=fits(s1,func,pin,notfixed,varargin)
%
% [sout,fitdata]=fits(s1,func,pin,notfixed,varargin)
%
% @SPEC1D/FIT Fits data in spec1d object s1 to a given fit function with
% starting parameters pin.
%
% Input parameters
%
% s        : Single or array of spec1d objects
% func     : Function of the form 'gauss', '@(x,p) gauss(x,p(1:4)) + p(5)*x',
%            @(x,p) gauss(x,p) + additionaldata.
% pin      : Parameters of func to be optimised. A double vector
% notfixed : Parameters to be varied. A logical vector.
% varargin : Additional options and parameters. See below...
%
% varargin are optional parameter, value pairs. Items marked with a * are 
% for the spec_lm optimiser. They can be the following:
% bounds      : Boundary conditions for pin. In the form of [Lower Upper]
%               column vectors size [length(pin) 2].
% dfdp*       : Function used to calculate the gradient of pin in func.
% confidence* : If errors can't be calculated through the jacobian a
%               confidence range is used.
% inequc*     : User supplied inequality relations.
% window      : Window to perform the fitting.
% parallel*   : Perform fits to the s1 vector in parallel
% verbose     : How much information is displayed.
%
% We can also vary the optimisation algorith and the criteria.
%
% The optional 'optimiser' command directs which optimiser to use. It can
% be any of the following:
% spec_lm    : General purpose implementation of the non-linear
%              Levenberg–Marquardt algorithm (Default)
% builtin    : Use the built in Levenberg–Marquardt algorithm. Will only
%              work if the Global Optimisation Toolbox is installed.
% sdext.getOptimisers : Running this command lists all available optimisation
%              algorithms. Default optimisation parameters cam be obtained 
%              by running options = algorithm('defaults'); where algorith is
%              an option in the list from sdext.getOptimisers. The options
%              can be passed into the fitting function by adding the
%              structure as the last varargin in the function call. 
%
% The optional 'criteria' can also be changed. Criteria can be found by
% running sdext.getCriteria. The default is least_squares.
%
% Note - Multifitting can also be performed. See the file
% multifit_example.m for worked examples.
%
% Version 4.3, Jan 2016
% Simon Ward, based on the work of Des McMorrow and Henrik Ronnow and
% Emmanuel FARHI

c = onCleanup(@cleanup);
global multifit_ind x_per_spec

%--- Parse the inputs and set defaults
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = true;
p.addRequired('spec1d',@(x) isa(x,'spec1d'))
if isa(func,'specfit')
    sf = true;
    f_in = func;
    func = f_in.func;
    pin = f_in.pin;
    notfixed = f_in.notfixed;
else
    sf = false;
end
p.addRequired('func',@(x) ischar(x) || isa(x,'function_handle'))
p.addRequired('pin',@(x) isnumeric(x) | iscell(x))
p.addRequired('fixed',@(x) isnumeric(x) | islogical(x))
fcp = [];
if length(varargin) ~= 1
    p.addParamValue('fcp',[0.0001 20 0.0001],@(x) length(x)==3) %#ok<*NVREPL>
else
    fcp = varargin{1};
end
p.addParamValue('bounds',[-Inf(length(pin),1) Inf(length(pin),1)],@(x) (size(x,1)==length(pin) & size(x,2)==2) | iscell(x))
p.addParamValue('dfdp','specdfdp_multi2',@(x) ischar(x) | isa(x,'function_handle'))
p.addParamValue('confidence',0.95,@(x) isnumeric(x) & length(x)==1 & x>0 & x<1)
p.addParamValue('inequc',{zeros(length(pin),0) zeros(0,1)},@iscell)
%     p.addParamValue('sep',cell(length(s1),1), @iscell)
p.addParamValue('optimiser','speclsqr',@(x) ischar(x) | isa(x,'function_handle'))
p.addParamValue('window',0,@(x) (isnumeric(x) && length(x)==2) || all(cellfun(@length,x)==2))
p.addParamValue('parallel',0,@(x) x==0 | x == 1)
p.addParamValue('criteria','least_square',@(x) ischar(x) | isa(x,'function_handle'))
p.addParamValue('verbose',0,@(x) isnumeric(x) | islogical(x))
if length(varargin) ~= 1
    p.parse(s1,func,pin,notfixed,varargin{:});
else
    p.parse(s1,func,pin,notfixed);
end
options = p.Results;
opt = p.Unmatched;

% Initial parameters
s1 = options.spec1d;
func = options.func;
pin = options.pin;
notfixed = options.fixed;
if any(notfixed > 1)
    options.multifit = 1;
    if ~iscell(pin)
        error('spec1d:fits:MalformedPin','The pin must be a cell vector.')
    end
else
    options.multifit = 0;
    if iscell(pin)
        error('spec1d:fits:MalformedPin','The pin must not be a cell vector.')
    end
end

% Do we need to change the functions?
pn = 1;
if ischar(func)
    if strcmp(func(1),'@')
        func = str2func(func);
        pn   = 0;
    end
end

if isempty(fcp)
    fcp = options.fcp;
end


f_in = specfit(func,pin,notfixed);
f_in.specID = [s1.ident];
f_in = f_in.setFitdata('Display',options.verbose,'MaxIter',options.fcp(2),...
                'TolGradCon',options.fcp(3),'Algorithm',options.optimiser,...
                'UseParallel',options.parallel,'TolX',options.fcp(1),...
                'CriteriaFcn',options.criteria,'JacobFcn',options.dfdp,...
                'Bounds',options.bounds,'TolCon',options.confidence);
f_in.userdata = struct('x_per_spec',[],'param_keep',[]);

% Is this a multifit?
if options.multifit
    if options.verbose
        warning('spec1d:fits:Multifit','This fit is a multifit.')
    end
%     [s1, pin, notfixed] = multifit_ini(s1,pin,notfixed);
%     f_in = f_in.copy;
%     f_in.specID = s1.ident;
%     f_in.pvals = pin;
%     f_in.evals = nan(size(pin));
%     f_in.notfixed = notfixed;
    [s1, f_in] = multifit_ini(s1,f_in);
    pin = f_in.pvals;
    notfixed = f_in.notfixed;
    options.inequc = {zeros(length(pin),0) zeros(0,1)};
    if iscell(options.bounds)
        options.bounds = multifit_bounds(options.fixed,options.bounds);
    else
        options.bounds = [-Inf(length(pin),1) Inf(length(pin),1)];
    end
    f_in = f_in.setFitdata('Bounds',options.bounds);
end

%----- Loop over the number of spec1d objects
sout = repmat(feval(class(s1(1))),size(s1));
% fitdata = struct;

% This is for parallel fitting.
if all([sdext.getpref('experimental').val, options.parallel, ~options.multifit])
    all_data = cell(length(s1),1);
    [all_data{:}] = get(s1,'x','y','e');
    for i = 1:length(all_data)
        all_data{i}{1}(isnan(all_data{i}{3}) | isinf(all_data{i}{3}) | all_data{i}{3}==0) = [];
        all_data{i}{2}(isnan(all_data{i}{3}) | isinf(all_data{i}{3}) | all_data{i}{3}==0) = [];
        all_data{i}{3}(isnan(all_data{i}{3}) | isinf(all_data{i}{3}) | all_data{i}{3}==0) = [];
        % Do we have a fitting window
        if ~isempty(options.window) && options.window ~= 0
            all_data{i}{3} = all_data{i}{3}(all_data{i}{1}>=options.window{i}(1) & all_data{i}{1}<=options.window{i}(2));
            all_data{i}{2} = all_data{i}{2}(all_data{i}{1}>=options.window{i}(1) & all_data{i}{1}<=options.window{i}(2));
            all_data{i}{1} = all_data{i}{1}(all_data{i}{1}>=options.window{i}(1) & all_data{i}{1}<=options.window{i}(2));
        end
    end
    % Do the fitting
    yfit = cell(length(all_data),1);
    RSq = cell(length(all_data),1);
    ra2 = cell(length(all_data),1);
    sig = cell(length(all_data),1);
    switch  lower(options.optimiser)
        case 'builtin'
            infun = @(x,xin) feval(func,xin,x);
            opt = optimset('MaxIter',fcp(2),'Display','Off','UseParallel',true);
            for n = 1:length(all_data)
                [p{n},resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit( @(p1, xdata) infun( interlace(pin(:), p1(:), ~notfixed), xdata),...
                    pin(notfixed),all_data{n}{1},all_data{n}{2}, options.bounds(notfixed,1), options.bounds(notfixed,2),opt);
                if exitflag ~= 1
                    warning('spec1d:fits:AlgorithmError','A minimum has not been found or an error has occoured. Error code %i. See documentation',exitflag)
                end
                p{n} = interlace( pin(:), p{n}(:), ~notfixed);
                yfit{n} = feval(func,all_data{n}{1},p{n});
                r = corrcoef([all_data{n}{2}(:),yfit{n}(:)]);
                RSq{n} = r(1,2).^2;
                wt = 1./all_data{n}{3}(:);
                Qinv = diag(wt.*wt);
                Q = diag((0*wt+1)./(wt.^2));
                m = length(all_data{n}{2});
                nn = sum(notfixed);
                covr = residual'*Qinv*residual*Q/(m-nn);                 %covariance of residuals
                Vy = 1/(1-nn/m)*covr;  % Eq. 7-13-22, Bard         %covariance of the data
                jtgjinv = pinv(jacobian'*Qinv*jacobian);
                covp = jtgjinv*jacobian'*Qinv*Vy*Qinv*jacobian*jtgjinv; % Eq. 7-5-13, Bard %cov of parm est
                sig{n} = interlace(zeros(size(pin)),sqrt(diag(covp)),~notfixed);
                ChiSq{n} = sum(((all_data{n}{2}-yfit{n})./all_data{n}{3}).^2 )/(length(all_data{n}{2})-sum(logical(notfixed)));
            end
        case 'speclsqr'
            parfor n = 1:length(all_data) %#ok<*PFUIX>
                [yfit{n},p{n},cvg,iter,corp,covp,covr,stdresid,Z,RSq{n},ra2{n},sig{n}] = speclsqr(all_data{n}{1},all_data{n}{2},all_data{n}{3},pin,notfixed,func,fcp,options); %#ok<*ASGLU>
                ChiSq{n} = sum(((all_data{n}{2}-yfit{n})./all_data{n}{3}).^2 )/(length(all_data{n}{2})-sum(logical(notfixed)));
            end
        otherwise
            error('spec1d:fits:InvalidOptimiser','This optimiser has not been implemented for parallel fitting.')
    end
    
    % Set the return
    for il = 1:length(all_data)
        % Get names of fit variables
        if pn
            if options.multifit
                pnames = {};
                for i = 1:length(s1)
                    [p_new, ~, ind] = multifitp2p(p{i},notfixed,i);
                    [~, ~, temp] = feval(func,all_data{il}{1}(ind),p_new,1);
                    pnames = {pnames{:}; temp{:}};
                end
            else
                [~,~,pnames] = feval(func,all_data{il}{1},p,1);
            end
        else
            for i = 1:numel(p)
                pnames{i} = num2str(i,'p%d');
            end
        end
        
        % make spec1d objects
        sloop.x = all_data{il}{1};
        sloop.y = all_data{il}{2};
        sloop.e = all_data{il}{3};
        sloop.x_label = s1(il).x_label;
        sloop.y_label = s1(il).y_label;
        sloop.datafile = s1(il).datafile;
        sloop.yfit = yfit{il};
        sout(il) = feval(class(sout(il)),sloop);

        fitdata(il) = specfit(p{il},sig{il},func,pnames,ChiSq{il},RSq{il});
        fitdata(il).notfixed = notfixed;
    end
    
else % This is for normal fitting and multi-fitting.
    for il=1:length(s1)
        
        s1(il) = clean(s1(il));
        if isempty(s1(il).x)
            error('spec1d:fits:NoDataPoints','There are no data points in the %s object, or all data points have been removed for being unsuitable for fitting',class(s1(il)))
        end
        
        % Do we have a fitting window
        if options.window
            s1(il) = cut(s1(il),options.window);
        end
        
        % When we just want to see the starting point
        if sum(notfixed) == 0
            f_out(i) = f_in.copy();
            yfit = feval(func,s1(il).x,pin);
            p = pin;
            r = corrcoef([s1(il).y(:),s1(il).yfit(:)]);
            RSq = r(1,2).^2;
            sig = zeros(size(p));
        else % The real fit
            if strcmpi(options.optimiser,'builtin') && exist('optim','dir')
                warning('spec1d:fits:InvalidOptimiser','Optimiser %s is invalid as the optimisation toolbox is not installed. Defaulting to spec_lm.',options.optimiser)
                options.optimiser = 'spec_lm';
            end
            %----- Fit data
            switch  options.optimiser
                case 'speclsqr'
%                     fitdata = speclsqr(f_in,s1(il));
                    [yfit,p,cvg,iter,corp,covp,covr,stdresid,Z,RSq,ra2,sig,f_out(il)] = speclsqr(s1(il),f_in,options);
                case 'builtin'
                    infun = @(x,xin) feval(func,xin,x);
                    opt = optimset('MaxIter',fcp(2),'Display','Off');
                    [p,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit( @(p1, xdata) infun( interlace(pin(:), p1(:), ~notfixed), xdata),...
                        pin(notfixed), s1(il).x, s1(il).y, options.bounds(notfixed,1), options.bounds(notfixed,2),opt);
                    if exitflag ~= 1
                        warning('spec1d:fits:AlgorithmError','A minimum has not been found or an error has occoured. Error code %i. See documentation',exitflag)
                    end
                    p = interlace( pin(:), p(:), ~notfixed);
                    yfit = feval(func,s1(il).x,p);
                    r = corrcoef([s1(il).y(:),yfit(:)]);
                    RSq = r(1,2).^2;
                    wt = 1./s1(il).e(:);
                    Qinv = diag(wt.*wt);
                    Q = diag((0*wt+1)./(wt.^2));
                    m = length(s1(il).y);
                    n = sum(notfixed);
                    covr = residual'*Qinv*residual*Q/(m-n);                 %covariance of residuals
                    Vy = 1/(1-n/m)*covr;  % Eq. 7-13-22, Bard         %covariance of the data
                    jtgjinv = pinv(jacobian'*Qinv*jacobian);
                    covp = jtgjinv*jacobian'*Qinv*Vy*Qinv*jacobian*jtgjinv; % Eq. 7-5-13, Bard %cov of parm est
                    sig = interlace(zeros(size(pin)),sqrt(diag(covp)),~notfixed);
                    f_out(il) = f_in.copy();
                otherwise
                    % Check that the criteria is possible.
                    if ~strcmpi(options.criteria,sdext.getCriteria)
                        error('spec1d:fits:InvalidCriteria','Criteria %s is invalid. See available criteria with sdext.getCriteria',options.criteria)
                    end
                    % Check that the optimiser is possible.
                    if ~strcmpi(options.optimiser,sdext.getOptimisers)
                        error('spec1d:fits:InvalidOptimiser','Optimiser %s is invalid. See available optimisers with sdext.getOptimisers',options.optimiser)
                    end
                    if strcmpi(options.criteria,'least_square')
                        criteria_func = @(pars) sum(feval(options.criteria,s1(il).y,s1(il).e,FitFuncEval(func,s1(il).x,pars,notfixed)));
                    else
                        criteria_func = @(pars) feval(options.criteria,s1(il).y,s1(il).e,FitFuncEval(func,s1(il).x,pars,notfixed));
                    end
                    % We may have specific options to pass to the
                    % optimiser.
                    if ~isfield(opt,'About')
                        opt = feval(options.optimiser,'defaults');
                    end
                    % Check for loose boundary conditions
                    if any(isinf(options.bounds(:)))
                        warning('spec1d:fits:BoundaryWarning','Some boundary conditions are infinate. These optimisers need well defined bounds. Taking +-2*pin')
                    else
                        constraints.min = options.bounds(:,1);
                        constraints.max = options.bounds(:,2);
                    end
                    constraints.fixed = ~notfixed;
                    % Do the fitting
                    [p, fval, flag, output ]= feval(options.optimiser,criteria_func,pin,setstructfields(opt,options),constraints);
                    if flag ~= 0
                        warning('spec1d:fits:AlgorithmError','A minimum has not been found or an error has occoured. Error code %i. See documentation',flag)
                    end
                    yfit = FitFuncEval(func,s1(il).x,p,notfixed);
                    sig = output.parsHessianUncertainty;
                    r = corrcoef([s1(il).y(:),yfit(:)]);
                    RSq = r(1,2).^2;
                    f_out(il) = f_in.copy();
                    f_out(il).fitdata = output;
%                     if nargout == 3
%                         optional = output;
%                     end
            end
        end
        
        % Get names of fit variables
        if pn
            if options.multifit
                pnames = {};
                for i = 1:length(s1)
                    [p_new, ~, ind] = multifitp2p(p,notfixed,i);
                    if nargin(func)<3
                        for j = 1:numel(p_new)
                            temp{j} = num2str(j,'p%d');
                        end
                    else
                        [~, ~, temp] = feval(func,s1(il).x(ind),p_new,1);
                    end
                    pnames = {pnames{:}; temp{:}};
                end
            else
                if nargin(func)<3
                    for i = 1:numel(p)
                        pnames{i} = num2str(i,'p%d');
                    end
                else
                    [~,~,pnames] = feval(func,s1(il).x,p,1);
                end
            end
        else
            for i = 1:numel(p)
                pnames{i} = num2str(i,'p%d');
            end
        end
        
        % Make  the return spec1d objects
        sout(il) = s1(il).copy;
        sout(il) = set(sout(il),'x',s1(il).x,'y',s1(il).y,'e',s1(il).e,'yfit',yfit);
        fitdata(il) = f_out(il);
        fitdata(il).specID = sout(il).ident;
        fitdata(il).pnames = pnames;  
%         fitdata(il) = specfit(p,sig,func,notfixed,sout(il).ident,pnames);
%         results = struct('pvals',p,'evals',sig,'func',func,'pnames',pnames,'chisq',ChiSq,'rsq',RSq,'notfixed',notfixed);
%         fitdata(il) = specfit(results);
        sout(il) = sout(il).setfitdata(fitdata(il));
    end
end

% If a multifit, convert back into single objects.
if options.multifit
    [sout, fitdata] = multifit_extract(sout,fitdata);
    x_per_spec = [];
    multifit_ind = [];
end
end

% Reset after finishing.
function cleanup
    global multifit_ind x_per_spec
    x_per_spec = [];
    multifit_ind = [];
end


function f = FitFuncEval(F,x,p,notfixed)
    global x_per_spec multifit_ind
    x_per_spec_local = x_per_spec;
    
    if isempty(x_per_spec_local)
        f = feval(F,x,p);
        f = f(:);
    else
        f = zeros(sum(x_per_spec_local),1);
        for i = 1:length(x_per_spec_local)
            [p_new, ~, ind] = multifitp2p(p,notfixed,i);
            multifit_ind = ind;
            f(ind) = feval(F,x(ind),p_new);
        end
    end
end


function a = interlace( a, x, fix )
    a(~fix) = x;
end