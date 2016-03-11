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
% running sdext.getCriteria. The lefault is least_squares.
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
p.addParamValue('dfdp','specdfdp_multi',@(x) ischar(x))
p.addParamValue('confidence',0.05,@(x) isnumeric(x) & length(x)==1 & x>0 & x<1)
p.addParamValue('inequc',{zeros(length(pin),0) zeros(0,1)},@iscell)
%     p.addParamValue('sep',cell(length(s1),1), @iscell)
p.addParamValue('optimiser','spec_lm',@ischar)
p.addParamValue('window',0,@(x) (isnumeric(x) && length(x)==2) || all(cellfun(@length,x)==2))
p.addParamValue('parallel',0,@(x) x==0 | x == 1)
p.addParamValue('criteria','least_square',@(x) ischar(x))
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

% Is this a multifit?
if options.multifit
    if options.verbose
        warning('spec1d:fits:Multifit','This fit is a multifit.')
    end
    [s1, pin, notfixed] = multifit_ini(s1,pin,notfixed);
    options.inequc = {zeros(length(pin),0) zeros(0,1)};
    if iscell(options.bounds)
        options.bounds = multifit_bounds(options.fixed,options.bounds);
    else
        options.bounds = [-Inf(length(pin),1) Inf(length(pin),1)];
    end
end

%----- Loop over the number of spec1d objects
sout = repmat(feval(class(s1(1))),size(s1));
fitdata = struct;

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
    parfor n = 1:length(all_data)
        switch  lower(options.optimiser)
            case 'spec_lm'
                [yfit{n},p{n},cvg,iter,corp,covp,covr,stdresid,Z,RSq{n},ra2{n},sig{n}] = speclsqr(all_data{n}{1},all_data{n}{2},all_data{n}{3},pin,notfixed,func,fcp,options); %#ok<*ASGLU>
            otherwise
                warning('spec1d:fits:InvalidOptimiser','This optimiser has not been implemented for parallel fitting.')
                [yfit{n},p{n},cvg,iter,corp,covp,covr,stdresid,Z,RSq{n},ra2{n},sig{n}] = speclsqr(all_data{n}{1},all_data{n}{2},all_data{n}{3},pin,notfixed,func,fcp,options);
        end
        ChiSq{n} = sum(((all_data{n}{2}-yfit{n})./all_data{n}{3}).^2 )/(length(all_data{n}{2})-sum(logical(notfixed)));
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
                    pnames = [pnames(:); temp(:)];
                end
            else
                [~,~,pnames] = feval(func,all_data{il}{1},p,1);
            end
        else
            pnames = cell(length(p),1);
            for i=1:length(p)
                pnames{i} = ['p' num2str(i)];
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
        
        fitdata(il).pvals = p{il};
        fitdata(il).evals = sig{il};
        fitdata(il).function = func;
        fitdata(il).pnames = pnames;
        fitdata(il).chisq = ChiSq{il};
        fitdata(il).rsq = RSq{il};
    end
    
else % This is for normal fitting and multi-fitting.
    for il=1:length(s1)
        x=s1(il).x; y=s1(il).y; e=s1(il).e;
        %----- Remove zeros from e
        x(e==0 | isinf(e) | isnan(e))=[];
        y(e==0 | isinf(e) | isnan(e))=[];
        e(e==0 | isinf(e) | isnan(e))=[];
        
        % Do we have a fitting window
        if options.window
            e = e(x>=options.window(1) & x<=options.window(2));
            y = y(x>=options.window(1) & x<=options.window(2));
            x = x(x>=options.window(1) & x<=options.window(2));
        end
        
        % When we just want to see the starting point
        if sum(notfixed) == 0
            yfit = feval(func,x,pin);
            p = pin;
            r = corrcoef([y(:),yfit(:)]);
            RSq = r(1,2).^2;
            sig = zeros(size(p));
        else % The real fit
            if strcmpi(options.optimiser,'builtin') && exist('optim','dir')
                warning('spec1d:fits:InvalidOptimiser','Optimiser %s is invalid as the optimisation toolbox is not installed. Defaulting to spec_lm.',options.optimiser)
                options.optimiser = 'spec_lm';
            end
            %----- Fit data
            switch  options.optimiser
                case 'spec_lm'
                    [yfit,p,cvg,iter,corp,covp,covr,stdresid,Z,RSq,ra2,sig] = speclsqr(x,y,e,pin,notfixed,func,fcp,options);
                case 'builtin'
                    infun = @(x,xin) feval(func,xin,x);
                    opt = optimset('MaxIter',fcp(2),'Display','Off');
                    [p,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit( @(p1, xdata) infun( interlace(pin(:), p1(:), ~notfixed), xdata),...
                        pin(notfixed), x, y, options.bounds(notfixed,1), options.bounds(notfixed,2),opt);
                    if exitflag ~= 1
                        warning('spec1d:fits:AlgorithmError','A minimum has not been found or an error has occoured. Error code %i. See documentation',exitflag)
                    end
                    p = interlace( pin(:), p(:), ~notfixed);
                    yfit = feval(func,x,p);
                    r = corrcoef([y(:),yfit(:)]);
                    RSq = r(1,2).^2;
                    wt = 1./e(:);
                    Qinv = diag(wt.*wt);
                    Q = diag((0*wt+1)./(wt.^2));
                    m = length(y);
                    n = sum(notfixed);
                    covr = residual'*Qinv*residual*Q/(m-n);                 %covariance of residuals
                    Vy = 1/(1-n/m)*covr;  % Eq. 7-13-22, Bard         %covariance of the data
                    jtgjinv = pinv(jacobian'*Qinv*jacobian);
                    covp = jtgjinv*jacobian'*Qinv*Vy*Qinv*jacobian*jtgjinv; % Eq. 7-5-13, Bard %cov of parm est
                    sig = interlace(zeros(size(pin)),sqrt(diag(covp)),~notfixed);
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
                        criteria_func = @(pars) sum(feval(options.criteria,y,e,FitFuncEval(func,x,pars,notfixed)));
                    else
                        criteria_func = @(pars) feval(options.criteria,y,e,FitFuncEval(func,x,pars,notfixed));
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
                    yfit = FitFuncEval(func,x,p,notfixed);
                    sig = output.parsHessianUncertainty;
                    r = corrcoef([y(:),yfit(:)]);
                    RSq = r(1,2).^2;
                    if nargout == 3
                        optional = output;
                    end
            end
        end
        % Goodness of fit
        v = length(y)-sum(logical(notfixed));
        ChiSq = sum(((y-yfit)./e).^2 )/v;
        
        % Get names of fit variables
        if pn
            if options.multifit
                pnames = {};
                for i = 1:length(s1)
                    [p_new, ~, ind] = multifitp2p(p,notfixed,i);
                    [~, ~, temp] = feval(func,x(ind),p_new,1);
                    pnames = [pnames(:); temp(:)];
                end
            else
                if nargin(func)<3
                    pnames = num2str((1:numel(p))','p%d');
                else
                    [~,~,pnames] = feval(func,x,p,1);
                end
            end
        else
            pnames = num2str((1:numel(p))','p%d');
        end
        
        % Make  the return spec1d objects
        sloop.x = x;
        sloop.y = y;
        sloop.e = e;
        sloop.x_label = s1(il).x_label;
        sloop.y_label = s1(il).y_label;
        sloop.datafile = s1(il).datafile;
        sloop.yfit = yfit;
        sout(il) = feval(class(sout(il)),sloop);
        
        fitdata(il).pvals = p;
        fitdata(il).evals = sig;
        fitdata(il).function = func;
        fitdata(il).pnames = pnames;
        fitdata(il).chisq = ChiSq;
        fitdata(il).rsq = RSq;
    end
end

% If a multifit, convert back into single objects.
if options.multifit
    [sout, fitdata] = multifit_extract(sout,fitdata,notfixed);
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