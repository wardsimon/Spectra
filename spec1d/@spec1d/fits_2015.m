function [sout,fitdata]=fits_2015(s1,func,pin,notfixed,varargin)
    %% function [sout,fitdata]=fits(s1,func,pin,notfixed,options)
    %
    % [sout,fitdata]=fits(s1,func,pin,notfixed,options)
    %
    % @SPEC1D/FIT Fits data in spec1d object s1 to MFIT function
    % specified in func. Also works for an array of spec1d objects.
    %
    % Version 5, September 2015
    % Simon Ward simon.ward@psi.ch
    % Loosly based on the work of Des McMorrow and Henrik Ronnow, with
    % substantial inspiration from iFit http://ifit.mccode.org/
    
    % Cleanup function
    c = onCleanup(@cleanup);
    
    %     global multifit_nd x_per_spec
    
    % Set the options structure
    opt = optimset;
    opt = setstructfields(setstructfields(opt,setFitDefaults()),setCriteriaDefaults());
    f = fieldnames(opt);
    %--- Parse the inputs and set defaults
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true;
    
    % set required inputs
    p.addRequired('spec1d',@(x) isa(x,'spec1d'))
    p.addRequired('func',@ischar)
    p.addRequired('parsStart',@(x) isnumeric(x) | iscell(x))
    p.addRequired('Fixed',@(x) isnumeric(x) | islogical(x))
    
    % Add the fitting arguemts
    for i = 1:length(f)
        p.addParameter(f{i},opt.(f{i}));
    end
    
    % Evaluate arguments
    p.parse(s1,func,pin,notfixed,varargin{:});
    opt = p.Resuts;
    
    % The fitting objects
    [s, options, criteria] = splitOptions(opt);
    
    
    % Fix anonomous functions
    pn = 1;
    if ischar(func)
        if strncmp(func,'@',1)
            func = str2func(func);
            pn   = 0;
        end
    end
    
    
    %----- Loop over the number of spec1d objects
    sout    = repmat(spec1d,size(s));
    fitdata = repmat(struct,size(s));
    
    for il=1:length(s)
        x=s(il).x; y=s(il).y; e=s(il).e;
        %----- Remove zeros from e
        x(e==0 | isinf(e) | isnan(e))=[];
        y(e==0 | isinf(e) | isnan(e))=[];
        e(e==0 | isinf(e) | isnan(e))=[];
        
        % Do we have a fitting window
        if criteria.Window
            e = e(x>=criteria.window(il,1) & x<=criteria.window(il,2));
            y = y(x>=criteria.window(il,1) & x<=criteria.window(il,2));
            x = x(x>=criteria.window(il,1) & x<=criteria.window(il,2));
        end
        
        if sum(notfixed)==0
            yfit = feval(func,x,criteria.parsStart);
            p = criteria.parsStart;
            r = corrcoef([y(:),yfit(:)]);
            RSq = r(1,2).^2;
            sig = zeros(size(p));
        else
            
            % Algorithm
            % Criteria
            % User Function
            
            % we need to call the optimization method with the eval_criteria as FUN
            % call minimizer ===============================================================
            if abs(nargin(options.Algorithm)) == 1 || abs(nargin(options.Algorithm)) >= 6
                [pars_out,criteria,message,output] = feval(options.Algorithm, ...
                    @(pars) eval_criteria(func, pars, options.criteria, a, varargin{:}), pars, options, constraints);
            else
                % Constraints not supported by optimizer
                [pars_out,criteria,message,output] = feval(options.optimizer, ...
                    @(pars) eval_criteria(func, pars, options.criteria, a, varargin{:}), pars, options);
            end
        end
        
        yfit = feval(func, p, x);
        output.corrcoef   = eval_corrcoef(y, e, yfit);
        output.residuals  = y - yfit;
        Rwp  = sqrt(sum((output.residuals./e).^2)/sum((y./e).^2));
        output.Rfactor    = Rwp;
        Rexp = sqrt((length(y) - length(pars))./sum((y./e).^2));
        GoF  = Rwp/Rexp;
        if strcmp(options.Verbose)
            fprintf(1, ' Correlation coefficient=%g (closer to 1 is better)\n',  output.corrcoef);
            fprintf(1, ' Weighted R-factor      =%g (Rwp, smaller that 0.2 is better)\n', output.Rfactor);
            fprintf(1, ' Experimental R-factor  =%g (Rexp)\n', Rexp);
        end
        
        %----- Goodness of fit
        v = length(y)-sum(logical(notfixed));
        ChiSq = sum(((y-yfit)./e).^2 )/v;
        
        %----- Get names of fit variables
        if pn
            %                 if options.multifit
            %                     pnames = {};
            %                     for i = 1:n
            %                         [p_new, ~, ind] = multifitp2p(p,notfixed,i);
            %                         [~, ~, temp] = feval(func,x(ind),p_new,1);
            %                         pnames = [pnames(:); temp(:)];
            %                     end
            %                 else
            [~,~,pnames] = feval(func,x,p,1);
            %                 end
        else
            pnames = cell(length(p),1);
            for i=1:length(p)
                pnames{i} = ['p' num2str(i)];
            end
        end
        
        %----- set return
        sloop.x = x;
        sloop.y = y;
        sloop.e = e;
        sloop.x_label = s1(il).x_label;
        sloop.y_label = s1(il).y_label;
        sloop.datafile = s1(il).datafile;
        sloop.yfit = yfit;
        sout(il) = spec1d(sloop);
        
        fitdata(il).pvals = p;
        fitdata(il).evals = sig;
        fitdata(il).function = func;
        fitdata(il).pnames = pnames;
        fitdata(il).chisq = ChiSq;
        fitdata(il).rsq = RSq;
    end
    %     end
    %     if options.multifit
    %         [sout, fitdata] = multifit_extract(sout,fitdata,notfixed);
    %         x_per_spec = [];
    %         multifit_ind = [];
    %     end
    % end
    
    
    % format output arguments ======================================================
    % pars_out = reshape(pars_out, [ 1 numel(pars_out) ]); % row vector
    % model.ParameterValues = pars_out;
    % if ~isempty(inputname(1))
    %     try
    %   assignin('caller',inputname(1),model); % update in original object
    %     end
    % end
    %
    % if nargout > 3 || (isfield(options,'Diagnostics') && (strcmp(options.Diagnostics, 'on') || any(options.Diagnostics == 1)))
    %   output.modelValue = feval(model, pars_out, a.Axes{:});
    %   index=find(~isnan(a.Signal) & ~isnan(output.modelValue));
    %   if ~isscalar(a.Error), e = a.Error(index); else e=a.Error; end
    %   output.corrcoef   = eval_corrcoef(a.Signal(index), e, output.modelValue(index));
    %   output.residuals  = a.Signal - output.modelValue;
    %   % Rwp
    %   Rwp  = sqrt(sum((output.residuals(index)./e).^2)/sum((a.Signal(index)./e).^2));
    %   output.Rfactor    = Rwp;
    %   % Rexp
    %   Rexp = sqrt((length(a.Signal(index)) - length(pars))./sum((a.Signal(index)./e).^2));
    %   % Goodness of fit
    %   GoF  = Rwp/Rexp;
    %   if strcmp(options.Display, 'iter') | strcmp(options.Display, 'final') | ...
    %     (isfield(options,'Diagnostics') && (strcmp(options.Diagnostics, 'on') || any(options.Diagnostics == 1)))
    %     fprintf(1, ' Correlation coefficient=%g (closer to 1 is better)\n',  output.corrcoef);
    %     fprintf(1, ' Weighted R-factor      =%g (Rwp, smaller that 0.2 is better)\n', output.Rfactor);
    %     fprintf(1, ' Experimental R-factor  =%g (Rexp)\n', Rexp);
    %   end
    %   if abs(output.corrcoef) < 0.6 && ~isscalar(a.Error)
    %     name = inputname(2);
    %     if isempty(name)
    %       name = 'a';
    %     end
    %     fprintf(1, ' WARNING: The fit result is BAD. You may improve it by setting %s.Error=1\n',...
    %       name);
    %   end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EMBEDDED FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this way 'options' is available in here...
    
    function c = eval_criteria(model, p, criteria, x, varargin)
        
        % criteria to minimize
        if nargin<5, varargin={}; end
        % then get model value
        Model  = feval(model, p, x, varargin{:}); % return model values
        
        % get actual parameters used during eval (in case of Constraints)
        p = model.ParameterValues;
        % send it back to input call
        Model  = iFunc_private_cleannaninf(Model);
        if isempty(Model)
            error([ 'iFunc:' mfilename ],[ 'The model ' model ' could not be evaluated (returned empty).' ]);
        end
        
        % compute criteria
        c = feval(criteria, a.Signal(:), a.Error(:), Model(:));
        % divide by the number of degrees of freedom
        % <http://en.wikipedia.org/wiki/Goodness_of_fit>
        if numel(a.Signal) > length(p)-1
            c = c/(numel(a.Signal) - length(p) - 1); % reduced 'Chi^2'
        end
        
        % overlay data and Model when in 'OutputFcn' mode
        if (isfield(options, 'OutputFcn') && ~isempty(options.OutputFcn) && ~isscalar(a.Signal) && ndims(a.Signal) <= 2)
            if ~isfield(options, 'updated')
                options.updated   = -clock;
                options.funcCount = 0;
            end
            
            options.funcCount = options.funcCount+1;
        end
        
    end % eval_criteria (embedded)
    
end % fits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = eval_corrcoef(Signal, Error, Model)
    % correlation coefficient between the data and the model
    
    % compute the correlation coefficient
    if isempty(Error) || isscalar(Error) || all(Error(:) == Error(end))
        wt = 1;
    else
        wt = 1./Error;
        wt(~isfinite(wt)) = 0;
    end
    r = corrcoef(Signal.*wt,Model.*wt);
    r = r(1,2);                                     % correlation coefficient
    if isnan(r)
        r = corrcoef(Signal,Model);
        r = r(1,2);
    end
end % eval_corrcoef

function cleanup
    global multifit_ind x_per_spec
    x_per_spec = [];
    multifit_ind = [];
end

function sout = setFitDefaults()
    sout.Display   = 'off';
    sout.Algorithm = 'fminpso';
    sout.UseMatlab = 0;
    sout.MultiFit   = 0;
end

function sout = setCriteriaDefaults()
    sout.Fixed      = 0;
    sout.parsStart  = 0;
    sout.Minimiser  = 'least_squares';
    sout.FuncCount  = 0;
    sout.Min        = [];
    sout.Max        = [];
    sout.Bounds     = [];
    sout.Window     = [-Inf Inf];
end

function [sout, options, criteria] = splitOptions(sin)
    % Spec1d
    sout = sin.spec1d;
    
    % Options
    options = setstructfields(setstructfields(optimset,setFitDefaults()),feval(sout.Algorithm,'defaults'));
    fnames = fieldnames(options);
    for i = 1:length(fnames)
        options.(fnames{i}) = sin.(fnames{i});
    end
    
    % Criteria
    criteria = setCriteriaDefaults();
    fnames = fieldnames(criteria);
    for i = 1:length(fnames)
        criteria.(fnames{i}) = sin.(fnames{i});
    end
    
    % MultiFit check and Bounds Fix
    if any(criteria.Fixed > 1)
        options.MultiFit = 1;
        if ~iscell(criteria.parsStart)
            error('The pin must be a cell')
        end
        warning('Just FYI, this is a multifit!')
        [sout, criteria.parsStart,criteria.Fixed] = multifit_ini(sout,criteria.parsStart,criteria.Fixed);
        if iscell(criteria.Bounds)
            criteria.Bounds = multifit_bounds(criteria.Fixed,criteria.Bounds);
        else
            criteria.Bounds = bsxfun(@times,Inf(length(criteria.parsStart),2),[-1 1]);
        end
    else
        if iscell(criteria.parsStart)
            error('The pin must not be a cell')
        end
        if isempty(criteria.Bounds)
            criteria.Bounds = bsxfun(@times,Inf(length(criteria.parsStart),2),[-1 1]);
        else
            if length(criteria.Bounds) ~= length(criteria.parsStart)
                warning('Boundaries changed to -Inf Inf for all parameters')
                criteria.Bounds = bsxfun(@times,Inf(length(criteria.parsStart),2),[-1 1]);
            end
        end
    end
    criteria.Min = criteria.Bounds(:,1);
    criteria.Max = criteria.Bounds(:,2);
    
    % Check for parallel and MultiFit incompatibility
    if options.MultiFit && options.UseParallel
        warning('Canceling parralel fit due to MultiFitting')
        options.UseParallel = 0;
    end
end