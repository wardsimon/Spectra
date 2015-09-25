function [sout,fitdata]=fits_2015(s1,func,pin,notfixed,varargin)
    %% function [sout,fitdata]=fits(s1,func,pin,notfixed,fcp,win)
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
    
    c = onCleanup(@cleanup);
    global multifit_ind x_per_spec
    
    %--- Parse the inputs and set defaults
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true;
    p.addRequired('spec1d',@(x) isa(x,'spec1d'))
    p.addRequired('func',@ischar)
    p.addRequired('pin',@(x) isnumeric(x) | iscell(x))
    p.addRequired('fixed',@(x) isnumeric(x) | islogical(x))
    fcp = [];
    if length(varargin) ~= 1
        p.addParamValue('fcp',[0.0001 20 0.0001],@(x) length(x)==3)
    else
        fcp = varargin{1};
    end
    p.addParamValue('bounds',[-Inf(length(pin),1) Inf(length(pin),1)],@(x) (size(x,1)==length(pin) & size(x,2)==2) | iscell(x))
    p.addParamValue('dfdp','specdfdp_multi',@(x) ischar(x))
    p.addParamValue('confidence',0.05,@(x) isnumeric(x) & length(x)==1 & x>0 & x<1)
    p.addParamValue('inequc',{zeros(length(pin),0) zeros(0,1)},@iscell)
    p.addParamValue('sep',cell(length(s1),1), @iscell)
    p.addParamValue('method','lsquare',@ischar)
    p.addParamValue('window',0,@(x) (isnumeric(x) && length(x)==2) || all(cellfun(@length,x)==2))
    p.addParamValue('parallel',0,@(x) x==0 | x == 1)
    p.addParamValue('criteria','least_square',@(x) ischar(x))
    p.addParamValue('relations',{}, @iscell)
    p.addParamValue('verbose',0,@(x) isnumeric(x) | islogical(x))
    
    if length(varargin) ~= 1
        p.parse(s1,func,pin,notfixed,varargin{:});
    else
        p.parse(s1,func,pin,notfixed);
    end
    options = p.Results;
    
    % Initial parameters
    s1 = options.spec1d;
    func = options.func;
    pin = options.pin;
    notfixed = options.fixed;
    if any(notfixed > 1)
        options.multifit = 1;
        if ~iscell(pin)
            error('The pin must be a cell')
        end
    else
        options.multifit = 0;
        if iscell(pin)
            error('The pin must not be a cell')
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
        warning('Just FYI, this is a multifit!')
        n = length(s1);
        [s1, pin, notfixed] = multifit_ini(s1,pin,notfixed,options.sep);
        options.inequc = {zeros(length(pin),0) zeros(0,1)};
        if iscell(options.bounds)
            options.bounds = multifit_bounds(options.fixed,options.bounds);
        else
            options.bounds = [-Inf(length(pin),1) Inf(length(pin),1)];
        end
    end
    
    %----- Loop over the number of spec1d objects
    sout = spec1d;
    fitdata = struct;
    
    if all([getpref('mtools','experimental'), options.parallel, ~options.multifit])
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
        spmd
            for n = labindex:numlabs:length(s1)
                switch  options.method
                    case 'lsquare'
                        [yfit{n},p,cvg,iter,corp,covp,covr,stdresid,Z,RSq,ra2{n},sig{n}] = speclsqr(all_data{n}{1},all_data{n}{2},all_data{n}{3},pin,notfixed,func,fcp,options);
                    case 'samin'
                        error('This has not been implemented yet, Sorry!')
                    otherwise
                        warning('Fitting method is not implemented or not understood. Using lsquare')
                        [yfit{n},p,cvg,iter,corp,covp,covr,stdresid,Z,RSq,ra2{n},sig{n}] = speclsqr(all_data{n}{1},all_data{n}{2},all_data{n}{3},pin,notfixed,func,fcp,options);
                end
                ChiSq{n} = sum(((all_data{n}{2}-yfit)./all_data{n}{3}).^2 )/(length(all_data{n}{2})-sum(logical(notfixed)));
            end
        end
        for il = 1:length(all_data)
            %----- Get names of fit variables
            if pn
                if options.multifit
                    pnames = {};
                    for i = 1:n
                        [p_new, ~, ind] = multifitp2p(p,notfixed,i);
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
            
            %----- set return
            sloop.x = all_data{il}{1};
            sloop.y = all_data{il}{2};
            sloop.e = all_data{il}{3};
            sloop.x_label = s1(il).x_label;
            sloop.y_label = s1(il).y_label;
            sloop.datafile = s1(il).datafile;
            sloop.yfit = yfit{il};
            sout(il) = spec1d(sloop);
            
            fitdata(il).pvals = p{il};
            fitdata(il).evals = sig{il};
            fitdata(il).function = func;
            fitdata(il).pnames = pnames;
            fitdata(il).chisq = ChiSq{il};
            fitdata(il).rsq = RSq{il};
        end
    else
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
            
            if sum(notfixed)==0
                yfit = feval(func,x,pin);
                p = pin;
                r = corrcoef([y(:),yfit(:)]);
                RSq = r(1,2).^2;
                sig = zeros(size(p));
            else
                %                 %----- Fit data
                %                 switch  options.method
                %                     case 'lsquare'
                %                         [yfit,p,cvg,iter,corp,covp,covr,stdresid,Z,RSq,ra2,sig] = speclsqr(x,y,e,pin,notfixed,func,fcp,options);
                %                     case 'samin'
                %                         error('This has not been implemented yet, Sorry!')
                %                     otherwise
                %                         warning('Fitting method is not implemented or not understood. Using lsquare')
                %                         [yfit,p,cvg,iter,corp,covp,covr,stdresid,Z,RSq,ra2,sig] = speclsqr(x,y,e,pin,notfixed,func,fcp,options);
                %                 end
                % we need to call the optimization method with the eval_criteria as FUN
                % call minimizer ===============================================================
                if abs(nargin(options.optimizer)) == 1 || abs(nargin(options.optimizer)) >= 6
                    [p,yfit,message,output] = feval(options.optimizer, ...
                        @(pars) eval_criteria(model, pars, options.criteria, a, varargin{:}), pars, options, constraints);
                else
                    % Constraints not supported by optimizer
                    [p,yfit,message,output] = feval(options.optimizer, ...
                        @(pars) eval_criteria(model, pars, options.criteria, a, varargin{:}), pars, options);
                end                
            end
            %----- Goodness of fit
            v = length(y)-sum(logical(notfixed));
            ChiSq = sum(((y-yfit)./e).^2 )/v;
            
            %----- Get names of fit variables
            if pn
                if options.multifit
                    pnames = {};
                    for i = 1:n
                        [p_new, ~, ind] = multifitp2p(p,notfixed,i);
                        [~, ~, temp] = feval(func,x(ind),p_new,1);
                        pnames = [pnames(:); temp(:)];
                    end
                else
                    [~,~,pnames] = feval(func,x,p,1);
                end
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
    end
    if options.multifit
        [sout, fitdata] = multifit_extract(sout,fitdata,notfixed);
        x_per_spec = [];
        multifit_ind = [];
    end
end


% format output arguments ======================================================
pars_out = reshape(pars_out, [ 1 numel(pars_out) ]); % row vector
model.ParameterValues = pars_out;
if ~isempty(inputname(1))  
    try
  assignin('caller',inputname(1),model); % update in original object
    end
end

if nargout > 3 || (isfield(options,'Diagnostics') && (strcmp(options.Diagnostics, 'on') || any(options.Diagnostics == 1)))
  output.modelValue = feval(model, pars_out, a.Axes{:});
  index=find(~isnan(a.Signal) & ~isnan(output.modelValue));
  if ~isscalar(a.Error), e = a.Error(index); else e=a.Error; end
  output.corrcoef   = eval_corrcoef(a.Signal(index), e, output.modelValue(index));
  output.residuals  = a.Signal - output.modelValue;
  % Rwp
  Rwp  = sqrt(sum((output.residuals(index)./e).^2)/sum((a.Signal(index)./e).^2));
  output.Rfactor    = Rwp;
  % Rexp
  Rexp = sqrt((length(a.Signal(index)) - length(pars))./sum((a.Signal(index)./e).^2));
  % Goodness of fit
  GoF  = Rwp/Rexp; 
  if strcmp(options.Display, 'iter') | strcmp(options.Display, 'final') | ...
    (isfield(options,'Diagnostics') && (strcmp(options.Diagnostics, 'on') || any(options.Diagnostics == 1)))
    fprintf(1, ' Correlation coefficient=%g (closer to 1 is better)\n',  output.corrcoef);
    fprintf(1, ' Weighted R-factor      =%g (Rwp, smaller that 0.2 is better)\n', output.Rfactor);
    fprintf(1, ' Experimental R-factor  =%g (Rexp)\n', Rexp);
  end
  if abs(output.corrcoef) < 0.6 && ~isscalar(a.Error)
    name = inputname(2);
    if isempty(name)
      name = 'a';
    end
    fprintf(1, ' WARNING: The fit result is BAD. You may improve it by setting %s.Error=1\n',...
      name);
  end
  
  if ~isempty(is_idata)
    % make it an iData
    b = is_idata;
    % fit(signal/monitor) but object has already Monitor -> we compensate Monitor^2
    if not(all(a.Monitor == 1 | a.Monitor == 0)) 
      output.modelValue    = bsxfun(@times,output.modelValue, a.Monitor); 
    end
    setalias(b,'Signal', output.modelValue, model.Name);
    b.Title = [ model.Name '(' char(b) ')' ];
    b.Label = b.Title;
    b.DisplayName = b.Title;
    setalias(b,'Error', 0);
    setalias(b,'Parameters', pars_out, [ model.Name ' model parameters for ' a.Name ]);
    setalias(b,'Model', model, model.Name);
    output.modelValue = b;
  else
    if length(a.Axes) == 1
      output.modelAxes  = a.Axes{1};
    else
      output.modelAxes  = a.Axes(:);
    end
  end
  output.model      = model;
  
  % set output/results
  if ischar(message) | ~isfield(output, 'message')
    output.message = message;
  else
    output.message = [ '(' num2str(message) ') ' output.message ];
  end
  output.parsNames  = model.Parameters;
  % final plot when in OutputFcn mode
  eval_criteria(model, pars_out, options.criteria, a, varargin{:});
end
if ~isempty(pars_isstruct)
  % first rebuild the model parameter structure
  pars_out = cell2struct(num2cell(pars_out(:)'), strtok(model.Parameters(:)'), 2);
  % then add initial additional fields
  if isstruct(pars_isstruct)
    f = fieldnames(pars_isstruct);
    for index=1:length(f)
      pars_out.(f{index}) = pars_isstruct.(f{index});
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EMBEDDED FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
% this way 'options' is available in here...

  function c = eval_criteria(model, p, criteria, a, varargin)

  % criteria to minimize
    if nargin<5, varargin={}; end
    % then get model value
    Model  = feval(model, p, a.Axes{:}, varargin{:}); % return model values

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
      
      if (options.funcCount < 50 && abs(etime(options.updated, clock)) > 0.5) ...
       || abs(etime(options.updated, clock)) > 2
        iFunc_private_fminplot(a,model,p,Model,options,c)
      end
    end
    
  end % eval_criteria (embedded)

end % fits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=eval_corrcoef(Signal, Error, Model)
% correlation coefficient between the data and the model
  
  % compute the correlation coefficient
  if isempty(Error) || isscalar(Error) || all(Error(:) == Error(end))
    wt = 1;
  else
    wt = 1./Error;
    wt(find(~isfinite(wt))) = 0;
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