function [sout,fitdata]=fits(s1,func,pin,notfixed,varargin)
    %% function [sout,fitdata]=fits(s1,func,pin,notfixed,fcp,win)
    %
    % [sout,fitdata]=fits(s1,func,pin,notfixed,fcp,win)
    %
    % @SPEC1D/FIT Fits data in spec1d object s1 to MFIT function
    % specified in func. Also works for an array of spec1d objects.
    %
    % Version 4.2, July 2015
    % Simon Ward, based on the work of Des McMorrow and Henrik Ronnow
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
%     p.addParamValue('sep',cell(length(s1),1), @iscell)
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
        [s1, pin, notfixed] = multifit_ini(s1,pin,notfixed);
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
                %----- Fit data
                switch  options.method
                    case 'lsquare'
                        [yfit,p,cvg,iter,corp,covp,covr,stdresid,Z,RSq,ra2,sig] = speclsqr(x,y,e,pin,notfixed,func,fcp,options);
                    case 'samin'
                        error('This has not been implemented yet, Sorry!')
                    otherwise
                        warning('Fitting method is not implemented or not understood. Using lsquare')
                        [yfit,p,cvg,iter,corp,covp,covr,stdresid,Z,RSq,ra2,sig] = speclsqr(x,y,e,pin,notfixed,func,fcp,options);
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
                    if nargin(func)<3
                        pnames = num2str((1:numel(p))','p%d');
                    else
                        [~,~,pnames] = feval(func,x,p,1);
                    end
                end
            else
            pnames = num2str((1:numel(p))','p%d');

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

function cleanup
    global multifit_ind x_per_spec
    x_per_spec = [];
    multifit_ind = [];
end

