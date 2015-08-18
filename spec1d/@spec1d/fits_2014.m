function [sout,fitdata]=fits(varargin)
    % function [sout,fitdata]=fits(s1,func,pin,notfixed,fcp,win)
    
    %
    % [sout,fitdata]=fits(s1,func,pin,notfixed,fcp,win)
    %
    % @SPEC1D/FIT Fits data in spec1d object s1 to MFIT function
    % specified in func. Also works for an array of spec1d objects.
    %
    % Notes:
    % 1. If pin='auto' (or is ommitted) only s1 and func are specified,
    %    and func is a simple peak shape,
    %    then an attempt is made to guess the start parameters.
    % 2. Retrun parameters are:
    %       r - original spectrum with fitted curve
    %       fitdata - structure with fields:
    %          fitdata.pvals - parameter values
    %          fitdata.evals - error values
    %          fitdata.pnames- parameter names
    %          fitdata.chisq - Chi^2 of fit
    %
    % Version 2.0, February 2001
    % Des McMorrow and Henrik Ronnow
    % Modified by Simon Ward.
    % Now with a better fitting routine and R^2 calculation
    
    %--- Set defaults
    func=varargin{2};
    pn=1;
    if ischar(func)
        if strcmp(func(1),'@')
            func=str2func(func);
            pn=0;
        end
    end
    if nargin==2,
        pin='auto';
        notfixed=ones(1,4);
    else
        pin=varargin{3};
    end
    if nargin==3,
        notfixed=ones(size(pin));
    elseif nargin>3
        notfixed=varargin{4};
    end
    if nargin<=4 || isa(varargin{5},'struct')
        fcp=[0.0001 20 0.0001];
    else
        fcp=varargin{5};
    end
    
    if isempty(pin), pin='auto'; end
    if isempty(notfixed), notfixed=ones(1,4); end
    
    options=[];
    j=1;
    for i=1:nargin
        if isa(varargin{i},'struct')
            options=varargin{i};
        end
        if isa(varargin{i},'spec1d')
            s1=varargin{i};
        end
    end
    
    options.dfdp = 'specdfdp_multi';
    %----- Loop over the number of spec1d objects
    
    sout=spec1d;
    fitdata=[];
    pinin=pin;
    
    global s_sep
    
    for il=1:length(s1)
        
        %----- Check for auto guess for simple peaks
        
        peaklist=strvcat('gauss','gauss2','lorz','lorz2','gauss_area');
        if ischar(pinin) %| nargin==2
            if ~isempty(strmatch(func,peaklist))
                stats=peakm(s1(il),1);
                pin=[stats(1) stats(2) stats(3) stats(5)];
            end
        end
        
        x=s1(il).x; y=s1(il).y; e=s1(il).e;
        
        %----- Check for window to remove points
        
        %    if nargin==6
        %      in=[];
        %      for n=1:size(win,1)
        %        in=[in x>win(n,1) & x<win(n,2)];
        %      end
        %      x=x(in);
        %      y=y(in);
        %      e=e(in);
        %    end
        
        %----- Remove zeros from e
        %    ezeros=find(e==0);
        x(e==0 | isinf(e) | isnan(e))=[];
        y(e==0 | isinf(e) | isnan(e))=[];
        e(e==0 | isinf(e) | isnan(e))=[];
        %    x(ezeros)=[]; y(ezeros)=[]; e(ezeros)=[];
        %----- Fit data
%         if isempty(s_sep)
            [yfit,p,cvg,iter,corp,covp,covr,stdresid,Z,RSq,RaSq,sig]=speclsqr(x,y,e,pin,notfixed,func,fcp,options);
%         else
%             for i = 1:length(s_sep)
%                 [pout dout] = p2plocal(pin,notfixed,s_sep);
%                 [yfit,p,cvg,iter,corp,covp,covr,stdresid,Z,RSq,RaSq,sig]=speclsqr(x,y,e,pout,dout,func,fcp,options);
%             end
%         end
        %----- Get names of fit variables
        
        v=length(y)-length(find(notfixed));
        ChiSq = sum(((y-yfit)./e).^2 )/v;
        
        if pn
            [dummy,dummy,pnames]=feval(func,x,p,1);
        else
            for i=1:length(p)
                pnames{i}=['p' num2str(i)];
            end
        end
        
        %----- set return
        sloop.x=x;
        sloop.y=y;
        sloop.e=e;
        sloop.x_label=s1(il).x_label;
        sloop.y_label=s1(il).y_label;
        sloop.datafile=s1(il).datafile;
        sloop.yfit=yfit;
        
        sout(il)=spec1d(sloop);
        
        fitdata(il).pvals=p;
        fitdata(il).evals=sig;
        fitdata(il).function = func;
        fitdata(il).pnames=pnames;
        fitdata(il).chisq=ChiSq;
        fitdata(il).rsq=RSq;
    end