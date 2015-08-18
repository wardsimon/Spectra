function [sout,fitdata]=fits2(varargin)
%
% [sout,fitdata]=fits(s1,func,pin,notfixed,fcp,win)
%
% @SPEC1D/FIT Fits data in spec1d object s1 to MFIT function
% specified in func. Also works for an array of spec1d objects.
%
% Simon Ward 2012.
% set maxima and minima by using [..., 'min', [], 'max', [] ].
% Either maxima, minima, both or none can be set.
% Backwards compatible with fits routine.
% FIXED Chi^2 that now WORKS!

%--- Set defaults

if ~areTheseToolboxesInstalled('Global Optimization Toolbox')
    warning('Global Optimization Toolbox must be installed, defaulting to fits routine!')
    j=1;
    skip=0;
    for i=1:length(varargin)
        if isa(varargin{i},'char') && i>2
            skip=1;
        elseif ~skip
            nvarargin{j}=varargin{i};
            j=j+1;
            skip=0;
        end
    end
    [sout,fitdata]=fits(nvarargin{:});
else
    s1=varargin{1};
    func=varargin{2};
    if length(varargin)==2
        pin='auto';
        notfixed=ones(1,4);
    else
        pin=varargin{3};
    end
    if length(varargin)>=5
        if length(varargin{5})==3
            fcp=varargin{5};
        end
    else
        fcp=[0.0001 50 0.0001];
    end
    lb=[];
    ub=[];
    for i=1:length(varargin)
        temp=varargin{i};
        if strcmp(temp,'min')
            lb=varargin{i+1};
        end
        if strcmp(temp,'max')
            ub=varargin{i+1};
        end
    end
    if isempty(lb)
        lb=-Inf(1,length(pin));
    end
    if isempty(ub)
        ub=Inf(1,length(pin));
    end
    
    if length(varargin)>=4
        if length(varargin{4})==length(varargin{3}),
            notfixed=varargin{4};
            notfixed=cast(notfixed,'logical');
            fixed=~notfixed;
        end
    else
        notfixed=ones(size(pin));
        fixed=~notfixed;
    end
    %----- Loop over the number of spec1d objects
    
    sout=spec1d;
    fitdata=[];
    pinin=pin;
    for il=1:length(s1)
        x=getfield(s1(il),'x'); y=getfield(s1(il),'y'); e=getfield(s1(il),'e');  wt=1./e(:);
        
        %----- Remove zeros from e
        ezeros=find(e==0);
        x(ezeros)=[]; y(ezeros)=[]; e(ezeros)=[]; wt(ezeros)=[];
        
        %----- Fit data
        infun = @(x,xin) feval(func,xin,x);
        options=optimset('MaxIter',fcp(2),'Display','Off');
        [p2,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit( @(p1, xdata) infun( interlace(pinin(:), p1(:), fixed), xdata),...
        pinin(notfixed), x, y, lb(notfixed), ub(notfixed),options);
        if exitflag < 0
            warning('Convergence not achieved')
        end
        p2 = interlace( pinin(:), p2(:), fixed);
        
        Qinv=diag(wt.*wt);
        Q=diag((0*wt+1)./(wt.^2));
        m=length(y);
        n=sum(notfixed);
        covr=residual'*Qinv*residual*Q/(m-n);                 %covariance of residuals
        Vy=1/(1-n/m)*covr;  % Eq. 7-13-22, Bard         %covariance of the data
        jtgjinv=pinv(jacobian'*Qinv*jacobian);
        covp=jtgjinv*jacobian'*Qinv*Vy*Qinv*jacobian*jtgjinv; % Eq. 7-5-13, Bard %cov of parm est
        for k=1:n,
            for j=k:n,
                corp(k,j)=covp(k,j)/sqrt(abs(covp(k,k)*covp(j,j)));
                corp(j,k)=corp(k,j);
            end;
        end;
        
        std=sqrt(diag(covp));
        
        j=1;
        dp=-notfixed*fcp(1);
        sig=zeros(size(p2));
        for i=1:length(std)
            while dp(j)==0
                j=j+1;
            end
            sig(j)=std(i);
            j=j+1;
        end
        
        
        %----- Fit function values
        yfit=feval(func,x,p2);
        
        ChiSq=sum(((y-yfit).^2)/sum((y-mean(y)).^2))*((length(y)-1)/(length(y)-sum(notfixed)+1));
        if areTheseToolboxesInstalled('Statistics Toolbox')
            chpdf=chi2pdf(ChiSq,sum(notfixed));
            if chpdf > 0.05
                warning('MATLAB:poorFit','This fit probably not valid!\nThe pdf is %d. the observed deviation from the null hypothesis is significant',chpdf);
            end
        end
        %----- Get names of fit variables
        [dummy,dummy,pnames]=feval(func,x,p2,1);
        
        %----- set return
        sloop.x=x;
        sloop.y=y;
        sloop.e=e;
        sloop.x_label=getfield(s1(il),'x_label');
        sloop.y_label=getfield(s1(il),'y_label');
        sloop.datafile=getfield(s1(il),'datafile');
        sloop.yfit=yfit;
        
        sout(il)=spec1d(sloop);
        
        fitdata(il).pvals=p2;
        fitdata(il).evals=sig;
        fitdata(il).pnames=pnames;
        fitdata(il).chisq=ChiSq;
    end
end

end

function a = interlace( a, x, fix )
a(~fix) = x;
end

function tf = areTheseToolboxesInstalled(requiredToolboxes)
%ARETHESETOOLBOXESINSTALLED takes a cell array of toolbox names and checks whether they are currently installed
% SYNOPSIS tf = areTheseToolboxesInstalled(requiredToolboxes)
%
% INPUT requiredToolboxes: cell array with toolbox names to test for. Eg.
%        {'MATLAB','Image Processing Toolbox'}
%
% OUTPUT tf: true or false if the required toolboxes are installed or not
%%%%%%%%%%%%%%%%%%%%%%%%%%

% get all installed toolbox names
v = ver;
[installedToolboxes{1:length(v)}] = deal(v.Name);
tf = all(ismember(requiredToolboxes,installedToolboxes));
end