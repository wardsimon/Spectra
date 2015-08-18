function varargout = confidence(s,fit,conf,varargin)
    %% Confidence
    %  Calculates the confidence for a given fit and interval.
    %  -------------------------------------------------------
    %  !!! Input !!!
    %  s    = spec1d object which has a fit attached
    %  fit  = results from spec1d's fitting procedure
    %  conf = confidence level (0 < conf < 1) OR n*sigma where n >=1
    %  OPTIONAL
    %  x   = 
    %  -------------------------------------------------------
    %  !!! Output !!!
    %  If no output vairables, plots to current figure
    %  If 1 output vairable, gives a column vector with [ub lb]
    %  If 2 output vairables, gives 2 column vectors ub and lb
    %
    % Simon Ward 05/08/2013
    % simon.ward@psi.ch
    
    for i = 1:length(s)
    
    if nargin == 4
        x = varargin{1};
    else
        x = s(i).x;
    end
    
    y0 = s(i).y;
    yf = s(i).yfit;
    
    if isempty(yf)
        error('Fit is needed!')
    end
    
    if conf < 0 || conf >= 1
        if conf >= 1
            conf = 1-erf(conf/sqrt(2));
        else
            error('Confidence level must be between 0 and 0.9999, or sigma >= 1')
        end
    end
    
    fprintf('Confidence level: %02.04f%% == %01.02f sigma. \n',100*(1-conf),sqrt(2)*erfinv(1-conf))

    resid = y0-yf;
    jac = feval(@dfdp,s(i).x,yf,fit(i).pvals,fit(i).evals,fit(i).function);
    jac(:,fit(i).evals==0)=[];
    
    p = fit(i).pvals(fit(i).evals~=0);
    
    delta=CalcConfidence(p,resid,jac,conf);
    
    p0 = fit(i).pvals;
    p1 = p0;
    p1(fit(i).evals~=0) = p0(fit(i).evals~=0) + delta;
    p2 = p0;
    p2(fit(i).evals~=0) = p0(fit(i).evals~=0) - delta;
    
    
    y1=feval(fit(i).function,x,p1);
    y2=feval(fit(i).function,x,p2);
    
    f = get(0,'CurrentFigure');
    
    if nargout == 1
        varargout{1}(:,1) = y1;
        varargout{1}(:,2) = y2;
    elseif nargout ==2
        varargout{1} = y1;
        varargout{2} = y2;
    else
        if ~isempty(f)
            figure(f)
            l = findobj(gcf,'Type','Line','-and','LineStyle','-','-not','Color',[.5 .5 .5]);
            hold on
            plot(x,y1,'Color',get(l(end-i+1),'Color'),'LineWidth',2,'LineStyle','--')
            plot(x,y2,'Color',get(l(end-i+1),'Color'),'LineWidth',2,'LineStyle','--')
        else
            error('You must give me a plot or a vairable to assign!')
        end
    end
    end
end



function delta = CalcConfidence(beta,resid,J,alpha)
    %   Confidence intervals for parameters in nonlinear regression.
    %   beta = pin;
    %   resid = y0-yf;
    %   J = jaccobian
    %   alpha = confidence bound.
    %   delta = confidence(pin,residules,jacobian,confidence)
    
    if ~isreal(beta) || ~isreal(J)
        error(message('Complex Jacobian or beta'));
    end
    
    % Remove missing values.
    resid = resid(:);
    missing = isnan(resid);
    if ~isempty(missing)
        resid(missing) = [];
    end
    n = length(resid);
    p = numel(beta);
    v = n-p;
    
    % Estimate covariance from J and residuals
    J(missing,:) = [];
    if size(J,1)~=n || size(J,2)~=p
        error('Size of Jacobian is wrong!');
    end
    
    % Approximation when a column is zero vector
    temp = find(max(abs(J)) == 0);
    if ~isempty(temp)
        J(:,temp) = sqrt(eps(class(J)));
    end
    
    % Calculate covariance matrix
    [~,R] = qr(J,0);
    Rinv = R\eye(size(R));
    diag_info = sum(Rinv.*Rinv,2);
    
    rmse = norm(resid) / sqrt(v);
    se = sqrt(diag_info) * rmse;
    
    % Calculate confidence interval
    delta = se * tinv(1-alpha/2,v);
    
end

function x=tinv(p,v)
    % Initialize Y to zero, or NaN for invalid d.f.
    if isa(p,'single') || isa(v,'single')
        x = NaN(size(p),'single');
    else
        x = NaN(size(p));
    end
    
    % The inverse cdf of 0 is -Inf, and the inverse cdf of 1 is Inf.
    x(p==0 & v > 0) = -Inf;
    x(p==1 & v > 0) = Inf;
    
    k0 = (0<p & p<1) & (v > 0);
    
    % Invert the Cauchy distribution explicitly
    k = find(k0 & (v == 1));
    if any(k)
        x(k) = tan(pi * (p(k) - 0.5));
    end
    
    % For small d.f., call betaincinv which uses Newton's method
    k = find(k0 & (v < 1000));
    if any(k)
        q = p(k) - .5;
        df = v(k);
        t = (abs(q) < .25);
        z = zeros(size(q),class(x));
        oneminusz = zeros(size(q),class(x));
        if any(t)
            % for z close to 1, compute 1-z directly to avoid roundoff
            oneminusz(t) = betaincinv(2.*abs(q(t)),0.5,df(t)/2,'lower');
            z(t) = 1 - oneminusz(t);
        end
        t = ~t; % (abs(q) >= .25);
        if any(t)
            z(t) = betaincinv(2.*abs(q(t)),df(t)/2,0.5,'upper');
            oneminusz(t) = 1 - z(t);
        end
        x(k) = sign(q) .* sqrt(df .* (oneminusz./z));
    end
    
    % For large d.f., use Abramowitz & Stegun formula 26.7.5
    % k = find(p>0 & p<1 & ~isnan(x) & v >= 1000);
    k = find(k0 & (v >= 1000));
    if any(k)
        xn = norminv(p(k));
        df = v(k);
        x(k) = xn + (xn.^3+xn)./(4*df) + ...
            (5*xn.^5+16.*xn.^3+3*xn)./(96*df.^2) + ...
            (3*xn.^7+19*xn.^5+17*xn.^3-15*xn)./(384*df.^3) +...
            (79*xn.^9+776*xn.^7+1482*xn.^5-1920*xn.^3-945*xn)./(92160*df.^4);
    end
end