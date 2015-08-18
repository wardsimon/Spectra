function [f,p,cvg,iter,corp,covp,covr,stdresid,Z,r2,ra2,std] = speclsqr(x,y,err,pin,dpin,F,fcp,options)
    %%
    %% Levenberg-Marquardt nonlinear regression of f(x,p) to y(x).
    %%
    %function [f,p,cvg,iter,corp,covp,covr,stdresid,Z,r2]=
    %                   leasqr(x,y,pin,F,{stol,niter,wt,dp,dFdp,options})
    %
    % Version 3.beta http://octave-optim.sourcearchive.com/documentation/1.0.12-1/files.html
    % Optional parameters are in braces {}.
    % x = vector or matrix of independent variables, 1 entry or row per
    %   observation.
    % y = vector of observed values, same length as x or as number of
    %   rows of x.
    % wt = vector (dim=length(y)) of statistical weights.  These should
    %   be set to be proportional to (sqrt of var(y))^-1; (That is, the
    %   covariance matrix of the data is assumed to be proportional to
    %   diagonal with diagonal equal to (wt.^2)^-1.  The constant of
    %   proportionality will be estimated.); default = ones(length(y),1).
    % pin = vec of initial parameters to be adjusted by leasqr.
    % dp = fractional increment of p for numerical partial derivatives;
    %   default = .001*ones(size(pin))
    %   dp(j) > 0 means central differences on j-th parameter p(j).
    %   dp(j) < 0 means one-sided differences on j-th parameter p(j).
    %   dp(j) = 0 holds p(j) fixed i.e. leasqr wont change initial guess: pin(j)
    % F = name of function in quotes or function handle; the function
    %   shall be of the form y=f(x,p), with y, x, p of the form y, x, pin
    %   as described above; the returned y must be a column vector.
    % dFdp = name of partial derivative function in quotes; default is 'dfdp', a
    %   slow but general partial derivatives function; the function shall be
    %   of the form prt=dfdp(x,f,p,dp,F[,bounds]). For backwards
    %   compatibility, the function will only be called with an extra
    %   'bounds' argument if the 'bounds' option is explicitely specified
    %   to leasqr (see dfdp.m).
    % stol = scalar tolerance on fractional improvement in scalar sum of
    %   squares = sum((wt.*(y-f))^2); default stol = .0001;
    % niter = scalar maximum number of iterations; default = 20;
    % options = structure, currently recognized fields are 'fract_prec',
    %   'max_fract_change', 'inequc', and 'bounds'. For backwards compatibility,
    %   'options' can also be a matrix whose first and second column
    %   contains the values of 'fract_prec' and 'max_fract_change',
    %   respectively.
    %   Field 'options.fract_prec': column vector (same length as 'pin') of
    %   desired fractional precisions in parameter estimates. Iterations
    %   are terminated if change in parameter vector (chg) on two
    %   consecutive iterations is less than their corresponding elements in
    %   'options.fract_prec'.  [ie. all(abs(chg*current parm est) <
    %   options.fract_prec) on two consecutive iterations.], default =
    %   zeros().
    %   Field 'options.max_fract_change': column vector (same length as
    %   'pin) of maximum fractional step changes in parameter vector.
    %   Fractional change in elements of parameter vector is constrained to
    %   be at most 'options.max_fract_change' between sucessive iterations.
    %   [ie. abs(chg(i))=abs(min([chg(i)
    %   options.max_fract_change(i)*current param estimate])).], default =
    %   Inf*ones().
    %   Field 'options.inequc': cell-array containing a matrix (say m)
    %   and a column vector (say v), specifying linear inequality
    %   constraints of the form `m.' * parameters + v >= 0'. If the
    %   constraints are just bounds, it is suggested to specify them in
    %   'options.bounds' instead, since then some sanity tests are
    %   performed, and since the function 'dfdp.m' is guarantied not to
    %   violate constraints during determination of the numeric gradient
    %   only for those constraints specified as 'bounds'.
    %   Field 'options.bounds': two-column-matrix, one row for each
    %   parameter in 'pin'. Each row contains a minimal and maximal value
    %   for each parameter. Default: [-Inf, Inf] in each row. If this
    %   field is used with an existing user-side function for 'dFdp'
    %   (see above) the functions interface might have to be changed.
    %   _Warning_: If constraints (or bounds) are set, returned guesses
    %   of corp, covp, and Z are generally invalid, even if no constraints
    %   are active for the final parameters.
    %
    %          OUTPUT VARIABLES
    % f = column vector of values computed: f = F(x,p).
    % p = column vector trial or final parameters. i.e, the solution.
    % cvg = scalar: = 1 if convergence, = 0 otherwise.
    % iter = scalar number of iterations used.
    % corp = correlation matrix for parameters.
    % covp = covariance matrix of the parameters.
    % covr = diag(covariance matrix of the residuals).
    % stdresid = standardized residuals.
    % Z = matrix that defines confidence region (see comments in the source).
    % r2 = coefficient of multiple determination, intercept form.
    %
    % Not suitable for non-real residuals.
    
    % The following two blocks of comments are chiefly from the original
    % version for Matlab. For later changes the logs of the Octave Forge
    % svn repository should also be consulted.
    
    % A modified version of Levenberg-Marquardt
    % Non-Linear Regression program previously submitted by R.Schrager.
    % This version corrects an error in that version and also provides
    % an easier to use version with automatic numerical calculation of
    % the Jacobian Matrix. In addition, this version calculates statistics
    % such as correlation, etc....
    %
    % Version 3 Notes
    % Errors in the original version submitted by Shrager (now called
    % version 1) and the improved version of Jutan (now called version 2)
    % have been corrected.
    % Additional features, statistical tests, and documentation have also been
    % included along with an example of usage.  BEWARE: Some the the input and
    % output arguments were changed from the previous version.
    %
    %     Ray Muzic     <rfm2@ds2.uh.cwru.edu>
    %     Arthur Jutan  <jutan@charon.engga.uwo.ca>
    
    % Richard I. Shrager (301)-496-1122
    % Modified by A.Jutan (519)-679-2111
    % Modified by Ray Muzic 14-Jul-1992
    %       1) add maxstep feature for limiting changes in parameter estimates
    %          at each step.
    %       2) remove forced columnization of x (x=x(:)) at beginning. x
    %          could be a matrix with the ith row of containing values of
    %          the independent variables at the ith observation.
    %       3) add verbose option
    %       4) add optional return arguments covp, stdresid, chi2
    %       5) revise estimates of corp, stdev
    % Modified by Ray Muzic 11-Oct-1992
    %  1) revise estimate of Vy.  remove chi2, add Z as return values
    %       (later remark: the current code contains no variable Vy)
    % Modified by Ray Muzic 7-Jan-1994
    %       1) Replace ones(x) with a construct that is compatible with versions
    %          newer and older than v 4.1.
    %       2) Added global declaration of verbose (needed for newer than v4.x)
    %       3) Replace return value var, the variance of the residuals
    %          with covr, the covariance matrix of the residuals.
    %       4) Introduce options as 10th input argument.  Include
    %          convergence criteria and maxstep in it.
    %       5) Correct calculation of xtx which affects coveraince estimate.
    %       6) Eliminate stdev (estimate of standard deviation of
    %          parameter estimates) from the return values.  The covp is a
    %          much more meaningful expression of precision because it
    %          specifies a confidence region in contrast to a confidence
    %          interval.. If needed, however, stdev may be calculated as
    %          stdev=sqrt(diag(covp)).
    %       7) Change the order of the return values to a more logical order.
    %       8) Change to more efficent algorithm of Bard for selecting epsL.
    %       9) Tighten up memory usage by making use of sparse matrices (if
    %          MATLAB version >= 4.0) in computation of covp, corp, stdresid.
    % Added linear inequality constraints with quadratic programming to
    % this file and special case bounds to this file and to dfdp.m, Olaf
    % Till 24-Feb-2010
    %
    % References:
    % Bard, Nonlinear Parameter Estimation, Academic Press, 1974.
    % Draper and Smith, Applied Regression Analysis, John Wiley and Sons, 1981.
    
    %% argument processing
    %%
    global verbose x_per_spec
    if isempty(verbose)
        verbose=0;       %This will not tell them the results
    end
    %     dFdp='dfdp';
    dFdp = options.dfdp;
    
    stol=fcp(3);
    niter=fcp(2);
    dp=-(dpin>0)*fcp(1);
    
    %%
    wt=1./err(:);
    y=y(:); wt=wt(:); pin=pin(:); dp=dp(:); %change all vectors to columns
    if (isvector (x))
        x = x(:);
    end
    %% check data vectors- same length?
    m=length(y); n=length(pin);
    if (size (x, 1) ~= m)
        error('input(x)/output(y) data must have same length ')
    end
    
    %% processing of 'options'
    pprec = zeros (n, 1);
    maxstep = Inf * ones (n, 1);
    mc = zeros (n, 0);
    vc = zeros (0, 1);
    bounds = cat (2, -Inf * ones (n, 1), Inf * ones (n, 1));
    if (nargin > 6)
        if ~isempty(options)
            if (isfield(options,'confidence'))
                conf_p=1-options.confidence;
            else
                conf_p=0.05;
            end
            if (isfield (options, 'fract_prec'))
                pprec = options.fract_prec;
                if (rows (pprec) ~= n || columns (pprec) ~= 1)
                    error ('fractional precisions: wrong dimensions');
                end
            end
            if (isfield (options, 'max_fract_change'))
                maxstep = options.max_fract_change;
                if (rows (maxstep) ~= n || columns (maxstep) ~= 1)
                    error ('maximum fractional step changes: wrong dimensions');
                end
            end
            if (isfield (options, 'inequc'))
                mc = options.inequc{1};
                vc = options.inequc{2};
                [rm, cm] = size (mc);
                [rv, cv] = size (vc);
                if (rm ~= n || cm ~= rv || cv ~= 1)
                    error ('inequality constraints: wrong dimensions');
                end
                if (any (mc.' * pin + vc < 0))
                    error ('initial parameters violate inequality constraints');
                end
            end
            if (isfield (options, 'bounds'))
                bounds = options.bounds;
                if (size(bounds,1) ~= n || size(bounds,2) ~= 2)
                    error ('bounds: wrong dimensions');
                end
                idx = bounds(:, 1) > bounds(:, 2);
                tp = bounds(idx, 2);
                bounds(idx, 2) = bounds(idx, 1);
                bounds(idx, 1) = tp;
                idx = bounds(:, 1) == bounds(:, 2);
                if (any (idx))
                    warning ('MATLAB_Fit_Error:LB_UB',...
                        'lower and upper bounds identical for some parameters, setting the respective elements of dp to zero');
                    dp(idx) = 0;
                end
                idx = pin < bounds(:, 1);
                if (any (idx))
                    warning ('MATLAB_Fit_Error:LB','some initial parameters set to lower bound');
                    pin(idx) = bounds(idx, 1);
                end
                idx = pin > bounds(:, 2);
                if (any (idx))
                    warning ('MATLAB_Fit_Error:UB','some initial parameters set to upper bound');
                    pin(idx) = bounds(idx, 2);
                end
                tp = eye (n);
                lidx = ~ isinf (bounds(:, 1));
                uidx = ~ isinf (bounds(:, 2));
                mc = cat (2, mc, tp(:, lidx), - tp(:, uidx));
                vc = cat (1, vc, - bounds(lidx, 1), bounds(uidx, 2));
            end
        else
            conf_p=0.05;
            bounds(:,1)=-Inf(1,length(pin));
            bounds(:,2)= Inf(1,length(pin));
        end
    end
    
    if (all (dp == 0))
        error ('No free parameters');
    end
    
    %% set up for iterations
    %%
    p = pin;
    if isempty(x_per_spec)
        f=feval(F,x,p);
    else
        f = zeros(sum(x_per_spec),1);
        for i = 1:length(x_per_spec)
            [p_new, d_new, ind] = multifitp2p(p,dpin,i);
            f(ind) = feval(F,x(ind),p_new);
        end
    end
    fbest=f;
    pbest=p;
    r=wt.*(y-f);
    if (~isreal (r))
        error ('Weighted residuals are not real');
    end
    ss = r.' * r;
    sbest=ss;
    chgprev=Inf*ones(n,1);
    cvg=0;
    epsLlast=1;
    epstab=[.1, 1, 1e2, 1e4, 1e6];
    
    %% for testing
    %% new_s = false;
    %% if (isfield (options, 'new_s'))
    %%   new_s = options.new_s;
    %% end
    
    nz = eps; % This is arbitrary. Constraint fuction will be regarded as
    % <= zero if less than nz.
    %% do iterations
    %%
    for iter=1:niter
        c_act = mc.' * p + vc < nz; % index of active constraints
        mca = mc(:, c_act);
        %     vca = vc(c_act);
        mcat = mca.';
        nrm = zeros (1, n);
        pprev=pbest;
        prt = feval(dFdp,x,fbest,p,dp,F,bounds);
        if any(isnan(prt(:)));
            error('Functional pin differential is NaN. Somethings gone wrong!')
        end
        r=wt.*(y-fbest);
        if (~isreal (r))
            error ('weighted residuals are not real');
        end
        sprev=sbest;
        sgoal=(1-stol)*sprev;
        msk = dp ~= 0;
        prt(:, msk) = prt(:, msk) .* wt(:, ones (1, sum (msk)));
        nrm(msk) = sumsq (prt(:, msk), 1);
        msk = nrm > 0;
        nrm(msk) = 1 ./ sqrt (nrm(msk));
        prt = prt .* nrm(ones (1, m), :);
        nrm = nrm.';
        % Matts change svd(prt,0) was the original. Now we can fit more
        % parameters than points!
        
        [prt,s,v]=svd(prt,'econ');
        s=diag(s);
        g = prt.' * r;
        for jjj=1:length(epstab)
            epsL = max(epsLlast*epstab(jjj),1e-7);
            %% printf ('epsL: %e\n', epsL); % for testing
            
            %% Usage of this 'ser' later is equivalent to pre-multiplying the
            %% gradient with a positive-definit matrix, but not with a
            %% diagonal matrix, at epsL -> Inf; so there is a fallback to
            %% gradient descent, but not in general to descent for each
            %% gradient component. Using the commented-out 'ser' ((1 / (1 +
            %% epsL^2)) * (1 ./ se + epsL * s)) would be equivalent to using
            %% Marquardts diagonal of the Hessian-approximation for epsL ->
            %% Inf, but currently this gives no advantages in tests, even with
            %% constraints.
            ser = 1 ./ sqrt((s.*s)+epsL);
            %% se=sqrt((s.*s)+epsL);
            %%if (new_s)
            %% %% for testing
            %% ser = (1 / (1 + epsL^2)) * (1 ./ se + epsL * s);
            %% else
            %% ser = 1 ./ se;
            %% end
            tp1 = (v * (g .* ser)) .* nrm;
            if (any (c_act))
                %% calculate chg by 'quadratic programming'
                idx = ones (1, size (mca, 2));
                nrme = nrm(:, idx);
                ser2 = ser .* ser;
                tp2 = nrme .* (v * (ser2(:, idx) .* (v.' * (nrme .* mca))));
                [lb, idx] = cpiv (mcat * tp1, mcat * tp2);
                chg = tp1 + tp2(:, idx) * lb;
                %% collect inactive constraints
                mcit = mc(:, ~ c_act).';
                vci = vc(~ c_act);
            else
                %% chg is the Levenberg/Marquardt step
                chg = tp1;
                %% inactive constraints consist of all constraints
                mcit = mc.';
                vci = vc;
            end
            %% apply inactive constraints (since this is a Levenberg/Marquardt
            %% algorithm, no line-search is performed here)
            hstep = mcit * chg;
            idx = hstep < 0;
            if (any (idx))
                k = min (1, min (- (vci(idx) + mcit(idx, :) * pprev) ./ ...
                    hstep(idx)));
                chg = k * chg;
            end
            %% check the maximal stepwidth and apply as necessary
            ochg=chg;
            idx = ~isinf(maxstep);
            limit = abs(maxstep(idx).*pprev(idx));
            chg(idx) = min(max(chg(idx),-limit),limit);
            if (verbose && any(ochg ~= chg))
                disp(['Change in parameter(s): ', ...
                    sprintf('%d ',find(ochg ~= chg)), 'maximal fractional stepwidth enforced']);
            end
            aprec=abs(pprec.*pbest);       %---
            %% ss=scalar sum of squares=sum((wt.*(y-f))^2).
            if (any(abs(chg) > 0.1*aprec))%---  % only worth evaluating
                % function if there is some non-miniscule
                % change
                p=chg+pprev;
                if isempty(x_per_spec)
                    f=feval(F,x,p);
                else
                    f = zeros(sum(x_per_spec),1);
                    for i = 1:length(x_per_spec)
                        [p_new, ~, ind] = multifitp2p(p,dpin,i);
                        f(ind) = feval(F,x(ind),p_new);
                    end
                end
                r=wt.*(y-f);
                if (~isreal (r))
                    error ('weighted residuals are not real');
                end
                ss = r.' * r;
                if (ss<sbest)
                    pbest=p;
                    fbest=f;
                    sbest=ss;
                end
                if (ss<=sgoal)
                    break;
                end
            end                          %---
        end
        %% printf ('epsL no.: %i\n', jjj); % for testing
        epsLlast = epsL;
        if (verbose)
            feval(dFdp,x,fbest,p,dp,F,bounds);
        end
        if (ss<eps)
            break;
        end
        aprec=abs(pprec.*pbest);
        %% [aprec, chg, chgprev]
        if (all(abs(chg) < aprec) && all(abs(chgprev) < aprec))
            cvg=1;
            if (verbose)
                fprintf('Parameter changes converged to specified precision\n');
            end
            break;
        else
            chgprev=chg;
        end
        if (ss>sgoal)
            break;
        end
    end
    
    %% set return values
    %%
    p=pbest;
    f=fbest;
    % ss=sbest;
    cvg=((sbest>sgoal)|(sbest<=eps)|cvg);
    
    resid=y-f;
    resid=resid(:);
    
    if (cvg ~= 1)
        %% runs test according to Bard. p 201.
        n1 = sum (resid > 0);
        n2 = sum (resid < 0);
        nrun=sum(abs(diff(resid > 0)))+1;
        fprintf('CONVERGENCE NOT ACHIEVED! \n')
        if ((n1 > 10) && (n2 > 10)) % sufficent data for test?
            zed=(nrun-(2*n1*n2/(n1+n2)+1)+0.5)/(2*n1*n2*(2*n1*n2-n1-n2)...
                /((n1+n2)^2*(n1+n2-1)));
            if (zed < 0)
                prob = erfc(-zed/sqrt(2))/2*100;
                fprintf('%0.03f%% chance of fewer than %i runs. \nTry changing fcp... \n',prob,nrun);
            else
                prob = erfc(zed/sqrt(2))/2*100;
                fprintf('%0.03f%% chance of greater than %i runs. \nTry changing fcp... \n',prob,nrun);
            end
        end
    end
    
    if (~(verbose || nargout > 4))
        return
    end
    
    %% CALC VARIANCE COV MATRIX AND CORRELATION MATRIX OF PARAMETERS
    %% re-evaluate the Jacobian at optimal values
    jac = feval(dFdp,x,fbest,p,dp,F,bounds);
    msk = dp ~= 0;
    n = sum(msk);           % reduce n to equal number of estimated parameters
    jac = jac(:, msk);    % use only fitted parameters
    
    %% following section is Ray Muzic's estimate for covariance and correlation
    %% assuming covariance of data is a diagonal matrix proportional to
    %% diag(1/wt.^2).
    %% cov matrix of data est. from Bard Eq. 7-5-13, and Row 1 Table 5.1
    
    %     tp = wt.^2;
    %     if exist('sparse')  % save memory
    %         Q = sparse (1:m, 1:m, 1 ./ tp);
    %         Qinv = sparse (1:m, 1:m, tp);
    %     else
    %         Q = diag (1 ./ tp);
    %         Qinv = diag (tp);
    %     end
    %
    %     %un-weighted residuals
    %     if (~isreal (r))
    %         error ('residuals are not real');
    %     end
    %     tp = resid.' * Qinv * resid;
    %     covr = (tp / m) * Q;    %covariance of residuals
    %
    %     %% Matlab compatibility and avoiding recomputation make the following
    %     %% logic clumsy.
    %     compute = 1;
    %     if (m <= n)
    %         compute = 0;
    %     else
    %         Qinv = ((m - n) / tp) * Qinv;
    %         %% simplified Eq. 7-5-13, Bard; cov of parm est, inverse; outer
    %         %% parantheses contain inverse of guessed covariance matrix of data    covpinv = jac.' * Qinv * jac;
    %         covpinv = jac.' * Qinv * jac;
    %         if (exist ('rcond','builtin')==5 && rcond(covpinv) <= eps)
    %             compute = 0;
    %             warning('MATLAB:Fit_Really_Bad_Fit','The covariance of the pin is malformed (possible error).\nErrors are calculated using a %3.0f%% confidence bound',100*0.95)
    %         elseif (rank (covpinv) < n)
    %             %% above test is not equivalent to 'rcond' and may unnecessarily
    %             %% reject some matrices
    %             compute = 0;
    %         end
    %     end
    %     if (compute)
    %         covp = inv (covpinv);
    %         d=sqrt(diag(covp));
    %         corp = covp ./ (d * d.');
    %     else
    %         covp = nan(n);
    %         corp = covp;
    %     end
    %
    %     if exist('sparse')
    %         covr=spdiags(covr,0);
    %     else
    %         covr=diag(covr);                 % convert returned values to
    %         % compact storage
    %     end
    %     stdresid = resid .* abs (wt) / sqrt (tp / m); % equivalent to resid ./
    %     % sqrt (covr)
    %     %
    
    % THis is the old stuff...
    try
        Qinv=diag(wt.*wt);
        Q=diag((0*wt+1)./(wt.^2));
        %[nrw ncw]=size(wt);
        %Q=ones(nrw,ncw)./wt; Q=diag(Q.*Q);
        resid=y-f;                                    %un-weighted residuals
        covr=resid'*Qinv*resid*Q/(m-n);                 %covariance of residuals
        Vy=1/(1-n/m)*covr;  % Eq. 7-13-22, Bard         %covariance of the data
        covr=diag(covr);                                %for compact storage
        Z=((m-n)*jac'*Qinv*jac)/(n*resid'*Qinv*resid);
        stdresid=resid./sqrt(diag(Vy));
        
        jtgjinv=pinv(jac'*Qinv*jac);
        covp=jtgjinv*jac'*Qinv*Vy*Qinv*jac*jtgjinv; % Eq. 7-5-13, Bard %cov of parm est
        for k=1:n,
            for j=k:n,
                corp(k,j)=covp(k,j)/sqrt(abs(covp(k,k)*covp(j,j)));
                corp(j,k)=corp(k,j);
            end;
        end;
        std=sqrt(diag(covp));
        compute = 1;
        
    catch
        compute = 0;
    end
    
    if compute
        j=1;
        sig=zeros(size(p));
        for i=1:length(std)
            while dp(j)==0
                j=j+1;
            end
            sig(j)=std(i);
            j=j+1;
        end
        std=sig;
    else
        % 95% Confidence region
        delta=confidence(p(msk),resid,jac,conf_p);
        std=zeros(size(p));
        std(msk)=std(msk)+delta;
    end
    
    if (~(verbose || nargout > 8))
        return;
    end
    
    if (m > n)
        Z = ((m - n) / (n * resid.' * Qinv * resid)) * covpinv;
    else
        Z = NaN * ones (n);
    end
    
    %%% alt. est. of cov. mat. of parm.:(Delforge, Circulation, 82:1494-1504, 1990
    %%disp('Alternate estimate of cov. of param. est.')
    %%acovp=resid'*Qinv*resid/(m-n)*inv(jac'*Qinv*jac);
    
    % %%Calculate R^2, intercept form
    % %%
    % % tp=sum((y-mean(y)).^2);
    % tp = sumsq (y - mean (y));
    % % te=sum((resid).^2);
    % if (tp > 0)
    %     r2 = 1 - sumsq (resid) / tp;
    % %     r2=1-te/tp;
    % else
    %     r2 = NaN;
    % end
    
    %Calculate R^2 (Ref Draper & Smith p.46)
    %
    r=corrcoef([y(:),f(:)]);
    r2=r(1,2).^2;
    ra2=1-(1-r2)*((length(y)-1)/(length(y)-sum(msk)-1));
    
    
    %% if someone has asked for it, let them have it
    %%
    if (verbose)
        feval(dFdp,x,fbest,p,dp,F,bounds);
        disp(' Least Squares Estimates of Parameters')
        disp(p.')
        disp(' Correlation matrix of parameters estimated')
        disp(corp)
        disp(' Covariance matrix of Residuals' )
        disp(covr)
        disp(' Correlation Coefficient R^2')
        disp(r2)
        sprintf(' 95%% conf region: F(0.05)(%.0f,%.0f)>= delta_pvec.''*Z*delta_pvec',n,m-n)
        disp(Z)
        %% runs test according to Bard. p 201.
        n1 = sum (resid > 0);
        n2 = sum (resid < 0);
        nrun=sum(abs(diff(resid > 0)))+1;
        if ((n1 > 10) && (n2 > 10)) % sufficent data for test?
            zed=(nrun-(2*n1*n2/(n1+n2)+1)+0.5)/(2*n1*n2*(2*n1*n2-n1-n2)...
                /((n1+n2)^2*(n1+n2-1)));
            if (zed < 0)
                prob = erfc(-zed/sqrt(2))/2*100;
                disp([num2str(prob),'% chance of fewer than ',num2str(nrun),' runs.']);
            else
                prob = erfc(zed/sqrt(2))/2*100;
                disp([num2str(prob),'% chance of greater than ',num2str(nrun),' runs.']);
            end
        end
    end
    
end

function a = sumsq(b,dim)
    % A = SUMSQ(B,DIM)
    %
    % Returns the sum of the square of each element in the matrix
    % If two arguments are specified, then it sums the squares in
    % a particular dimension.
    
    if nargin==1, a = sum(b.*b); return; end
    if ndims(b) < dim, error('DIM exceeds number of dimensions of B.'); end
    a = sum(b.*b,dim);
end

function [lb, idx] = cpiv (v, m)
    
    %% [lb, idx] = cpiv (v, m)
    %%
    %% v: column vector; m: matrix. length (v) must equal rows (m). m must
    %% be positive definit, which is not be explicitely checked. Finds
    %% column vectors w and l with w == v + m * l, w >= 0, l >= 0, l.' * w
    %% == 0. lb: column vector of components of l for which the
    %% corresponding components of w are zero; idx: logical index of these
    %% components in l. This is called solving the 'complementary pivot
    %% problem' (Cottle, R. W. and Dantzig, G. B., 'Complementary pivot
    %% theory of mathematical programming', Linear Algebra and Appl. 1,
    %% 102--125. References for the current algorithm: Bard, Y.: Nonlinear
    %% Parameter Estimation, p. 147--149, Academic Press, New York and
    %% London 1974; Bard, Y., 'An eclectic approach to nonlinear
    %% programming', Proc. ANU Sem. Optimization, Canberra, Austral. Nat.
    %% Univ.).
    
    n = length (v);
    if (n > size (v, 1))
        error ('first argument is no column vector'); % the most typical mistake
    end
    m = cat (2, m, v);
    id = ones (n, 1);
    nz = -eps; % This is arbitrary; components of w and -l are regarded as
    % non-negative if >= nz.
    nl = 100 * n; % maximum number of loop repeats, after that give up
    ready = false;
    while (~ ready && nl > 0)
        [vm, idm] = min (id .* m(:, end));
        if (vm >= nz)
            ready = true;
        else
            id(idm) = -id(idm);
            m = gjp (m, idm);
            nl = nl - 1;
        end
    end
    if (~ ready)
        error ('not successful');
    end
    idx = id < 0;
    lb = -m(idx, end);
end

function m = gjp (m, k, l)
    
    %% m = gjp (m, k[, l])
    %%
    %% m: matrix; k, l: row- and column-index of pivot, l defaults to k.
    %%
    %% Gauss-Jordon pivot as defined in Bard, Y.: Nonlinear Parameter
    %% Estimation, p. 296, Academic Press, New York and London 1974. In
    %% the pivot column, this seems not quite the same as the usual
    %% Gauss-Jordan(-Clasen) pivot. Bard gives Beaton, A. E., 'The use of
    %% special matrix operators in statistical calculus' Research Bulletin
    %% RB-64-51 (1964), Educational Testing Service, Princeton, New Jersey
    %% as a reference, but this article is not easily accessible. Another
    %% reference, whose definition of gjp differs from Bards by some
    %% signs, is Clarke, R. B., 'Algorithm AS 178: The Gauss-Jordan sweep
    %% operator with detection of collinearity', Journal of the Royal
    %% Statistical Society, Series C (Applied Statistics) (1982), 31(2),
    %% 166--168.
    
    if (nargin < 3)
        l = k;
    end
    
    p = m(k, l);
    
    if (p == 0)
        error ('pivot is zero');
    end
    
    %% This is a case where I really hate to remain Matlab compatible,
    %% giving so many indices twice.
    m(k, [1:l-1, l+1:end]) = m(k, [1:l-1, l+1:end]) / p; % pivot row
    m([1:k-1, k+1:end], [1:l-1, l+1:end]) = ... % except pivot row and col
        m([1:k-1, k+1:end], [1:l-1, l+1:end]) - ...
        m([1:k-1, k+1:end], l) * m(k, [1:l-1, l+1:end]);
    m([1:k-1, k+1:end], l) = - m([1:k-1, k+1:end], l) / p; % pivot column
    m(k, l) = 1 / p;
end

function delta = confidence(beta,resid,J,alpha)
    %   Confidence intervals for parameters in nonlinear regression.
    %   delta = confidence(pin,residules,jacobian,confidence)
    %   Thanks to matlab statistics toolbox!
    
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
    [temp,R] = qr(J,0);
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


%!demo
%! % Define functions
%! leasqrfunc = @(x, p) p(1) * exp (-p(2) * x);
%! leasqrdfdp = @(x, f, p, dp, func) [exp(-p(2)*x), -p(1)*x.*exp(-p(2)*x)];
%!
%! % generate test data
%! t = [1:10:100]';
%! p = [1; 0.1];
%! data = leasqrfunc (t, p);
%!
%! rnd = [0.352509; -0.040607; -1.867061; -1.561283; 1.473191; ...
%!        0.580767;  0.841805;  1.632203; -0.179254; 0.345208];
%!
%! % add noise
%! % wt1 = 1 /sqrt of variances of data
%! % 1 / wt1 = sqrt of var = standard deviation
%! wt1 = (1 + 0 * t) ./ sqrt (data);
%! data = data + 0.05 * rnd ./ wt1;
%!
%! % Note by Thomas Walter <walter@pctc.chemie.uni-erlangen.de>:
%! %
%! % Using a step size of 1 to calculate the derivative is WRONG !!!!
%! % See numerical mathbooks why.
%! % A derivative calculated from central differences need: s
%! %     step = 0.001...1.0e-8
%! % And onesided derivative needs:
%! %     step = 1.0e-5...1.0e-8 and may be still wrong
%!
%! F = leasqrfunc;
%! dFdp = leasqrdfdp; % exact derivative
%! % dFdp = dfdp;     % estimated derivative
%! dp = [0.001; 0.001];
%! pin = [.8; .05];
%! stol=0.001; niter=50;
%! minstep = [0.01; 0.01];
%! maxstep = [0.8; 0.8];
%! options = [minstep, maxstep];
%!
%! global verbose;
%! verbose = 1;
%! [f1, p1, kvg1, iter1, corp1, covp1, covr1, stdresid1, Z1, r21] = ...
%!    leasqr (t, data, pin, F, stol, niter, wt1, dp, dFdp, options);

%!demo
%!  %% Example for linear inequality constraints.
%!  %% model function:
%!  F = @ (x, p) p(1) * exp (p(2) * x);
%!  %% independents and dependents:
%!  x = 1:5;
%!  y = [1, 2, 4, 7, 14];
%!  %% initial values:
%!  init = [.25; .25];
%!  %% other configuration (default values):
%!  tolerance = .0001;
%!  max_iterations = 20;
%!  weights = ones (5, 1);
%!  dp = [.001; .001]; % bidirectional numeric gradient stepsize
%!  dFdp = 'dfdp'; % function for gradient (numerical)
%!
%!  %% linear constraints, A.' * parametervector + B >= 0
%!  A = [1; -1]; B = 0; % p(1) >= p(2);
%!  options.inequc = {A, B};
%!
%!  %% start leasqr, be sure that 'verbose' is not set
%!  global verbose; verbose = false;
%!  [f, p, cvg, iter] = ...
%!      leasqr (x, y, init, F, tolerance, max_iterations, ...
%!          weights, dp, dFdp, options)