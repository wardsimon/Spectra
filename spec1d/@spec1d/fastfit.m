function varargout = fastfit(s1,func,pin,notfixed)
%% [pin_out sig_out] = fastfit(s_in,func,pin,notfixed)
%
% This fitting routine for spec1d objects implements the bare essencials in
% order to get the results in the fastest possible time. This is achieved
% by reducing features, simplifying the algorithms and reducing function
% calls. The following should be noted:
% 
% There is no pre-checking of the data. Inf, NaN and zero errors will
% throw errors
% Features like setting optimisers/algorithms, multifitting and selecting
% windows etc. have been sacraficed. 
% The output is NOT a spec1d or fit object, just the parameters and/or
% errors.
%
% Brief testing has shown that pin_out is the same as the fits
% LM implementation, but the errors are different. Speedup on easy
% functions is ~15.
%
% THIS IMPLEMENTATION IS FULLY PARALLEL. There are no worries on setting
% this off on workers as in the fits function.

xf = zeros(length(pin),length(s1));
J = cell(1,length(s1));
options = flmsolve_fast('default'); % Stops the repeated fn call.

for il = 1:length(s1)
    criteria_func = @(pars) feval(@least_square,s1(il).y,s1(il).e,feval(func,s1(il).x,pars));
    [xf(:,il), ~, ~, flag, J{il}] = feval(@flmsolve_fast,criteria_func,pin,notfixed,options);
end
% These are the fit parameters
varargout{1} = xf;

% These are the errors, calculate if needed. 
if nargout == 2
    ste = zeros(length(pin),length(s1));
    for il = 1:length(s1)
        Jl = J{il};
        wt = 1./s1(il).e;
        m  = length(s1(il).y);
        n  = sum(notfixed);
        
        Q       = diag(1./(wt.^2));
        Qinv    = diag(wt.^2);
        jtgjinv = pinv(Jl'*Qinv*Jl);
        
        resid = s1(il).y - feval(func,s1(il).x,xf(:,il));
        covr  = resid'*Qinv*resid*Q/(m-n);
        Vy    = 1/(1-n/m)*covr;
        
        covp = jtgjinv*Jl'*Qinv*Vy*Qinv*Jl*jtgjinv;
        ste(:,il) = sqrt(diag(covp));
    end
    varargout{2} = ste;
end
