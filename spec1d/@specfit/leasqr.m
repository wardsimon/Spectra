function [f,p,cvg,iter,corp,covp,covr,stdresid,Z,r2] = leasqr(f_in,s_in)
%% The following two blocks of comments are chiefly from the original
%% version for Matlab. For later changes the logs of the Octave Forge
%% svn repository should also be consulted.

%% A modified version of Levenberg-Marquardt
%% Non-Linear Regression program previously submitted by R.Schrager.
%% This version corrects an error in that version and also provides
%% an easier to use version with automatic numerical calculation of
%% the Jacobian Matrix. In addition, this version calculates statistics
%% such as correlation, etc....
%%
%% Version 3 Notes
%% Errors in the original version submitted by Shrager (now called
%% version 1) and the improved version of Jutan (now called version 2)
%% have been corrected.
%% Additional features, statistical tests, and documentation have also been
%% included along with an example of usage.  BEWARE: Some the the input and
%% output arguments were changed from the previous version.
%%
%%     Ray Muzic     <rfm2@ds2.uh.cwru.edu>
%%     Arthur Jutan  <jutan@charon.engga.uwo.ca>

%% Richard I. Shrager (301)-496-1122
%% Modified by A.Jutan (519)-679-2111
%% Modified by Ray Muzic 14-Jul-1992
%%       1) add maxstep feature for limiting changes in parameter estimates
%%          at each step.
%%       2) remove forced columnization of x (x=x(:)) at beginning. x
%%          could be a matrix with the ith row of containing values of
%%          the independent variables at the ith observation.
%%       3) add verbose option
%%       4) add optional return arguments covp, stdresid, chi2
%%       5) revise estimates of corp, stdev
%% Modified by Ray Muzic 11-Oct-1992
%%	1) revise estimate of Vy.  remove chi2, add Z as return values
%%       (later remark: the current code contains no variable Vy)
%% Modified by Ray Muzic 7-Jan-1994
%%       1) Replace ones(x) with a construct that is compatible with versions
%%          newer and older than v 4.1.
%%       2) Added global declaration of verbose (needed for newer than v4.x)
%%       3) Replace return value var, the variance of the residuals
%%          with covr, the covariance matrix of the residuals.
%%       4) Introduce options as 10th input argument.  Include
%%          convergence criteria and maxstep in it.
%%       5) Correct calculation of xtx which affects coveraince estimate.
%%       6) Eliminate stdev (estimate of standard deviation of
%%          parameter estimates) from the return values.  The covp is a
%%          much more meaningful expression of precision because it
%%          specifies a confidence region in contrast to a confidence
%%          interval.. If needed, however, stdev may be calculated as
%%          stdev=sqrt(diag(covp)).
%%       7) Change the order of the return values to a more logical order.
%%       8) Change to more efficent algorithm of Bard for selecting epsL.
%%       9) Tighten up memory usage by making use of sparse matrices (if
%%          MATLAB version >= 4.0) in computation of covp, corp, stdresid.
%% Modified by Francesco Potortì
%%       for use in Octave
%% Added linear inequality constraints with quadratic programming to
%% this file and special case bounds to this file and to dfdp.m
%% (24-Feb-2010) and later also general inequality constraints
%% (12-Apr-2010) (Reference: Bard, Y., 'An eclectic approach to
%% nonlinear programming', Proc. ANU Sem. Optimization, Canberra,
%% Austral. Nat. Univ.). Differences from the reference: adaption to
%% svd-based algorithm, linesearch or stepwidth adaptions to ensure
%% decrease in objective function omitted to rather start a new
%% overall cycle with a new epsL, some performance gains from linear
%% constraints even if general constraints are specified. Equality
%% constraints also implemented. Olaf Till
%% Now split into files leasqr.m and __lm_svd__.m.
%% Copyright (C) 1992-1994 Richard Shrager
%% Copyright (C) 1992-1994 Arthur Jutan <jutan@charon.engga.uwo.ca>
%% Copyright (C) 1992-1994 Ray Muzic <rfm2@ds2.uh.cwru.edu>
%% Copyright (C) 2010-2016 Olaf Till <i7tiol@t-online.de>
%%
%% This program is free software; you can redistribute it and/or modify it under
%% the terms of the GNU General Public License as published by the Free Software
%% Foundation; either version 3 of the License, or (at your option) any later
%% version.
%%
%% This program is distributed in the hope that it will be useful, but WITHOUT
%% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
%% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
%% details.
%%
%% You should have received a copy of the GNU General Public License along with
%% this program; if not, see <http://www.gnu.org/licenses/>.

%% -*- texinfo -*-
%% @deftypefn  {Function File} {} leasqr (@var{x}, @var{y}, @var{pin}, @var{F})
%% @deftypefnx {Function File} {} leasqr (@var{x}, @var{y}, @var{pin}, @var{F}, @var{stol})
%% @deftypefnx {Function File} {} leasqr (@var{x}, @var{y}, @var{pin}, @var{F}, @var{stol}, @var{niter})
%% @deftypefnx {Function File} {} leasqr (@var{x}, @var{y}, @var{pin}, @var{F}, @var{stol}, @var{niter}, @var{wt})
%% @deftypefnx {Function File} {} leasqr (@var{x}, @var{y}, @var{pin}, @var{F}, @var{stol}, @var{niter}, @var{wt}, @var{dp})
%% @deftypefnx {Function File} {} leasqr (@var{x}, @var{y}, @var{pin}, @var{F}, @var{stol}, @var{niter}, @var{wt}, @var{dp}, @var{dFdp})
%% @deftypefnx {Function File} {} leasqr (@var{x}, @var{y}, @var{pin}, @var{F}, @var{stol}, @var{niter}, @var{wt}, @var{dp}, @var{dFdp}, @var{options})
%% @deftypefnx {Function File} {[@var{f}, @var{p}, @var{cvg}, @var{iter}, @var{corp}, @var{covp}, @var{covr}, @var{stdresid}, @var{Z}, @var{r2}] =} leasqr (@dots{})
%% Levenberg-Marquardt nonlinear regression.
%%
%% Input arguments:
%%
%% @table @var
%% @item x
%% Vector or matrix of independent variables.
%%
%% @item y
%% Vector or matrix of observed values.
%%
%% @item pin
%% Vector of initial parameters to be adjusted by leasqr.
%%
%% @item F
%% Name of function or function handle. The function must be of the form
%% @code{y = f(x, p)}, with y, x, p of the form @var{y}, @var{x}, @var{pin}.
%%
%% @item stol
%% Scalar tolerance on fractional improvement in scalar sum of squares, i.e.,
%% @code{sum ((@var{wt} .* (@var{y}-@var{f}))^2)}.  Set to 0.0001 if
%% empty or not given;
%%
%% @item niter
%% Maximum number of iterations.  Set to 20 if empty or not given.
%%
%% @item wt
%% Statistical weights (same dimensions as @var{y}).  These should be
%% set to be proportional to @code{sqrt (@var{y}) ^-1}, i.e., the
%% covariance matrix of the data is assumed to be proportional to
%% diagonal with diagonal equal to @code{(@var{wt}.^2)^-1}.  The constant of
%% proportionality will be estimated.  Set to @code{ones (size
%% (@var{y}))} if empty or not given.
%%
%% @item dp
%% Fractional increment of @var{p} for numerical partial derivatives.  Set
%% to @code{0.001 * ones (size (@var{pin}))} if empty or not given.
%%
%% @itemize @bullet
%% @item dp(j) > 0 means central differences on j-th parameter p(j).
%% @item dp(j) < 0 means one-sided differences on j-th parameter p(j).
%% @item dp(j) = 0 holds p(j) fixed, i.e., leasqr won't change initial guess: pin(j)
%% @end itemize
%%
%% @item dFdp
%% Name of partial derivative function in quotes or function handle. If
%% not given or empty, set to @code{dfdp}, a slow but general partial
%% derivatives function. The function must be of the form @code{prt =
%% dfdp (x, f, p, dp, F [,bounds])}.  For backwards compatibility, the
%% function will only be called with an extra 'bounds' argument if the
%% 'bounds' option is explicitly specified to leasqr (see dfdp.m).
%%
%% @item options
%% Structure with multiple options. The following fields are recognized:
%%
%% @table @asis
%% @item @qcode{fract_prec}
%% Column vector (same length as @var{pin})
%% of desired fractional precisions in parameter estimates.
%% Iterations are terminated if change in parameter vector (chg)
%% relative to current parameter estimate is less than their
%% corresponding elements in 'fract_prec', i.e.,
%% @code{all (abs (chg) < abs (options.fract_prec .* current_parm_est))} on two
%% consecutive iterations. Defaults to @code{zeros (size (@var{pin}))}.
%%
%% @item @qcode{max_fract_change}
%% Column vector (same length as @var{pin}) of maximum fractional step
%% changes in parameter vector.
%% Fractional change in elements of parameter vector is constrained to
%% be at most 'max_fract_change' between sucessive iterations, i.e.,
%% @code{abs (chg(i)) = abs (min([chg(i), options.max_fract_change(i) * current param estimate]))}.
%% Defaults to @code{Inf * ones (size (@var{pin}))}.
%%
%% @item @qcode{inequc}
%% Cell-array containing up to four entries,
%% two entries for linear inequality constraints and/or one or two
%% entries for general inequality constraints.  Initial parameters
%% must satisfy these constraints.  Either linear or general
%% constraints may be the first entries, but the two entries for
%% linear constraints must be adjacent and, if two entries are given
%% for general constraints, they also must be adjacent.  The two
%% entries for linear constraints are a matrix (say m) and a vector
%% (say v), specifying linear inequality constraints of the form
%% `m.' * parameters + v >= 0'. If the constraints are just bounds,
%% it is suggested to specify them in 'options.bounds' instead,
%% since then some sanity tests are performed, and since the
%% function 'dfdp.m' is guarantied not to violate constraints during
%% determination of the numeric gradient only for those constraints
%% specified as 'bounds' (possibly with violations due to a certain
%% inaccuracy, however, except if no constraints except bounds are
%% specified). The first entry for general constraints must be a
%% differentiable vector valued function (say h), specifying general
%% inequality constraints of the form `h (p[, idx]) >= 0'; p is the
%% column vector of optimized paraters and the optional argument idx
%% is a logical index. h has to return the values of all constraints
%% if idx is not given, and has to return only the indexed
%% constraints if idx is given (so computation of the other
%% constraints can be spared). If a second entry for general
%% constraints is given, it must be a function (say dh) which
%% returnes a matrix whos rows contain the gradients of the
%% constraint function h with respect to the optimized parameters.
%% It has the form jac_h = dh (vh, p, dp, h, idx[, bounds]); p is
%% the column vector of optimized parameters, and idx is a logical
%% index --- only the rows indexed by idx must be returned (so
%% computation of the others can be spared). The other arguments of
%% dh are for the case that dh computes numerical gradients: vh is
%% the column vector of the current values of the constraint
%% function h, with idx already applied. h is a function h (p) to
%% compute the values of the constraints for parameters p, it will
%% return only the values indexed by idx. dp is a suggestion for
%% relative step width, having the same value as the argument 'dp'
%% of leasqr above. If bounds were specified to leasqr, they are
%% provided in the argument bounds of dh, to enable their
%% consideration in determination of numerical gradients. If dh is
%% not specified to leasqr, numerical gradients are computed in the
%% same way as with 'dfdp.m' (see above). If some constraints are
%% linear, they should be specified as linear constraints (or
%% bounds, if applicable) for reasons of performance, even if
%% general constraints are also specified.
%%
%% @item @qcode{bounds}
%% Two-column-matrix, one row for each
%% parameter in @var{pin}. Each row contains a minimal and maximal value
%% for each parameter. Default: [-Inf, Inf] in each row. If this
%% field is used with an existing user-side function for 'dFdp'
%% (see above) the functions interface might have to be changed.
%%
%% @item @qcode{equc}
%% Equality constraints, specified the same
%% way as inequality constraints (see field 'options.inequc').
%% Initial parameters must satisfy these constraints.
%% Note that there is possibly a certain inaccuracy in honoring
%% constraints, except if only bounds are specified.
%% @emph{Warning}: If constraints (or bounds) are set, returned guesses
%% of @var{corp}, @var{covp}, and @var{Z} are generally invalid, even if
%% no constraints
%% are active for the final parameters. If equality constraints are
%% specified, @var{corp}, @var{covp}, and @var{Z} are not guessed at all.
%%
%% @item @qcode{cpiv}
%% Function for complementary pivot algorithm
%% for inequality constraints. Defaults to cpiv_bard.  No different
%% function is supplied.
%%
%% @end table
%%
%% For backwards compatibility, @var{options} can also be a matrix whose
%% first and second column contains the values of @qcode{fract_prec} and
%% @qcode{max_fract_change}, respectively.
%%
%% @end table
%%
%% Output:
%%
%% @table @var
%% @item f
%% Column vector of values computed: f = F(x,p).
%%
%% @item p
%% Column vector trial or final parameters, i.e, the solution.
%%
%% @item cvg
%% Scalar: = 1 if convergence, = 0 otherwise.
%%
%% @item iter
%% Scalar number of iterations used.
%%
%% @item corp
%% Correlation matrix for parameters.
%%
%% @item covp
%% Covariance matrix of the parameters.
%%
%% @item covr
%% Diag(covariance matrix of the residuals).
%%
%% @item stdresid
%% Standardized residuals.
%%
%% @item Z
%% Matrix that defines confidence region (see comments in the source).
%%
%% @item r2
%% Coefficient of multiple determination, intercept form.
%%
%% @end table
%%
%% Not suitable for non-real residuals.
%%
%% References:
%% Bard, Nonlinear Parameter Estimation, Academic Press, 1974.
%% Draper and Smith, Applied Regression Analysis, John Wiley and Sons, 1981.
%%
%% @end deftypefn
%   __plot_cmds__ (); % flag persistent variables invalid

verbose = f_in.fitoptions.Display;

%% argument processing
%%

%   if (nargin < 4)
%     print_usage ();
%   end

%   if (nargin > 8 && ~ isempty (dFdp))
%     if (ischar (dFdp))
%       dfdp = str2func (dFdp);
%     else
%       dfdp = dFdp;
%     end
%   end
%
%   if (nargin <= 7 || isempty (dp)) dp=.001*(pin*0+1); end %DT
%   if (nargin <= 6 || isempty (wt)) wt = ones (size (y)); end %SMB modification
%   if (nargin <= 5) niter = []; end
%   if (nargin == 4 || isempty (stol)) stol=.0001; end
%   if (ischar (F)) F = str2func (F); end
%%

options = f_in.fitoptions;

x = s_in.x(:);
y = s_in.y(:);
wt = 1./s_in.e(:);

pin = f_in.pvals;
dpin = f_in.notfixed;
fcp(1) = options.TolFun;

if (any (size (y) ~= size (wt)))
    error ('dimensions of observations and weights do not match');
end

wtl = wt(:);
pin = pin(:);
dp = -(dpin>0)*fcp(1); %change all vectors to columns
[rows_y, cols_y] = size (y);
m = rows_y * cols_y; n=length(pin);
f_pin = feval(f_in.func,x, pin);
if (any (size (f_pin) ~= size (y)))
    error ('dimensions of returned values of model function and of observations do not match');
end
f_pin = y - f_pin;

dFdp = @ (p, dfdp_hook) - dfdp (x, y(:) - dfdp_hook.f, p, dp, F);

%% processing of 'options'
pprec = zeros (n, 1);
maxstep = Inf * ones (n, 1);
have_gencstr = false; % no general constraints
have_genecstr = false; % no general equality constraints
n_gencstr = 0;
mc = zeros (n, 0);
vc = zeros (0, 1); rv = 0;
emc = zeros (n, 0);
evc = zeros (0, 1); erv = 0;
bounds = cat (2, -Inf * ones (n, 1), Inf * ones (n, 1));
pin_cstr.inequ.lin_except_bounds = [];
pin_cstr.inequ.gen = [];
pin_cstr.equ.lin = [];
pin_cstr.equ.gen = [];
dfdp_bounds = {};
cpiv = @ cpiv_bard;
eq_idx = []; % numerical index for equality constraints in all
% constraints, later converted to
% logical index
if (isnumeric (options)) % backwards compatibility
    tp = options;
    options = struct ('fract_prec', tp(:, 1));
    if (columns (tp) > 1)
        options.max_fract_change = tp(:, 2);
    end
end
if (isfield (options, 'cpiv') && ~isempty (options.cpiv))
    %% As yet there is only one cpiv function distributed with leasqr,
    %% but this may change; the algorithm of cpiv_bard is said to be
    %% relatively fast, but may have disadvantages.
    if (ischar (options.cpiv))
        cpiv = str2func (options.cpiv);
    else
        cpiv = options.cpiv;
    end
end

if (isfield (options, 'fract_prec'))
    pprec = options.fract_prec;
    if (any (size (pprec) ~= [n, 1]))
        error ('fractional precisions: wrong dimensions');
    end
end

if (isfield (options, 'max_fract_change'))
    maxstep = options.max_fract_change;
    if (any (size (maxstep) ~= [n, 1]))
        error ('maximum fractional step changes: wrong dimensions');
    end
end

if (isfield (options, 'inequc'))
    inequc = options.inequc;
    if (isnumeric (inequc{1}))
        mc = inequc{1};
        vc = inequc{2};
        if (length (inequc) > 2)
            have_gencstr = true;
            f_gencstr = inequc{3};
            if (length (inequc) > 3)
                df_gencstr = inequc{4};
            else
                df_gencstr = @ dcdp;
            end
        end
    else
        lid = 0; % no linear constraints
        have_gencstr = true;
        f_gencstr = inequc{1};
        if (length (inequc) > 1)
            if (isnumeric (inequc{2}))
                lid = 2;
                df_gencstr = @ dcdp;
            else
                df_gencstr = inequc{2};
                if (length (inequc) > 2)
                    lid = 3;
                end
            end
        else
            df_gencstr = @ dcdp;
        end
        if (lid)
            mc = inequc{lid};
            vc = inequc{lid + 1};
        end
    end
    
    if (have_gencstr)
        if (ischar (f_gencstr))
            f_gencstr = str2func (f_gencstr);
        end
        tp = f_gencstr (pin);
        n_gencstr = length (tp);
        f_gencstr = @ (p, idx) tf_gencstr (p, idx, f_gencstr);
        if (ischar (df_gencstr))
            df_gencstr = str2func (df_gencstr);
        end
        if (strcmp (func2str (df_gencstr), 'dcdp'))
            df_gencstr = @ (f, p, dp, idx, db) ...
                df_gencstr (f(idx), p, dp, ...
                @ (tp) f_gencstr (tp, idx), db{:});
        else
            df_gencstr = @ (f, p, dp, idx, db) ...
                df_gencstr (f(idx), p, dp, ...
                @ (tp) f_gencstr (tp, idx), idx, db{:});
        end
    end
    [rm, cm] = size (mc);
    [rv, cv] = size (vc);
    if (rm ~= n || cm ~= rv || cv ~= 1)
        error ('linear inequality constraints: wrong dimensions');
    end
    pin_cstr.inequ.lin_except_bounds = mc.' * pin + vc;
    if (have_gencstr)
        pin_cstr.inequ.gen = tp;
    end
end

if (isfield (options, 'equc'))
    equc = options.equc;
    if (isnumeric (equc{1}))
        emc = equc{1};
        evc = equc{2};
        if (length (equc) > 2)
            have_genecstr = true;
            f_genecstr = equc{3};
            if (length (equc) > 3)
                df_genecstr = equc{4};
            else
                df_genecstr = @ dcdp;
            end
        end
    else
        lid = 0; % no linear constraints
        have_genecstr = true;
        f_genecstr = equc{1};
        if (length (equc) > 1)
            if (isnumeric (equc{2}))
                lid = 2;
                df_genecstr = @ dcdp;
            else
                df_genecstr = equc{2};
                if (length (equc) > 2)
                    lid = 3;
                end
            end
        else
            df_genecstr = @ dcdp;
        end
        if (lid)
            emc = equc{lid};
            evc = equc{lid + 1};
        end
    end
    if (have_genecstr)
        if (ischar (f_genecstr))
            f_genecstr = str2func (f_genecstr);
        end
        tp = f_genecstr (pin);
        n_genecstr = length (tp);
        f_genecstr = @ (p, idx) tf_gencstr (p, idx, f_genecstr);
        if (ischar (df_genecstr))
            df_genecstr = str2func (df_genecstr);
        end
        if (strcmp (func2str (df_genecstr), 'dcdp'))
            df_genecstr = @ (f, p, dp, idx, db) ...
                df_genecstr (f, p, dp, ...
                @ (tp) f_genecstr (tp, idx), db{:});
        else
            df_genecstr = @ (f, p, dp, idx, db) ...
                df_genecstr (f, p, dp, ...
                @ (tp) f_genecstr (tp, idx), idx, db{:});
        end
    end
    [erm, ecm] = size (emc);
    [erv, ecv] = size (evc);
    if (erm ~= n || ecm ~= erv || ecv ~= 1)
        error ('linear equality constraints: wrong dimensions');
    end
    pin_cstr.equ.lin = emc.' * pin + evc;
    if (have_genecstr)
        pin_cstr.equ.gen = tp;
    end
end

if (isfield (options, 'bounds'))
    bounds = options.bounds;
    if (any (size (bounds) ~= [n, 2]))
        error ('bounds: wrong dimensions');
    end
    idx = bounds(:, 1) > bounds(:, 2);
    tp = bounds(idx, 2);
    bounds(idx, 2) = bounds(idx, 1);
    bounds(idx, 1) = tp;
    %% It is possible to take this decision here, since this frontend
    %% is used only with one certain backend. The backend will check
    %% this again; but it will not reference 'dp' in its message,
    %% thats why the additional check here.
    idx = bounds(:, 1) == bounds(:, 2);
    if (any (idx))
        warning ('leasqr:constraints', 'lower and upper bounds identical for some parameters, setting the respective elements of dp to zero');
        dp(idx) = 0;
    end
    %%
    tp = eye (n);
    lidx = ~isinf (bounds(:, 1));
    uidx = ~isinf (bounds(:, 2));
    mc = cat (2, mc, tp(:, lidx), - tp(:, uidx));
    vc = cat (1, vc, - bounds(lidx, 1), bounds(uidx, 2));
    [rm, cm] = size (mc);
    [rv, cv] = size (vc);
    dfdp_bounds = {bounds};
    dFdp = @ (p, dfdp_hook) - dfdp (x, y(:) - dfdp_hook.f, p, dp, ...
        F, bounds);
end

%% concatenate inequality and equality constraint functions, mc, and
%% vc; update eq_idx, rv, n_gencstr, have_gencstr
if (erv > 0)
    mc = cat (2, mc, emc);
    vc = cat (1, vc, evc);
    eq_idx = rv + 1 : rv + erv;
    rv = rv + erv;
end
if (have_genecstr)
    eq_idx = cat (2, eq_idx, ...
        rv + n_gencstr + 1 : rv + n_gencstr + n_genecstr);
    nidxi = 1 : n_gencstr;
    nidxe = n_gencstr + 1 : n_gencstr + n_genecstr;
    n_gencstr = n_gencstr + n_genecstr;
    if (have_gencstr)
        f_gencstr = @ (p, idx) cat (1, ...
            f_gencstr (p, idx(nidxi)), ...
            f_genecstr (p, idx(nidxe)));
        df_gencstr = @ (f, p, dp, idx, db) ...
            cat (1, ...
            df_gencstr (f(nidxi), p, dp, idx(nidxi), db), ...
            df_genecstr (f(nidxe), p, dp, idx(nidxe), db));
    else
        f_gencstr = f_genecstr;
        df_gencstr = df_genecstr;
        have_gencstr = true;
    end
end

if (have_gencstr)
    nidxl = 1:rv;
    nidxh = rv+1:rv+n_gencstr;
    f_cstr = @ (p, idx) ...
        cat (1, mc(:, idx(nidxl)).' * p + vc(idx(nidxl), 1), ...
        f_gencstr (p, idx(nidxh)));
    %% in the case of this interface, diffp is already zero at fixed;
    %% also in this special case, dfdp_bounds can be filled in directly
    %% --- otherwise it would be a field of hook in the called function
    df_cstr = @ (p, idx, dfdp_hook) ...
        cat (1, mc(:, idx(nidxl)).', ...
        df_gencstr (dfdp_hook.f(nidxh), p, dp, ...
        idx(nidxh), ...
        dfdp_bounds));
else
    f_cstr = @ (p, idx) mc(:, idx).' * p + vc(idx, 1);
    df_cstr = @ (p, idx, dfdp_hook) mc(:, idx).';
end



%% in a general interface, check for all(fixed) here

%% passed constraints
hook.mc = mc; % matrix of linear constraints
hook.vc = vc; % vector of linear constraints
hook.f_cstr = f_cstr; % function of all constraints
hook.df_cstr = df_cstr; % function of derivatives of all constraints
hook.n_gencstr = n_gencstr; % number of non-linear constraints
hook.eq_idx = false (size (vc, 1) + n_gencstr, 1);
hook.eq_idx(eq_idx) = true; % logical index of equality constraints in
% all constraints
hook.lbound = bounds(:, 1); % bounds, subset of linear inequality
% constraints in mc and vc
hook.ubound = bounds(:, 2);

%% passed values of constraints for initial parameters
hook.pin_cstr = pin_cstr;

%% passed derivative of model function
hook.dfdp = dFdp;

%% passed function for complementary pivoting
hook.cpiv = cpiv;

%% passed value of residual function for initial parameters
hook.f_pin = f_pin;

%% passed options
hook.max_fract_change = maxstep;
hook.fract_prec = pprec;
hook.TolFun = stol;
hook.MaxIter = niter;
hook.weights = wt;
hook.fixed = dp == 0;
if (verbose)
    hook.Display = 'iter';
    hook.plot_cmd = @ (f) 0; % `plot_cmd' is deprecated
    %     hook.user_interaction = ...
    %         {@ (p, v, s) ...
    %          {ifelse(strcmp(s, 'iter'),
    %                  false(__plot_cmds__(x, y, y - v.residual)),
    %                  false),
    %           []}{:}};
else
    hook.Display = 'off';
    hook.plot_cmd = @ (f) 0; % `plot_cmd' is deprecated
    hook.user_interaction = {};
end

%% only preliminary, for testing
hook.testing = false;
hook.new_s = false;

[p, resid, cvg, outp] = lm_svd (@ (p) y - F (x, p), pin, hook);
f = y - resid;
iter = outp.niter;
cvg = cvg > 0;

if (~cvg) disp(' CONVERGENCE NOT ACHIEVED~ '); end

if (~(verbose || nargout > 4)) return; end

yl = y(:);
f = f(:);

%% CALC VARIANCE COV MATRIX AND CORRELATION MATRIX OF PARAMETERS
%% re-evaluate the Jacobian at optimal values
jac = dFdp (p, struct ('f', f));
msk = ~hook.fixed;
n = sum(msk);           % reduce n to equal number of estimated parameters
jac = jac(:, msk);	% use only fitted parameters

%% following section is Ray Muzic's estimate for covariance and correlation
%% assuming covariance of data is a diagonal matrix proportional to
%% diag(1/wt.^2).
%% cov matrix of data est. from Bard Eq. 7-5-13, and Row 1 Table 5.1

tp = wtl.^2;
if (exist('sparse'))  % save memory
    Q = sparse (1:m, 1:m, 1 ./ tp);
    Qinv = sparse (1:m, 1:m, tp);
else
    Q = diag (ones (m, 1) ./ tp);
    Qinv = diag (tp);
end
resid = resid(:); % un-weighted residuals
if (~isreal (resid)); error ('residuals are not real'); end
tp = resid.' * Qinv * resid;
covr = (tp / m) * Q;    %covariance of residuals

%% Matlab compatibility and avoiding recomputation make the following
%% logic clumsy.
compute = 1;
if (m <= n || any (eq_idx))
    compute = 0;
else
    Qinv = ((m - n) / tp) * Qinv;
    %% simplified Eq. 7-5-13, Bard; cov of parm est, inverse; outer
    %% parantheses contain inverse of guessed covariance matrix of data
    covpinv = jac.' * Qinv * jac;
    if (exist ('rcond'))
        if (rcond (covpinv) <= eps)
            compute = 0;
        end
    elseif (rank (covpinv) < n)
        %% above test is not equivalent to 'rcond' and may unnecessarily
        %% reject some matrices
        compute = 0;
    end
end
if (compute)
    covp = inv (covpinv);
    d=sqrt(diag(covp));
    corp = covp ./ (d * d.');
else
    covp = NA * ones (n);
    corp = covp;
end

if (exist('sparse'))
    covr=spdiags(covr,0);
else
    covr=diag(covr);                 % convert returned values to
    % compact storage
end
covr = reshape (covr, rows_y, cols_y);
stdresid = resid .* abs (wtl) / sqrt (tp / m); % equivalent to resid ./
% sqrt (covr)
stdresid = reshape (stdresid, rows_y, cols_y);

if (~(verbose || nargout > 8)) return; end

if (m > n && ~any (eq_idx))
    Z = ((m - n) / (n * resid.' * Qinv * resid)) * covpinv;
else
    Z = NA * ones (n);
end

%%% alt. est. of cov. mat. of parm.:(Delforge, Circulation, 82:1494-1504, 1990
%%disp('Alternate estimate of cov. of param. est.')
%%acovp=resid'*Qinv*resid/(m-n)*inv(jac'*Qinv*jac);

if (~(verbose || nargout > 9)) return; end

%%Calculate R^2, intercept form
%%
tp = sumsq (yl - mean (yl));
if (tp > 0)
    r2 = 1 - sumsq (resid) / tp;
else
    r2 = NA;
end

%% if someone has asked for it, let them have it
%%
if (verbose)
    %     __plot_cmds__ (x, y, f);
    disp(' Least Squares Estimates of Parameters')
    disp(p.')
    disp(' Correlation matrix of parameters estimated')
    disp(corp)
    disp(' Covariance matrix of Residuals' )
    disp(covr)
    disp(' Correlation Coefficient R^2')
    disp(r2)
    fprintf(' 95%% conf region: F(0.05)(%.0f,%.0f)>= delta_pvec.%s*Z*delta_pvec\n', n, m - n, char (39)); % works with ' and '
    Z
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

function ret = tf_gencstr (p, idx, f)

%% necessary since user function f_gencstr might return [] or a row
%% vector

ret = f (p, idx);
if (isempty (ret))
    ret = zeros (0, 1);
elseif (size (ret, 2) > 1)
    ret = ret(:);
end

end

function [lb, idx, ridx, m] = cpiv_bard (v, m, incl)
%  [lb, idx, ridx, mv] = cpiv_bard (v, m[, incl])
%
%  v: column vector; m: matrix; incl (optional): index. length (v)
%  must equal rows (m). Finds column vectors w and l with w == v + m *
%  l, w >= 0, l >= 0, l.' * w == 0. Chooses idx, w, and l so that
%  l(~idx) == 0, l(idx) == -inv (m(idx, idx)) * v(idx), w(idx) roughly
%  == 0, and w(~idx) == v(~idx) + m(idx, ~idx).' * l(idx). idx indexes
%  at least everything indexed by incl, but l(incl) may be < 0. lb:
%  l(idx) (column vector); idx: logical index, defined above; ridx:
%  ~idx & w roughly == 0; mv: [m, v] after performing a Gauss-Jordan
%  'sweep' (with gjp.m) on each diagonal element indexed by idx.
%  Except the handling of incl (which enables handling of equality
%  constraints in the calling code), this is called solving the
%  'complementary pivot problem' (Cottle, R. W. and Dantzig, G. B.,
%  'Complementary pivot theory of mathematical programming', Linear
%  Algebra and Appl. 1, 102--125. References for the current
%  algorithm: Bard, Y.: Nonlinear Parameter Estimation, p. 147--149,
%  Academic Press, New York and London 1974; Bard, Y., 'An eclectic
%  approach to nonlinear programming', Proc. ANU Sem. Optimization,
%  Canberra, Austral. Nat. Univ.).

n = length (v);
if (n > size (v, 1))
    error ('first argument is no column vector'); % the most typical mistake
end
if (nargin < 3)
    incl = [];
elseif (islogical (incl))
    incl = find (incl);
end
nincl = 1:n;
nincl(incl) = [];
sgn = ones (n, 1);
if (length (incl) == n)
    sgn = - sgn;
    m = cat (2, inv(m), m\v);
else
    m = cat (2, m, v);
    for id = incl(:).'
        sgn(id) = -sgn(id);
        m = gjp (m, id);
    end
end
nz = eps; % This is arbitrary; components of w and -l are regarded as
% non-negative if >= -nz.
nl = 100 * n; % maximum number of loop repeats, after that give up
if (isempty (nincl))
    ready = true;
else
    ready = false;
    while (~ready && nl > 0)
        [vm, idm] = min (sgn(nincl) .* m(nincl, end));
        if (vm >= -nz)
            ready = true;
        else
            idm = nincl(idm);
            sgn(idm) = -sgn(idm);
            m = gjp (m, idm);
            nl = nl - 1;
        end
    end
end
if (~ready)
    error ('not successful');
end
idx = sgn < 0;
lb = -m(idx, end);
ridx = ~idx & abs (m(:, end)) <= nz;


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

% This is a case where I really hate to remain Matlab compatible,
% giving so many indices twice.
m(k, [1:l-1, l+1:end]) = m(k, [1:l-1, l+1:end]) / p; % pivot row
m([1:k-1, k+1:end], [1:l-1, l+1:end]) = ... % except pivot row and col
    m([1:k-1, k+1:end], [1:l-1, l+1:end]) - ...
    m([1:k-1, k+1:end], l) * m(k, [1:l-1, l+1:end]);
m([1:k-1, k+1:end], l) = - m([1:k-1, k+1:end], l) / p; % pivot column
m(k, l) = 1 / p;

end