function varargout = samin(varargin)
%%  samin: simulated annealing minimization of a function.
%   [x, obj, convergence] = samin('f', {args}, {control})
%   Arguments:
%   'f'         : function name (string)
%   {args}      : a cell array that holds all arguments of the function,
%   {control}   : a cell array with 11 elements
%                   LB      - vector of lower bounds
%                   UB      - vector of upper bounds
%                   nt      - integer: # of iterations between temperature reductions
%                   ns      - integer: # of iterations between bounds adjustments
%                   rt      - 0 < rt <1: temperature reduction factor
%                   maxevals- integer: limit on function evaluations
%                   neps    - integer:  # number of values final result is compared to
%                   functol -   > 0: the required tolerance level for function value comparisons
%                   paramtol-  > 0: the required tolerance level for parameters
%                   verbosity - scalar: 0, 1, or 2.
%                                       0 = no screen output
%                                       1 = summary every temperature change
%                                       2 = only final results to screen
%                   minarg - integer: which of function args is minimization over?
%
%   Returns:
%   x: the minimizer
%   obj: the value of f() at x
%   convergence: 1 if normal conv, other values if not
%
%% Copyright (C) 2004   Michael Creel   <michael.creel@uab.es>
%%
%%  This program is free software; you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation; either version 2 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program; if not, write to the Free Software
%%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%% simann.cc (c) 2004 Michael Creel <michael.creel@uab.es>
%% References:
%%
%% The code follows this article:
%% Goffe, William L. (1996) "SIMANN: A Global Optimization Algorithm
%%	using Simulated Annealing " Studies in Nonlinear Dynamics & Econometrics
%%  Oct96, Vol. 1 Issue 3.
%%
%% The code uses the same names for control variables,
%% for the most part. A notable difference is that the initial
%% temperature is found automatically to ensure that the active
%% bounds when the temperature begins to reduce cover the entire
%% parameter space (defined as a n-dimensional
%% rectangle that is the Cartesian product of the
%% (lb_i, ub_i), i = 1,2,..n
%%
%% Also of note:
%% Corana et. al., (1987) "Minimizing Multimodal Functions of Continuous
%%	Variables with the "Simulated Annealing" Algorithm",
%% 	ACM Transactions on Mathematical Software, V. 13, N. 3.
%%
%% Goffe, et. al. (1994) "Global Optimization of Statistical Functions
%% 	with Simulated Annealing", Journal of Econometrics,
%% 	V. 60, N. 1/2.

if nargin==0
    msg={...
        'Example: A really stupid way to calculate pi\n'...
        'f(x) cos(x) + 0.01*(x-pi)^2;\n'...
        '\n'...
        'Set up the controls:\n'...
        'ub = 20;\n'...
        'lb = -ub;\n'...
        'nt = 20;\n'...
        'ns = 5;\n'...
        'rt = 0.5;\n'...
        'maxevals = 1e10;\n'...
        'neps = 5;\n'...
        'functol = 1e-10;\n'...
        'paramtol = 1e-5;\n'...
        'verbosity = 2;\n'...
        'minarg = 1;\n'...
        '\n'...
        'Put them in a cell array:\n'...
        'control = {lb, ub, nt, ns, rt, maxevals,\n'...
        'neps, functol, paramtol, verbosity, minarg};\n'...
        '\n'...
        'Call the minimizer (function argument also in cell array):\n'...
        'samin(''f'', {-8}, control)\n'...
        '\n'...
        'The result:\n'...
        '================================================\n'...
        'SAMIN final results\n'...
        'Sucessful convergence to tolerance 0.000010\n'...
        '\n'...
        'Obj. fn. value -1.000000\n'...
        '       parameter        search width\n'...
        '        3.141597            0.002170\n'...
        '================================================\n'...
        'ans = 3.1416\n'};
    fprintf([msg{:}])
elseif nargin ~=3
    error('You must supply 3 arguments')
    return
else
    pn=1;
    obj_fn=varargin{1};
    f_args=varargin{2};
    control=varargin{3};
    
    if strcmp(obj_fn(1),'@')
        obj_fn=eval(obj_fn);
        pn=0;
    end
    
    % 	int m, i, j, h, n, nacc, func_evals;
    % 	int nup, nrej, nnew, ndown, lnobds;
    % 	int converge, test;
    
    lb=control{1};
    ub=control{2};
    nt=control{3};
    ns=control{4};
    rt=control{5};
    maxevals=control{6};
    neps=control{7};
    functol=control{8};
    paramtol=control{9};
    verbosity=control{10};
    minarg=control{11};
    
    
    % 	double f, fp, p, pp, fopt, rand_draw, ratio, t;
    
    % 	%% type checking for minimization parameter done here, since we don't know minarg
    % 	%% until now
    % 	if (!(f_args(minarg - 1).is_real_matrix() || (f_args(minarg - 1).is_real_scalar())))
    % 	{
    % 		error("samin: minimization must be with respect to a column vector");
    % 		return octave_value_list();
    % 	}
    % 	if ((f_args(minarg - 1).is_real_matrix()) && (f_args(minarg - 1).columns() != 1))
    % 	{
    %         	error("samin: minimization must be with respect to a column vector");
    %         	return octave_value_list();
    % 	}
    
    x  = f_args(minarg);
    bounds = ub(:) - lb(:);
    n = length(x);
    xopt = x;
    xp=zeros(n,1);    % Might be needed
    
    %%  Set initial values
    nacc = 0;
    func_evals = 0;
    
    fstar=1e20*ones(neps,1);
    
    %% check for out-of-bounds starting value
    for i=(1:n)
        if(( x(i) > ub(i)))
            warning('Starting point out of bounds, setting to upper bound')
            x(i)=ub(i);
        elseif(x(i) < lb(i))
            warning('Starting point out of bounds, setting to lower bound')
            x(i)=lb(i);
        end
    end
    
    %% Initial obj_value
    c_args{1} = obj_fn;
    c_args{2} = f_args;
    f_return  = feval(obj_fn,f_args);
    f         = f_return;
    
    fopt = f;
    fstar(1) = f;
    
    %% First stage: find initial temperature so that
    %% bounds grow to cover parameter space
    t = 1000;
    converge = 0;
    while (converge==0)
        nup = 0;
        nrej = 0;
        nnew = 0;
        ndown = 0;
        lnobds = 0;
        
        %% repeat nt times then adjust temperature
        for m = (1:nt)
            nacp=zeros(n,1);  % Might be needed
            %% repeat ns times, then adjust bounds
            for j= (1:ns)
                
                %% generate new point by taking last
                %% and adding a random value to each of elements,
                %% in turn
                for h = (1:n)
                    xp = x;
                    rand_draw = rand;
                    xp(h) = x(h) + (2.0 * rand_draw - 1.0) * bounds(h);
                    if (xp(h) < lb(h)) || (xp(h) > ub(h))
                        xp(h) = lb(h) + (ub(h) - lb(h))*rand_draw;
                        lnobds = lnobds + 1;
                    end
                    
                    %% Evaluate function at new point
                    f_args(minarg) = xp;
                    c_args{2} = f_args;
                    f_return = feval(obj_fn,f_args);
                    fp = f_return;
                    func_evals = func_evals + 1;
                    
                    %%  If too many function evaluations occur, terminate the algorithm.
                    if(func_evals >= maxevals)
                        warning('NO CONVERGENCE: MAXEVALS exceeded before initial temparature found');
                        if(verbosity >= 1)
                            fprintf('\n================================================\n');
                            fprintf('SAMIN results\n');
                            fprintf('NO CONVERGENCE: MAXEVALS exceeded\n');
                            fprintf('Stage 1, increasing temperature\n');
                            fprintf('\nObj. fn. value %f\n', fopt);
                            fprintf('	   parameter	    search width\n');
                            for i =(1:n)
                                fprintf('%20f%20f\n', xopt(i), bounds(i));
                            end
                            fprintf('================================================\n');
                        end
                        varargout{1} = xopt;
                        varargout{2} = fopt;
                        varargout{3} = 0;
                        return
                    end
                    
                    %%  Accept the new point if the function value decreases
                    if(fp <= f)
                        x = xp;
                        f = fp;
                        nacc = nacc + 1;
                        nacp(h) = nacp(h) + 1;
                        nup = nup + 1;
                        
                        %%  If greater than any other point, record as new optimum.
                        if(fp < fopt)
                            xopt = xp;
                            fopt = fp;
                            nnew = nnew + 1;
                        end
                        
                        %% If the point is higher, use the Metropolis criteria to decide on
                        %% acceptance or rejection.
                    else
                        p = exp(-(fp - f) / t);
                        if(rand < p)
                            x = xp;
                            f = fp;
                            nacc = nacc + 1;
                            nacp(h) = nacp(h) + 1;
                            ndown = ndown + 1;
                        else
                            nrej = nrej + 1;
                        end
                    end
                end
            end
            %%  Adjust bounds so that approximately half of all evaluations are accepted.
            test = 0;
            for i = (1:n)
                ratio = nacp(i) / ns;
                if(ratio > 0.6)
                    bounds(i) = bounds(i) * (1.0 + 2.0 * (ratio - 0.6) / 0.4);
                elseif(ratio < 0.4)
                    bounds(i) = bounds(i) / (1.0 + 2.0 * ((0.4 - ratio) / 0.4));
                end
                %% keep within initial bounds
                if(bounds(i) >= (ub(i) - lb(i)))
                    test = test + 1; %% when this gets to n, we're done with fist stage
                    bounds(i) = ub(i) - lb(i);
                end
            end
            if ~converge
                converge = (test == n);
            end
        end
        if(verbosity == 1)
            fprintf('First stage: Increasing temperature to cover parameter space');
            fprintf('\nTemperature  %e', t);
            fprintf('\nmin function value so far %f', fopt);
            fprintf('\ntotal evaluations so far %d', func_evals);
            fprintf('\ntotal moves since temp change %d', nup + ndown + nrej);
            fprintf('\ndownhill  %d', nup);
            fprintf('\naccepted uphill %d', ndown);
            fprintf('\nrejected uphill %d', nrej);
            fprintf('\nout of bounds trials %d', lnobds);
            fprintf('\nnew minima this temperature %d', nnew);
            fprintf('\nparameter	search width\n');
            for i = (1:n)
                fprintf('%20f%20f\n', xopt(i), bounds(i));
            end
            fprintf('\n');
            %% Increase temperature quickly
            t = t*t;
%             for i =(neps:-1:2)
%                 fstar(i) = fstar(i-1);
%             end
            f = fopt;
            x = xopt;
        end
    end
    %% Second stage: temperature reduction loop
    converge = 0;
    while (converge==0)
        nup = 0;
        nrej = 0;
        nnew = 0;
        ndown = 0;
        lnobds = 0;
        
        %% repeat nt times then adjust temperature
        for m = (1:nt)
            %% repeat ns times, then adjust bounds
            for j=(1:ns)
                %% generate new point by taking last
                %% and adding a random value to each of elements,
                %% in turn
                for h=(1:n)
                    xp = x;
                    f_return = rand;
                    rand_draw = f_return;
                    xp(h) = x(h) + (2.0 * rand_draw - 1.0) * bounds(h);
                    if((xp(h) < lb(h)) || (xp(h) > ub(h)))
                        xp(h) = lb(h) + (ub(h) - lb(h)) * rand_draw;
                        lnobds = lnobds + 1;
                    end
                    
                    %% Evaluate function at new point
                    f_args(minarg) = xp;
                    c_args{2} = f_args;
                    f_return = feval(obj_fn,f_args);
                    fp = f_return;
                    func_evals = func_evals + 1;
                    
                    %% If too many function evaluations occur, terminate the algorithm
                    if(func_evals >= maxevals)
                        warning('NO CONVERGENCE: maxevals exceeded');
                        if(verbosity >= 1)
                            fprintf('\n================================================\n');
                            fprintf('SAMIN results\n');
                            fprintf('NO CONVERGENCE: MAXEVALS exceeded\n');
                            fprintf('Stage 2, decreasing temperature\n');
                            fprintf('\nObj. fn. value %f\n', fopt);
                            fprintf('	   parameter	    search width\n');
                            for i=(1:n)
                                fprintf('%20f%20f\n', xopt(i), bounds(i));
                            end
                            fprintf('================================================\n');
                        end
                        varargout{1} = xopt;
                        varargout{2} = fopt;
                        varargout{3} = 0;
                        return
                    end
                    
                    %%  Accept the new point if the function value decreases
                    if(fp <= f)
                        x = xp;
                        f = fp;
                        nacc = nacc + 1;
                        nacp(h) = nacp(h) + 1;
                        nup = nup + 1;
                        %%  If greater than any other point, record as new optimum
                        if(fp < fopt)
                            xopt = xp;
                            fopt = fp;
                            nnew = nnew + 1;
                        end           
                        
                        %% If the point is higher, use the Metropolis criteria to decide on
                        %% acceptance or rejection.
                    else
                        p = exp(-(fp - f) / t);
                        f_return = rand(1);
                        rand_draw = f_return;
                        if(rand_draw < p)
                            x = xp;
                            f = fp;
                            nacc = nacc + 1;
                            nacp(h) = nacp(h) + 1;
                            ndown = ndown + 1;
                        else
                            nrej = nrej + 1;
                        end
                    end
                end
            end
            
            %%  Adjust bounds so that approximately half of all evaluations are accepted
            for i = (1:n)
                ratio = nacp(i) / ns;
                if(ratio > 0.6)
                    bounds(i) = bounds(i) * (1.0 + 2.0 * (ratio - 0.6) / 0.4);
                elseif(ratio < .4)
                    bounds(i) = bounds(i) / (1.0 + 2.0 * ((0.4 - ratio) / 0.4));
                end
                %% keep within initial bounds
                if(bounds(i) > (ub(i) - lb(i)))
                    bounds(i) = ub(i) - lb(i);
                end
                nacp=zeros(n,1); % Possibly??
            end
            if(verbosity == 1)
                fprintf('\nIntermediate results before next temperature reduction');
                fprintf('\nTemperature  %e', t);
                fprintf('\nmin function value so far %f', fopt);
                fprintf('\ntotal evaluations so far %d', func_evals);
                fprintf('\ntotal moves since last temp reduction  %d', nup + ndown + nrej);
                fprintf('\ndownhill  %d', nup);
                fprintf('\naccepted uphill %d', ndown);
                fprintf('\nrejected uphill %d', nrej);
                fprintf('\nout of bounds trials %d', lnobds);
                fprintf('\nnew minima this temperature %d', nnew);
                fprintf('\n\n	       parameter	search width\n');
                for i=(1:n)
                    fprintf('%20f%20f\n', xopt(i), bounds(i));
                end
                fprintf('\n');
            end
            
            %% Check for convergence
            %% current function value must be within "tol"
            %% of last "neps" (an integer) function values,
            %% AND the last "neps" function values
            %% must be withing tol of overall best
            fstar(1) = f;
            test = 0;
            for i=(1:neps)
                test = test + (abs(f - fstar(i))>functol);
            end
            test = (test > 0); %% if different from zero, function conv. has failed
            test
            if (abs(fopt - f) <= functol) && (~test)
                %% check for bound narrow enough for parameter convergence
                converge = 1;
                for i=(1:n)
                    if(bounds(i) > paramtol)
                        converge = 0; %% no conv. if bounds too wide
                        break;
                    end
                end
            end
            
            %% check if too close to bounds, and change convergence message if so
            if (converge)
                if (lnobds > 0)
                    converge = 2;
                end
            end
            %% Are we done yet?
            if(converge>0)
                if(verbosity >= 1)
                    fprintf('\n================================================\n');
                    fprintf('SAMIN final results\n');
                    if (converge == 1)
                        fprintf('NORMAL CONVERGENCE\n\n');
                    end
                    if (converge == 2)
                        fprintf('WARNING: last point satisfies conv. criteria, \n\but is too close to bounds of parameter space\n');
                        fprintf('%f \% of last round evaluations out-of-bounds\n', 100*lnobds/(nup+ndown+nrej));
                        fprintf('Expand bounds and re-run\n');
                    end
                    fprintf('Func. tol. %e	Param. tol. %e\n', functol, paramtol);
                    fprintf('Obj. fn. value %f\n\n', fopt);
                    fprintf('	   parameter	    search width\n');
                    for i=(1:n)
                        fprintf('%20f%20f\n', xopt(i), bounds(i));
                    end
                    fprintf('================================================\n');
                end
                varargout{1} = xopt;
                varargout{2} = fopt;
                if (lnobds > 0)
                    converge = 2;
                end
                varargout{3} = converge;
                return
            end
            
            %% Reduce temperature, record current function value in the
            %% list of last "neps" values, and loop again
            t = rt * t;
            for i = neps:-1:2
                fstar(i) = fstar(i-1);
            end
            f = fopt;
            x = xopt;
        end
        varargout{1} = xopt;
        varargout{2} = fopt;
        varargout{3} = converge;
        return
    end
end