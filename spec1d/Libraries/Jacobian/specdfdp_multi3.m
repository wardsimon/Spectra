function prt = specdfdp_multi3 (s_in,f_in,f,add_dat)
% numerical partial derivatives (Jacobian) df/dp for use with leasqr
% --------INPUT VARIABLES---------
% x=vec or matrix of indep var(used as arg to func) x=[x0 x1 ....]
% f=func(x,p) vector initialsed by user before each call to dfdp
% p= vec of current parameter values
% dp= fractional increment of p for numerical derivatives
%      dp(j)>0 central differences calculated
%      dp(j)<0 one sided differences calculated
%      dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
% func=string naming the function (.m) file
%      e.g. to calc Jacobian for function expsum prt=dfdp(x,f,p,dp,'expsum')
% bounds=two-column-matrix of lower and upper bounds for parameters
%      If no 'bounds' options is specified to leasqr, it will call
%      dfdp without the 'bounds' argument.
%----------OUTPUT VARIABLES-------
% prt= Jacobian Matrix prt(i,j)=df(i)/dp(j)
%================================

global param_keep x_per_spec multifit_ind

x = s_in.x;
func = f_in.func;
bounds = f_in.fitdata.Bounds;
dp = -(f_in.notfixed(:)>0)*f_in.fitdata.TolX;

len = f_in.lengthP();
f_out = f_in.copy;

mf_f = add_dat.mf_f;
add_dat = rmfield(add_dat,'mf_f');

m=size(x,1); if (m==1), m=size(x,2); end  %# PAK: in case #cols > #rows
n=len;      %dimensions

if length(dp) ~= len
    error('Fixed and Pin need to be the same length')
end
prt=zeros(m,n);       % initialise Jacobian to Zero
del = dp .* f_in.pvals; %cal delx=fract(dp)*param value(p)
idx = f_in.pvals == 0;
del(idx) = dp(idx); %if param=0 delx=fraction
idx = dp > 0;
del(idx) = abs (del(idx)); % not for one-sided intervals, changed
% direction of intervals could change
% behavior of optimization without bounds
min_del = min (abs (del), bounds(:, 2) - bounds(:, 1));

if isempty(x_per_spec)
    for j=1:n
        ps = f_in.pvals;
        if (dp(j)~=0)
            if (dp(j) < 0)
                ps(j) = f_in.pvals(j) + del(j);
                if (ps(j) < bounds(j, 1) || ps(j) > bounds(j, 2))
                    t_del1 = max (bounds(j, 1) - f_in.pvals(j), - abs (del(j))); %
                    %non-positive
                    t_del2 = min (bounds(j, 2) - f_in.pvals(j), abs (del(j))); %
                    %non-negative
                    if (- t_del1 > t_del2)
                        del(j) = t_del1;
                    else
                        del(j) = t_del2;
                    end
                    ps(j) = p(j) + del(j);
                end
                f_out.pvals = ps;
                if mf_f
                    f1 = feval(func,x,f_out.p_vals);
                else
                    f1 = feval(func,x,f_out.p_vals,add_dat);
                end
                prt(:, j) = (f1(:) - f) / del(j);
            else
                if (f_in.pvals(j) - del(j) < bounds(j, 1))
                    tp = bounds(j, 1);
                    ps(j) = tp + min_del(j);
                elseif (f_in.pvals(j) + del(j) > bounds(j, 2))
                    ps(j) = bounds(j, 2);
                    tp = ps(j) - min_del(j);
                else
                    ps(j) = f_in.pvals(j) + del(j);
                    tp = f_in.pvals(j) - del(j);
                    min_del(j) = 2 * del(j);
                end
                f_out.pvals = ps;
                if mf_f
                    f1 = feval(func,x,f_out.p_vals);
                else
                    f1 = feval(func,x,f_out.p_vals,add_dat);
                end
                ps(j) = tp;
                f_out.pvals = ps;
                if mf_f
                    f2 = feval(func,x,f_out.p_vals);
                else
                    f2 = feval(func,x,f_out.p_vals,add_dat);
                end
                prt(:, j) = (f1(:) - f2(:)) / min_del(j);
            end
        end
    end
else
    for j=1:param_keep(length(param_keep))
        ps = f_in.pvals;
        par_ind=find(param_keep==j);
        del_temp = dp(par_ind) .* f_in.pvals(par_ind);        %cal delx=fract(dp)*param value(p)
        idx = f_in.pvals(par_ind) == 0;
        dp_temp = dp(par_ind);
        del_temp(idx) = dp_temp(idx); %if param=0 delx=fraction
        idx = dp(par_ind) > 0;
        del_temp(idx) = abs (del_temp(idx));
        del = del_temp;
        if all(dp(par_ind)~=0)
            % Now split into singular and multi
            if length(par_ind) == 1
                if (dp(par_ind) < 0)
                    ps(par_ind) = f_in.pvals(par_ind) + del(1);
                    if (ps(par_ind) < bounds(par_ind, 1) || ps(par_ind) > bounds(par_ind, 2))
                        t_del1 = max (bounds(par_ind, 1) - f_in.pvals(par_ind), - abs (del(1))); %
                        %non-positive
                        t_del2 = min (bounds(par_ind, 2) - f_in.pvals(par_ind), abs (del(1))); %
                        %non-negative
                        if (- t_del1 > t_del2)
                            del(1) = t_del1;
                        else
                            del(1) = t_del2;
                        end
                        ps(par_ind) = f_in.pvals(par_ind) + del(1);
                    end
                    f1 = zeros(sum(x_per_spec),1);
                    for i = 1:length(x_per_spec)
                        [p_new, ~, ind] = multifitp2p(ps,zeros(size(ps)),i);
                        multifit_ind = ind;
                        f_out.pvals = p_new;
                        add_dat.multifit_ind = ind;
                        if mf_f
                            f1(ind) = feval(func,x(ind),f_out.p_vals);
                        else
                            f1(ind) = feval(func,x(ind),f_out.p_vals,add_dat);
                        end
                    end
                    prt(:, par_ind) = (f1 - f) / del(1);
                else
                    if (f_in.pvals(par_ind) - del(1) < bounds(par_ind, 1))
                        tp = bounds(par_ind, 1);
                        ps(par_ind) = tp + min_del(par_ind);
                    elseif (f_in.pvals(par_ind) + del(1) > bounds(par_ind, 2))
                        ps(par_ind) = bounds(par_ind, 2);
                        tp = ps(par_ind) - min_del(par_ind);
                    else
                        ps(par_ind) = f_in.pvals(par_ind) + del(1);
                        tp = f_in.pvals(j) - del(1);
                        min_del(par_ind) = 2 * del(1);
                    end
                    f1 = zeros(sum(x_per_spec),1);
                    for i = 1:length(x_per_spec)
                        [p_new, ~, ind] = multifitp2p(ps,zeros(size(ps)),i);
                        f_out.pvals = p_new;
                        multifit_ind = ind;
                        add_dat.multifit_ind = ind;
                        if mf_f
                            f1(ind) = feval(func,x(ind),f_out.p_vals);
                        else
                            f1(ind) = feval(func,x(ind),f_out.p_vals,add_dat);
                        end
                    end
                    ps(par_ind) = tp;
                    f2 = zeros(sum(x_per_spec),1);
                    for i = 1:length(x_per_spec)
                        [p_new, ~, ind] = multifitp2p(ps,zeros(size(ps)),i);
                        multifit_ind = ind;
                        f_out.pvals = p_new;
                        add_dat.multifit_ind = ind;
                        if mf_f
                            f2(ind) = feval(func,x(ind),f_out.p_vals);
                        else
                            f2(ind) = feval(func,x(ind),f_out.p_vals,add_dat);
                        end
                    end
                    prt(:, par_ind) = (f1 - f2) / min_del(par_ind);
                end
            else
                lo_ind=1;
                for jj=1:length(par_ind)
                    jjj=par_ind(jj);
                    ma_ind=sum(x_per_spec(1:jj));
                    if (dp(jjj) < 0)
                        ps(jjj) = p(jjj) + del(jj);
                        if (ps(jjj) < bounds(jjj, 1) || ps(jjj) > bounds(jjj, 2))
                            t_del1 = max (bounds(jjj, 1) - f_in.pvals(jjj), - abs (del(jj))); %
                            %non-positive
                            t_del2 = min (bounds(jjj, 2) - f_in.pvals(jjj), abs (del(jj))); %
                            %non-negative
                            if (- t_del1 > t_del2)
                                del(jj) = t_del1;
                            else
                                del(jj) = t_del2;
                            end
                            ps(jjj) = f_in.pvals(jjj) + del(jj);
                        end
                        f1 = zeros(sum(x_per_spec),1);
                        for i = 1:length(x_per_spec)
                            [p_new, ~, ind] = multifitp2p(ps,zeros(size(ps)),i);
                            multifit_ind = ind;
                            f_out.pvals = p_new;
                            add_dat.multifit_ind = ind;
                            if mf_f
                                f1(ind) = feval(func,x(ind),f_out.p_vals);
                            else
                                f1(ind) = feval(func,x(ind),f_out.p_vals,add_dat);
                            end
                        end
                        prt((lo_ind:ma_ind), jjj) = (f1(lo_ind:ma_ind) - f(lo_ind:ma_ind)) / del(jj);
                    else
                        if (p(jjj) - del(jj) < bounds(jjj, 1))
                            tp = bounds(jjj, 1);
                            ps(jjj) = tp + min_del(jjj);
                        elseif (p(par_ind) + del(jjj) > bounds(jjj, 2))
                            ps(jjj) = bounds(jjj, 2);
                            tp = ps(jjj) - min_del(jjj);
                        else
                            ps(jjj) = p(jjj) + del(jjj);
                            tp = p(j) - del(jjj);
                            min_del(jjj) = 2 * del(jjj);
                        end
                        f1 = zeros(sum(x_per_spec),1);
                        for i = 1:length(x_per_spec)
                            [p_new, ~, ind] = multifitp2p(ps,zeros(size(ps)),i);
                            multifit_ind = ind;
                            f_out.pvals = p_new;
                            add_dat.multifit_ind = ind;
                            if mf_f
                                f1(ind) = feval(func,x(ind),f_out.p_vals);
                            else
                                f1(ind) = feval(func,x(ind),f_out.p_vals,add_dat);
                            end
                        end
                        ps(jjj) = tp;
                        f2 = zeros(sum(x_per_spec),1);
                        for i = 1:length(x_per_spec)
                            [p_new, ~, ind] = multifitp2p(ps,zeros(size(ps)),i);
                            f_out.pvals = p_new;
                            multifit_ind = ind;
                            add_dat.multifit_ind = ind;
                            if mf_f
                                f2(ind) = feval(func,x(ind),f_out.p_vals);
                            else
                                f2(ind) = feval(func,x(ind),f_out.p_vals,add_dat);
                            end
                        end
                        prt((lo_ind:ma_ind), jjj) = (f1(lo_ind:ma_ind) - f2(lo_ind:ma_ind)) / min_del(jjj);
                    end
                end
            end
        end
    end
end
end