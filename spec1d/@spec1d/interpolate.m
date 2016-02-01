function s_out = interpolate(s_in,x_new,varargin)
%
% function s_out = interpolate(s1,s2,...x_new)
%
% @SPEC1D/INTERPOLATE function to interpolate a given spectra to a new x-axis.
%
% Simon Ward 01/02/2016

p = inputParser;
p.addRequired('s',@(x)isa(x,'spec1d'));
p.addRequired('x_new',@(x) (isnumeric(x) && isreal(x)) || (isa(x,'spec1d') && length(x)==1))
p.addOptional('order',-1,@(x) (isnumeric(x) && isreal(x)))
p.addOptional('method','weightedpoly',@ischar)
p.KeepUnmatched = 1;

p.parse(s_in,x_new,varargin{:})

s = p.Results.s;
x_new = p.Results.x_new;
method = p.Results.method;
varargin = struct2cell(p.Unmatched);

c = onCleanup(@() warning('on','MATLAB:nearlySingularMatrix'));
warning('off','MATLAB:nearlySingularMatrix')
for i = 1:length(s)
    
    x = s(i).x(:);
    y = s(i).y(:);
    e = s(i).e(:);
    
    switch lower(method)
        case 'weightedpoly'
            if p.Results.order < 1
                var_poly = zeros(length(x),1);
                for j = 1:(length(x)-1)
                    [pp, S] = sdinterp.weightedpoly(x,y,e,j);
                    [~, e_new] = polyval(pp,x_new,S);
                    var_poly(j) = sum((e_new).^2)/(length(x)-j-1);
                    if j > 2
                        if var_poly(j) > 1E3 || abs(diff(var_poly(j-1:j))) > 10
                            break
                        end
                    end
                end
                var_poly(var_poly==0) = NaN;
                [~, ind] = min(var_poly);
                [pp, S] = sdinterp.weightedpoly(x,y,e,ind);
                [y_new, e_new] = polyval(pp,x_new,S);
            else
                [pp, S] = sdinterp.weightedpoly(x,y,e,p.Results.order);
                [y_new, e_new] = polyval(pp,x_new,S);
            end
        case 'builtin'
            y_new = interp1(x,y,x_new);
            e_new = sqrt(interp1(x,e.^2,x_new));
        otherwise
            try
                y_new = feval(sprintf('sdinterp.%s',method),x,y,x_new,varargin{:});
                e_new = sqrt(interp1(x,e.^2,x_new));
            catch
                help('sdinterp')
                error('spec1d:interpolate:ErrorInInterpolation','The function %s has thrown an error or doesnt exist. See documentation',method)
            end
    end
    
    % Make output object
    r = s(i);
    r.x = x_new;
    r.y = y_new;
    r.e = e_new;
    s_out(i) = spec1d(r);
end