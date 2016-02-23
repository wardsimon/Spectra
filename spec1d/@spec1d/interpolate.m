function s_out = interpolate(s_in,x_new,varargin)
%
% function s_out = interpolate(s1,x_new,options)
%
% @SPEC1D/INTERPOLATE function to interpolate a given spectra to a new x-axis.
%
% Interpolate the given spectra for the x-ranges given by the vector x_new.
%
% Depending on the optional 'method', points are interpolated as
% 'weightedpoly' : Adaptive polynomial interpolation where errors are calculated and
%                  minimised for a n-length polynomial function.
% 'linear'	     : Nearest neighbout interpolation where values are
%                  weighted by inverse distance (including errors)
% 'builtin'	     : Use the MATLAB interp1 function. Optional arguments can
%                  be passed
% 'FUNCTION'     : The value FUNCTION is a function which can be explained
%                  with the documentation "help('sdinterp')". Optional
%                  arguments can be passed.
% Default is 'linear'. 'builtin' or 'linear' is suggested for speed.
%
% If the method 'weightedpoly' is selected, the optional 'order', paramter
% is available. The value is a positive integer and corresponds to the
% order of the fitted polynomial. Default is -1, which optimises for least
% error.
%
% s1 can be single spectra or arrays of spectra.
% x_new is a vecor of arbitary length. Where min(x_new) < min(s1.x) and
% max(x_new) < max(s1.x)
%
% Example:
% Interpolate s1 to vector -4:0.1:4.
% >r = interpolate(s1,-4:0.1:4) % standard interpolate
% >r = interpolate(s1,-4:0.1:4,'order',5) % interpolate with a 5th order
%                                           polynomial
% >r = interpolate(s1,-4:0.1:4,'method','linear') % linear interpolation
% >r = interpolate(s1,-4:0.1:4,'method','builtin')% linear interpolation
%                                                   using the builtin MATLAB function
% >r = interpolate(s1,-4:0.1:4,'method','cbezier')% Interpolation using a
%                                                   function in the
%                                                   sdinterp library
% !!Note!!
% For data where e is the sqrt of y this is OK. For arbitary data, it is not
% correct but an approximation.
%
% Simon Ward 02/02/2016

p = inputParser;
p.addRequired('s',@(x)isa(x,'spec1d'));
p.addRequired('x_new',@(x) (isnumeric(x) && isreal(x)) || (isa(x,'spec1d') && length(x)==1))
p.addOptional('order',-1,@(x) (isnumeric(x) && isreal(x)))
p.addOptional('method','linear',@ischar)
p.KeepUnmatched = 1;

p.parse(s_in,x_new,varargin{:})

s = p.Results.s;
x_new = p.Results.x_new;
method = p.Results.method;
varargin = struct2cell(p.Unmatched);

for i = 1:length(s)
    
    x = s(i).x(:);
    y = s(i).y(:);
    e = s(i).e(:);
    
    if (min(x_new) < min(x)) || (max(x) < max(x_new))
        warning('spec1d:interpolate:MinMaxXnewValues','The supplied interpolation range is out of the spectra range.')
        if any(strcmpi(method,{'linear','weightedpoly'}))
            x_new(x_new>max(x) | x_new<min(x)) = [];
        end
    end
    switch lower(method)
        case 'linear'
            temp = mat2cell(bsxfun(@minus,x(:,ones(1,length(x_new))),x_new),length(x),ones(size(x_new)));
            [rind1,~] = cellfun(@(x) find(x <= 0,1,'last'),temp);
            [rind2,~] = cellfun(@(x) find(x >= 0,1,'first'),temp);
            y_new = arrayfun(@(x_new,n1,n2) y(n1)*(x(n2)-x_new)/(x(n2)-x(n1))+...
                y(n2)*(x(n1)-x_new)/(x(n1)-x(n2)),x_new,rind1,rind2);
            e_new = arrayfun(@(x_new,n1,n2) sqrt((e(n1)*(x(n2)-x_new)/(x(n2)-x(n1)))^2+...
                (e(n2)*(x(n1)-x_new)/(x(n1)-x(n2)))^2),x_new,rind1,rind2);
        case 'weightedpoly'
            if p.Results.order < 1
                c = onCleanup(@() [warning('on','MATLAB:nearlySingularMatrix'), warning('on','MATLAB:polyval:ZeroDOF')]);
                warning('off','MATLAB:nearlySingularMatrix')
                warning('off','MATLAB:polyval:ZeroDOF')
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
                [y_new, e_new] = polyval(pp,x_new,S); % This is the error for the interpolation. We need to add on the error for the points.
            else
                [pp, S] = sdinterp.weightedpoly(x,y,e,p.Results.order);
                [y_new, e_new] = polyval(pp,x_new,S);
            end
            temp = mat2cell(bsxfun(@minus,x(:,ones(1,length(x_new))),x_new),length(x),ones(size(x_new)));
            [rind1,~] = cellfun(@(x) find(x <= 0,1,'last'),temp);
            [rind2,~] = cellfun(@(x) find(x >= 0,1,'first'),temp);
            e_new = e_new + arrayfun(@(x_new,n1,n2) sqrt((e(n1)*(x(n2)-x_new)/(x(n2)-x(n1)))^2+...
                (e(n2)*(x(n1)-x_new)/(x(n1)-x(n2)))^2),x_new,rind1,rind2);
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
    
    % Remove all bad data
    rind = isnan(y_new) | isnan(e_new) | isinf(y_new)| isinf(e_new);
    
    % Check for imaginary
    if any(~isreal(y_new) | ~isreal(y_new))
        warning('spec1d:interpolate:ImaginaryReturnValues','The result of interpolation is imaginary. Taking the real value')
        y_new = real(y_new);
        e_new = real(e_new);
    end
    
    % Make output object
    r = s(i);
    r.x = x_new(~rind);
    r.y = y_new(~rind);
    r.e = e_new(~rind);
    s_out(i) = spec1d(r);
end