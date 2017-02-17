function varargout = multifit_ini(s,pin,flag)
%% Multifit initialisation procedure
% s = arrray of spec1d objects
% pin = parameters for fitting
% flag = 0, 1, 2 
% sep = separation of spec1d objects.
% Based on the work of MM


%----- Define global variables

global  x_per_spec param_keep

% if length(sep)~=length(s)
%     error('There must be the same number of spec1d objects as separation points')
% end
% 
% s_sep = sep;

%----- Create a big spec1d file that contains, in order, all the smaller
%      ones. Also, memorize the number of points per spec1d file so that
%      this operation stays reversible!

if nargin == 2 && isa(pin,'specfit')
   f_in = pin;
   pin = f_in.pvals;
   flag = f_in.notfixed;
else
    f_in = [];
end

x_per_spec = zeros(1,length(s));

x = []; y = []; e = [];

for il = 1:length(s)
    [xil, yil, eil] = extract(s(il)); 
%----- Remove zeros from e
    xil(eil==0) = []; 
    yil(eil==0) = []; 
    eil(eil==0) = [];
    x_per_spec(il) = length(xil);
%----- Fold all data in the same array
    x((length(x)+1):length(x)+length(xil)) = xil(:);
    y((length(y)+1):length(y)+length(yil)) = yil(:);
    e((length(e)+1):length(e)+length(eil)) = eil(:);
end

% Handle the basic data.
s_out = s(1).copy;
s_out.x = x;
s_out.y = y;
s_out.e = e;
s_out.yfit = [];

% This deals with additional data in each spectra.
add_dat = {s.userdata};
for i = 1:length(add_dat)
    add_dat{i} = rmfield(add_dat{i},'rind');
    if isfield(add_dat{i},'combine')
        if i == 1
            ns = mergestruct(rmfield(add_dat{i},'combine'),add_dat{i}.combine);
        else
            f = fieldnames(add_dat{i}.combine);
           for j = 1:length(f)
               if isfield(ns,f{j})
                    ns.(f{j}) = vertcat(ns.(f{j}),add_dat{i}.combine.(f{j}));
               else
                   ns.(f{j}) = add_dat{i}.combine.(f{j});
               end
           end
        end
    end
end
ns.rind = s(1).userdata.rind;
s_out.userdata = ns;

s_out = feval(class(s(1)),s_out);

%---- Take care of the parameters!
pin_old=pin(:); 
dpin_old=flag(:);

%----- dpin now has 3 possible states: 0=fixed, 1=the same for all spec1d
%      objects, 2=fitted independently in all spec1d objects. I have now to
%      create the actual pin and dp vectors

%----- Also, param_keep will store the order of the parameters

pout = zeros(1,(length(pin)+length(find(flag>1))*(length(s)-1)));
dpin = pout;
param_keep = dpin;

ll = 1;
for j = 1:length(dpin_old)
    if dpin_old(j) == 2
        pout(ll:(ll+length(x_per_spec)-1)) = pin_old{j};
        dpin(ll:(ll+length(x_per_spec)-1)) = 1;
        param_keep(ll:(ll+length(x_per_spec)-1)) = j;
        ll = ll+length(x_per_spec);
    elseif dpin_old(j) == 3
        pout(ll:(ll+length(x_per_spec)-1)) = pin_old{j};
        dpin(ll:(ll+length(x_per_spec)-1)) = 0;
        param_keep(ll:(ll+length(x_per_spec)-1)) = j;
        ll = ll+length(x_per_spec);
    elseif (dpin_old(j) == 1) || (dpin_old(j) == 0)
        pout(ll) = pin_old{j};
        dpin(ll) = dpin_old(j);
        param_keep(ll) = j;
        ll = ll+1;
    else
        disp('Please use following notation: 0 = fixed, 1 = same for all objects, 2 = independently fitted in all objects, 3 = fixed and different for all objects')
    end
end

if ~isempty(f_in)
    f_in = f_in.copy;
    f_in.specID = s_out.ident;
    f_in.pvals = pout;
    f_in.evals = nan(size(pout));
    f_in.notfixed = dpin;
    f_in.userdata.param_keep = param_keep;
    f_in.userdata.x_per_spec = x_per_spec;
end

if nargout == 2
    varargout{1} = s_out;
    varargout{2} = f_in;
else
    varargout{1} = s_out;
    varargout{2} = pout;
    varargout{3} = dpin;
end
