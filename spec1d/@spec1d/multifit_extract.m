function [sout, vout] = multifit_extract(s,w)

global x_per_spec param_keep

wout = cell(param_keep(length(param_keep)),1);
eout = wout;
dout = zeros(1,param_keep(length(param_keep)));

for j = 1:param_keep(length(param_keep))
    index = find(param_keep == j);
    wout{j} = w.pvals(index);
    eout{j} = w.evals(index); 
    dout(j) = w.notfixed(index(1));
end

try
    x = s.x;
    y = s.y;
    e = s.e;
    yfit = s.yfit;
    rind = s.userdata.rind;
    x = x(rind);
    y = y(rind);
    e = e(rind);
    yfit = yfit(rind);
catch
    [x, y, e]=extract(s);
    yfit=getfield(s,'yfit');
end

for il=1:length(x_per_spec)
    
    if il == 1
       ind = 1:x_per_spec(1); 
    else
        ind = (sum(x_per_spec(1:(il-1)))+1):sum(x_per_spec(1:(il)));
    end

    sloop.x = x(ind);
    sloop.y = y(ind);
    sloop.e = e(ind);
    sloop.yfit = yfit(ind);
    
    sout(il) = feval(class(s),sloop);
    pout = zeros(1,param_keep(length(param_keep)));
    eout = zeros(1,param_keep(length(param_keep)));
       
    for h = 1:param_keep(length(param_keep))
        index = (find(param_keep == h));
        if length(index) == 1
            pout(h) = w.pvals(index);
            eout(h) = w.evals(index);
        else
            pp      = w.pvals(index);
            pout(h) = pp(il);
            ee      = w.evals(index);
            eout(h) = ee(il);
        end
    end
    
    vout(il) = specfit(pout,eout,w.func,dout,sout(il).ident,w.pnames);
%         results = struct('pvals',p,'evals',sig,'func',func,'pnames',pnames,'chisq',ChiSq,'rsq',RSq,'notfixed',notfixed);
%         fitdata(il) = specfit(results);
    sout(il) = sout(il).setfitdata(vout(il));
 
%     vout(il).pvals    = pout;
%     vout(il).evals    = eout;
%     vout(il).func = w.func;
%     vout(il).pnames   = w.pnames;
%     % I think this is the correct correction!
%     vout(il).chisq    = (length(w.pvals)-sum(dpin==0))*w.chisq/param_keep(length(param_keep));
%     vout(il).rsq      = w.rsq;
end