function [sout, vout] = multifit_extract(s,w,dpin)

global x_per_spec param_keep

wout = cell(param_keep(length(param_keep)),1);
eout = wout;

for j = 1:param_keep(length(param_keep))
    index = find(param_keep == j);
    wout{j}=w.pvals(index);
    eout{j}=w.evals(index); 
end


[x, y, e]=extract(s);
yfit=getfield(s,'yfit');

for il=1:length(x_per_spec)
    
    if il == 1
       ind = 1:x_per_spec(1); 
    else
        ind = (sum(x_per_spec(1:(il-1)))+1):sum(x_per_spec(1:(il)));
    end

    sloop.x=x(ind);
    sloop.y=y(ind);
    sloop.e=e(ind);
    sloop.yfit=yfit(ind);
    
    sout(il)=spec1d(sloop);
    pout=zeros(1,param_keep(length(param_keep)));
    eout=zeros(1,param_keep(length(param_keep)));
        
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
 
    vout(il).pvals    = pout(:);
    vout(il).evals    = eout(:);
    vout(il).function = w.function;
    vout(il).pnames   = w.pnames{((il-1)*(length(w.pnames)/length(x_per_spec))+1):(il*length(w.pnames)/length(x_per_spec))};
    % I think this is the correct correction!
    vout(il).chisq    = (length(w.pvals)-sum(dpin==0))*w.chisq/param_keep(length(param_keep));
    vout(il).rsq      = w.rsq;
end