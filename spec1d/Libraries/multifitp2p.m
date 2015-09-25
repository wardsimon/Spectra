function [pout, dout, ind] = multip2p(pin,din,i)


global param_keep x_per_spec

if i == 1
    ind = 1:x_per_spec(1);
else
    ind = (sum(x_per_spec(1:(i-1)))+1):(sum(x_per_spec(1:i)));
end

% Initialises pout so that it is as long as the original p-vector

pout=zeros(1,param_keep(length(param_keep)));
dout=zeros(1,param_keep(length(param_keep)));

for h=1:param_keep(length(param_keep))
    index=(find(param_keep==h));
    if length(index)==1
        pout(h)=pin(index);
        dout(h)=din(index);
    else
        pp=pin(index);
        pout(h)=pp(i);
        dd=din(index);
        dout(h)=dd(i);
    end
end

end