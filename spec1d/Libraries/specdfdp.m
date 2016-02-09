function prt=specdfdp(x,f,p,dp,func)
% MFIT function prt=mf_dfdp(x,f,p,dp,func)
%	Called by mf_lsqr.m
%
m=length(x);
n=length(p);                % dimensions
ps=p; 
prt=zeros(m,n);     % Jacobian
del=zeros(n,1);    % initialise Jacobian to Zero
for j=1:n
    del(j)=dp(j) .*p(j);                %cal delx=fract(dp)*param value(p)
    if p(j)==0
        del(j)=dp(j);                   %if param=0 delx=fraction
    end
    p(j)=ps(j) + del(j);
    if del(j)~=0, 
        f1=feval(func,x,p);
        if dp(j) < 0,
            prt(:,j)=(f1-f)./del(j);
        else
            p(j)=ps(j)- del(j);
            prt(:,j)=(f1-feval(func,x,p))./(2 .*del(j));
        end
    end
    p(j)=ps(j);                        %restore p(j)
end
return
