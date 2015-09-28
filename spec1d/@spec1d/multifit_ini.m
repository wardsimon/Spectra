function [sout, pout, dpin]=multifit_ini(s,pin,flag,sep)
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

x_per_spec=zeros(1,length(s));

x=[]; y=[]; e=[];

for il=1:length(s)
    [xil, yil, eil] = extract(s(il)); 
%----- Remove zeros from e
    xil(eil==0)=[]; 
    yil(eil==0)=[]; 
    eil(eil==0)=[];
    x_per_spec(il)=length(xil);
%----- Fold all data in the same array
    x((length(x)+1):length(x)+length(xil))=xil(:);
    y((length(y)+1):length(y)+length(yil))=yil(:);
    e((length(e)+1):length(e)+length(eil))=eil(:);
end

sout=spec1d(x,y,e);

%---- Take care of the parameters!
pin_old=pin(:); 
dpin_old=flag(:);

%----- dpin now has 3 possible states: 0=fixed, 1=the same for all spec1d
%      objects, 2=fitted independently in all spec1d objects. I have now to
%      create the actual pin and dp vectors

%----- Also, param_keep will store the order of the parameters

pout=zeros(1,(length(pin)+length(find(flag>1))*(length(s)-1)));
dpin=pout;
param_keep=dpin;

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

