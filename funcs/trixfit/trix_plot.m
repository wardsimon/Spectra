function trix_plot(s,f,initialise,x_label,y_label,tit);
%
% TRIXFIT function to plot fit results plus parameters
%
% Des McMorrow 10 May 2003

%----- Create new figure

h=figure;
set(h,'units','centimeters')
set(h,'position',[2 2 15 20])    
axes('position',[0.1 0.5 0.4 0.4])
axis square

%----- Extract data from s and plot graph

x=getfield(s,'x');
y=getfield(s,'y');
err=getfield(s,'e');
yfit=getfield(s,'yfit');
x=x(:); y=y(:); err=err(:);

ytop = (y + err)';
ybot = (y - err)';
tee = (max(x)-min(x))/100;  % make tee .02 x-distance for error bars
%tee = 0;
xleft = (x-tee)';
xright = (x+tee)';
nnan=NaN*ones(size(x'));
xb=[xleft; xright; nnan; x'; x'; nnan; xleft; xright; nnan];
yb=[ybot; ybot; nnan; ybot; ytop; nnan; ytop; ytop; nnan];
n=9*size(xb,2);
xb=reshape(xb,n,1);
yb=reshape(yb,n,1);
hle=line(xb,yb,'color','b','linewidth',2);

hll=line(x,y,'color','b',...
         'LineStyle','None','Marker','o',...
         'MarkerSize',7,'MarkerFaceColor','b');

hll=line(x,yfit,'color','r',...
         'LineStyle','-','Marker','none',...
         'LineWidth',2);

xlabel(x_label)
ylabel(y_label)
title(tit)

grid on; box on

%----- Create label arrays and print

axes('position',[0 0 1 1])
axis off

npars=size([f.pnames],1);

pnames=f.pnames;
pvals=f.pvals;
evals=f.evals;

h=text(0,0,'x','position',[0.65 0.70]);
set(h,'units','points');
basepos=get(h,'position');
delete(h);

fontheight=12;

mat=strvcat(...
[sprintf('%10s','Method') '=' sprintf('%s',initialise.resolution_method)],...     
[sprintf('%10s','RescalPars') '=' sprintf('%s',initialise.rescal_pars)],...          
[sprintf('%10s','InstConfig') '=' sprintf('%s',initialise.popovici_pars)],...
' ',...
[sprintf('%10s','XsecFile') '=' sprintf('%s',initialise.xsec_file)],...
[sprintf('%10s','BckgFile') '=' sprintf('%s',initialise.bkgd_file)],...          
[sprintf('%10s','FitPNames') '=' sprintf('%s',initialise.pnam_file)],...            
[sprintf('%10s','ICorrFile') '=' sprintf('%s',initialise.corr_file)],...
' ',...
[sprintf('%10s','MonFlag') '=' sprintf('%i',initialise.monitor_flag)],...                   
[sprintf('%10s','MCSamples') '=' sprintf('%7.2e',initialise.monte_carlo_samples)],...     
' '...
);

for i=1:4
   str=[sprintf('%10s',pnames(i,1:end)) '=' sprintf('%7.2e',pvals(i))];   
   mat=strvcat(mat,str);
end

mat=strvcat(mat,' ');

for i=5:npars
   str=[sprintf('%10s',pnames(i,1:end)) '=' sprintf('%7.2e',pvals(i)) '+-' sprintf('%7.2e',evals(i))];   
   mat=strvcat(mat,str);
end

mat=strvcat(mat,' ');
mat=strvcat(mat,[sprintf('%10s','Chisq') '=' sprintf('%7.2e',f.chisq)]);

nstrs=size(mat,1);

for i=1:nstrs

pos=basepos(1:2)-[0 (i-1)*fontheight/1.1];
h=text(0,0,mat(i,1:end),'units','points','position',pos,'Fontsize',fontheight,'fontname','fixedwidth','interpreter','none');

end   

set(gcf,'Name',tit)
wysiwyg