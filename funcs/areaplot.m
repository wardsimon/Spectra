function [p1 p2]=areaplot(x,y,pt,zg,sh,sht)

gr=diff(pt(:,2))/diff(pt(:,1));
ic=pt(1,2)-pt(1,1)*gr;
ylin=gr*x + ic;
xcpt=x(round(length(x)/2));
ycpt=ylin(round(length(x)/2));

[X Y]=meshgrid(min(x):((max(x)-min(x))/1000):max(x),min(y):((max(y)-min(y))/1000):max(y));
zb=linspace(-1,1,1000+1);

stepx=floor(length(x)/1000);
Z=zeros(size(X));

for i=1:length(X(:,1))
    z=  zb'.*(Y(:,i)<=y(((i-1)*stepx)+1)).*(Y(:,i)>=ylin(((i-1)*stepx)+1));
    z=z+zb'.*(Y(:,i)>=y(((i-1)*stepx)+1)).*(Y(:,i)<=ylin(((i-1)*stepx)+1));
    Z(:,i)=z';
    Z((Z(:,i))==0,i)=NaN;
    if ~zg
                Z(:,i)=Z(:,i)-ylin(((i-1)*stepx)+1);
    end
end

set(0,'CurrentFigure',gcf)
hold on
p1=plot3(x,y,1.1*ones(1,length(x)),'k-','LineWidth',2);
p2=pcolor(X,Y,Z,'LineStyle','None','FaceColor','interp','FaceLighting','phong');
if ~zg
    set(gca,'CLim',[-max(abs(ylin)) max(abs(ylin))])
end
axis tight
view([0 90])

if nargin==6
    if strcmp(sh,'alpha')
        alpha('color')
        set(p2,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
        alphamap(sht)
    end
elseif nargin==5
    if strcmp(sh,'alpha')
        alpha('color')
        set(p2,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
        alphamap('rampup')
    end
end

hold off