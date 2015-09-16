function [hout,hbout hfout]=areaplot(varargin)
% areaplot     : Plots a spec1d object with a gradient fill
% [PointsHandle ErrorBarHandles SurfHandels]=areaplot(s,zeroline,transparency)
% s            : spec1d object or array of objects
% zeroline     : [xstart ystart; xend yend] to describe sero line
amap=0;
pt=[];
for i=1:length(varargin)
    if isa(varargin{i},'spec1d')
        s(i)=varargin{i};
    elseif isa(varargin{i},'double')
        if length(varargin{i})>1
            pt=varargin{i};
        elseif length(varargin{i})==1
            if varargin{i}==1 || varargin{i}==0
                zg=varargin{i};
                npt=1000;
            else
                npt=varargin{i};
            end
        end
    elseif isa(varargin{i},'char')
        if amap==1
            sht=varargin{i};
        end
        if strcmp(varargin{i},'alpha') && amap==0
            amap=1;
        else
            amap=0;
        end
    end
end

% % Test Vairables
% s
% pt
% zg
% amap
% npt

set(0,'CurrentFigure',gcf)
hold on

k=1;
for i=1:length(s)
    x=s(i).x;   y=s(i).y;
    if isempty(pt)
        [xmin ind1]=min(x); [xmax ind2]=max(x);
        pt=[xmin y(ind1); xmax y(ind2)];
    end
    gr=diff(pt(:,2))/diff(pt(:,1));
    ic=pt(1,2)-pt(1,1)*gr;

    si(i)=interpolate(s(i),min(x):((max(x)-min(x))/(npt+1)):max(x));
    x=si(i).x;
    y=si(i).yfit;
    e=s(i).e;
    if isempty(y)
        warning('Spec1d object is not fitted! Using Y-data')
        y=si(i).y;
    end
    
    ylin=gr*x + ic;
    
    [X Y]=meshgrid(min(x):((max(x)-min(x))/npt):max(x),min(y):((max(y)-min(y))/npt):max(y));
    zb=linspace(-1,1,npt+1);
    
    stepx=floor(length(x)/(npt+1));
    Z=zeros(size(X));
    
    for l=1:length(X(:,1))
        z=  zb'.*(Y(:,l)<=y(((l-1)*stepx)+1)).*(Y(:,l)>=ylin(((l-1)*stepx)+1));
        z=z+zb'.*(Y(:,l)>=y(((l-1)*stepx)+1)).*(Y(:,l)<=ylin(((l-1)*stepx)+1));
        Z(:,l)=z';
        Z((Z(:,l))==0,l)=NaN;
        if ~zg
            Z(:,l)=Z(:,l)-ylin(((l-1)*stepx)+1);
        end
    end
    y=s(i).y;    x=s(i).x;
    tee=min(diff(x))/25;
    hout(i)=plot3(x,y,1.1*ones(1,length(x)),'Color','k','LineStyle','None','Marker','x','LineWidth',2);
    for j=1:length(x)
        hbout(k)=plot3([x(j) x(j)],[y(j)-e(j) y(j)+e(j)],1.1*ones(2,1),'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',1);
        k=k+1;
        hbout(k)=plot3([x(j)-tee x(j)+tee],[y(j)+e(j) y(j)+e(j)],1.1*ones(2,1),'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',1);
        k=k+1;
        hbout(k)=plot3([x(j)-tee x(j)+tee],[y(j)-e(j) y(j)-e(j)],1.1*ones(2,1),'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',1);
        k=k+1;
    end
    hfout(i)=surf(X,Y,Z,'LineStyle','None','FaceColor','interp','FaceLighting','phong','EdgeColor',[0.5 0.5 0.5]);
    if ~zg
        if amap
            alpha('color')
            set(hfout(i),'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
            alphamap(sht)
        else
            alpha('color')
            set(hfout(i),'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
            alphamap('rampup')
        end
    end
end
Cmap=[linspace(0,1,100)' zeros(100,2)];
colormap(gca,Cmap)
if ~zg
    set(gca,'CLim',[-1 1])
end
axis tight
view([0 90])
hold off