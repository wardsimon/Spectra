%%%%%%%%%%%%%%%%%%%%

addpath('C:\Arbeitsverzeichnis\Experimente\2008\LiHoErF4_RITAII\neuhdf')

%%%(100)
files=[6090];
n=1;
fn=['hdf_data\rita22008n00',num2str(files(n)),'.hdf'];
data=readhdfdata(fn);
Magf(n)=double(getfield(data,'magnetic_field'));
Temp(n)=mean(double(getfield(data,'sample_temperature')));
ScanCommand=getfield(data,'scancommand');

s=hdf2spec(fn,'Ql',0);
[sf,f]=fits(s,'lorzs',[10 0 1 0 1],[1 1 1 1 1])
plot(sf(5))
formatpars(f(5),2)
title([num2str(files(n)),' ; ',ScanCommand])


%%%(100) diff
files=[6090 6103];
for n=1:2
    fn=['hdf_data\rita22008n00',num2str(files(n)),'.hdf'];
    data=readhdfdata(fn);
    Magf(n)=double(getfield(data,'magnetic_field'));
    Temp(n)=mean(double(getfield(data,'sample_temperature')));
    ScanCommand=getfield(data,'scancommand');
    s(n,:)=hdf2spec(fn,'Ql',0);
end
for n=1:7
    sdiff(n)=s(1,n)-s(2,n)*2;
end
plot(sdiff)
[sf,f]=fits(sdiff,'lorzs',[100 0 0.1 0 1],[1 1 1 1 0])
plot(sf(n))
formatpars(f(5),2)
title([num2str(files(1)),' minus ',num2str(files(2)),' ; ',ScanCommand])

plot(sf(4))
xlabel('Ql' )
ylabel('counts per 150000 mon' )
text(-0.9,500,['lorz with=',num2str(abs(f(4).pvals(3)),'%1.2f')],'Fontsize',15)
text(-0.9,600,['(100) 0 T minus 1 T'],'Fontsize',15)
print -depsc plots/corr100.eps


%%%(003)
files=[6094 6098];
for n=1:2
    fn=['hdf_data\rita22008n00',num2str(files(n)),'.hdf'];
    data=readhdfdata(fn);
    Magf(n)=double(getfield(data,'magnetic_field'));
    Temp(n)=mean(double(getfield(data,'sample_temperature')));
    ScanCommand=getfield(data,'scancommand');
    s(n,:)=hdf2spec(fn,'Qh',0);
end
for n=1:7
    sdiff(n)=s(1,n)-s(2,n);
end
plot(sdiff)
[sf,f]=fits(sdiff,'lorzs',[100 0 0.1 0 1],[1 1 1 1 1])
plot(sf(n))
formatpars(f(5),2)
title([num2str(files(1)),' minus ',num2str(files(2)),' ; ',ScanCommand])

%%%(200) long
files=[6095 6104]; %6102
for n=1:length(files)
    fn=['hdf_data\rita22008n00',num2str(files(n)),'.hdf'];
    data=readhdfdata(fn);
    Magf(n)=double(getfield(data,'magnetic_field'));
    Temp(n)=mean(double(getfield(data,'sample_temperature')));
    ScanCommand=getfield(data,'scancommand');
    s(n,:)=hdf2spec(fn,'Ql',0);
end
for n=1:7
    sfeld(n)=s(2,n);
    x=getfield(sfeld(n),'x');
    sdiff(n)=cut(s(1,n),[min(x)-0.02,max(x)+0.01])-sfeld(n);
end
plot(sfeld(4),s(1,4))
plot(sdiff(4))
pm=peakm(s(1,4));
pm=peakm(sfeld(4));
plot(sdiff(4))
set(gca,'YLim',[-300,300])

%%%(200) short
files=[6107 6110];
for n=1:length(files)
    fn=['hdf_data\rita22008n00',num2str(files(n)),'.hdf'];
    data=readhdfdata(fn);
    Magf(n)=double(getfield(data,'magnetic_field'));
    Temp(n)=mean(double(getfield(data,'sample_temperature')));
    ScanCommand=getfield(data,'scancommand');
    s(n,:)=hdf2spec(fn,'Ql',0);
end
pm1=peakm(s(1,4));
pm2=peakm(s(2,4));
pm1(2)-pm2(2);
[pm1-pm2]./pm2;
for n=1:7
    s(1,n)=s(1,n)-1*[0.001,0]
    sdiff(n)=s(2,n)-s(1,n);
end
plot(s(2,4),s(1,4))
for n=1:7
plot(sdiff(n),s(1,n)*0.03)
%waitforbuttonpress()
end
h=plot(sdiff(4),s(2,4)*0.03)
xlabel('Ql' )
ylabel('counts' )
legend(h,['0 T minus 1 T   ';'(200) times 0.03'])
print -depsc plots/peak200.eps

%%%(101)
files=[6092 6103]; %6100
for n=1:length(files)
    fn=['hdf_data\rita22008n00',num2str(files(n)),'.hdf'];
    data=readhdfdata(fn);
    Magf(n)=double(getfield(data,'magnetic_field'));
    Temp(n)=mean(double(getfield(data,'sample_temperature')));
    ScanCommand=getfield(data,'scancommand');
    s(n,:)=hdf2spec(fn,'Ql',0);
end
for n=1:7
    sfeld(n)=s(2,n);
    x=getfield(sfeld(n),'x');
    sdiff(n)=cut(s(1,n)/2,[min(x)-0.01,max(x)+0.01])-sfeld(n);
end
plot(sfeld(4),s(1,4)/2)
plot(sdiff(2))
pm=peakm(sfeld(4));
plot(sdiff(4))
set(gca,'YLim',[-300,300])

%%%(004) longitudinal
files=[6079 6084];
for n=1:length(files)
    fn=['hdf_data\rita22008n00',num2str(files(n)),'.hdf'];
    data=readhdfdata(fn);
    Magf(n)=double(getfield(data,'magnetic_field'));
    Temp(n)=mean(double(getfield(data,'sample_temperature')));
    ScanCommand=getfield(data,'scancommand');
    s(n,:)=hdf2spec(fn,'Ql',0);
end
for n=1:7
    sfeld(n)=s(2,n);
    x=getfield(sfeld(n),'x');
    sdiff(n)=cut(s(1,n),[min(x)-0.01,max(x)+0.01])-sfeld(n);
end
plot(sfeld(5),s(1,5))
plot(sdiff(5))
pm=peakm(sfeld(4));
plot(sdiff(4)/pm(1))

%%%(004) transversal low stat
files=[6078 6083];
for n=1:length(files)
    fn=['hdf_data\rita22008n00',num2str(files(n)),'.hdf'];
    data=readhdfdata(fn);
    Magf(n)=double(getfield(data,'magnetic_field'));
    Temp(n)=mean(double(getfield(data,'sample_temperature')));
    ScanCommand=getfield(data,'scancommand');
    s(n,:)=hdf2spec(fn,'Qh',0);
end
for n=1:7
    sfeld(n)=s(2,n);
    x=getfield(sfeld(n),'x');
    sdiff(n)=cut(s(1,n),[min(x)-0.01,max(x)+0.01])-sfeld(n);
end
plot(sfeld(5),s(1,5))
pm=peakm(sfeld(5));
plot(sdiff(5)/pm(1))

%%%(004) transversal high stat
files=[6108 6111];
for n=1:length(files)
    fn=['hdf_data\rita22008n00',num2str(files(n)),'.hdf'];
    data=readhdfdata(fn);
    Magf(n)=double(getfield(data,'magnetic_field'));
    Temp(n)=mean(double(getfield(data,'sample_temperature')));
    ScanCommand=getfield(data,'scancommand');
    s(n,:)=hdf2spec(fn,'Qh',0);
end
for n=1:7
    s(1,n)=s(1,n)+1*[0.001,0];
    sdiff(n)=s(1,n)-s(2,n);
end
plot(s(1,4),s(2,4))
h=plot(sdiff(4),s(2,4)*0.03)
xlabel('Qh' )
ylabel('counts' )
legend(h,['0 T minus 1 T   ';'(004) times 0.03'])
print -depsc plots/peak004.eps

%Fieldscans
files=[6080 6085 6097 6109];
monitor=[10 50 150 150]; %times 1000
for n=1:length(files)
    fn=['hdf_data\rita22008n00',num2str(files(n)),'.hdf'];
    data=readhdfdata(fn);
    Temp(n)=mean(double(getfield(data,'sample_temperature')));
    s(n,:)=hdf2spec(fn,'mf2',1);
    %%binning
    for m=1:7
       % s(n,m)=combine(0.01,s(n,m));
    end
end
plot(s(:,5))
plot(s(1,4))
plot(s(2,2))
plot(s(1,1))
plot(s(4,4),s(2,4)*3.2)
xlabel('Field along c, T' )
ylabel('normalized to monitor' )
text(0.2,1.93,'(200)','Fontsize',15)
text(0.6,1.9,'(004) times 3.2','Fontsize',15)
print -depsc plots/fieldscans.eps
plot(s(4,1),s(2,1)*1.2)
plot(s(1,5),s(2,5))

