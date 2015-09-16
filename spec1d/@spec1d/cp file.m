d= dlmread(['LiErF42','.dat'],'\t');
spc(1)=spec1d(d(:,1),d(:,2),0*d(:,2)+0.000001);
plot(spc(1))
s1=cut(spc(1),[0,0.7]); %cut away high T tail
plot(s1)
tcest=0.3689
%tcest0.3674
tcest=0.3675
tcest=0.3685
tcest=0.37

sc1=cut(s1 ,[0.871 1.03]*tcest);
plot(sc1)
%sc2=cut(sc1 ,[1.004 0.9959]*tcest);
sc2=cut(sc1 ,[1.005 0.993]*tcest);
plot(sc2)
%[scpf,cpf]=fits(sc2,'cpfitfunction',[-0.001 tcest 0.47 1.92 -0.001 0.00048],[1 0 1 0 1 0]);
%[scpf,cpf]=fits(sc2,'cpfitfunction',[-0.001 tcest 0.47 0.000522 -0.001 0],[1 0 1 1 1 1]);
[scpf,cpf]=fits(sc2,'cpfitfunction',[-0.001 tcest 0.29 0.000494 -0.001 0],[1 0 0 1 1 0]);
figure(1)
clf
plot(s1,sc2)
t=0.1:0.0002:1;
line(t,cpfitfunction(t,cpf.pvals),'color','r','linewidth',2)
line([cpf.pvals([2 2])],[0 0.00045],'color','k')

scold=cut(s1,[0 cpf.pvals(2)]);
shot=cut(s1,[cpf.pvals(2) 2]);

sllc=-1*setfield(scold,'x',abs(getfield(scold,'x')/cpf.pvals(2)-1))+cpf.pvals(4)+cpf.pvals(6);
sllh=-1*setfield(shot,'x',abs(getfield(shot,'x')/cpf.pvals(2)-1))+cpf.pvals(4);
figure(2)
clf
plot(sllc,sllh,'loglog')
line(abs(t/cpf.pvals(2)-1),cpfitfunction(t,cpf.pvals.*[-1 1 1 0 -1 0]'));

clf
plot(s1,scpf)


%===== play with confoluting the power-law function to see how it rounds
t=0.35:0.0002:0.4;
c=cpfitfunction(t,cpf.pvals);
d=0.000235;
nd=-10:10;
ed=exp(-0.5*(nd/4).^2);
c2=0*c;
for n=1:length(nd)
  c2=c2+ed(n)*cpfitfunction(t+nd(n)*d,cpf.pvals);
end
c2=c2/sum(ed);
sf=spec1d(t,c,1e-6);
sf2=spec1d(t,c2,1e-6);
clf
h=plot(sf,smooth(sf(1),0.001),sf2);
set(h,'marker','.','markersize',1,'linestyle','-')


%===== fit to a smoothed powerlaw function. p(7) is the smoothing width

[scpf,cpf]=fits(sc1,'cpfitfunction3',[-0.001 tcest 0.29 0.000494 -0.001 0 0.0001],[1 1 1 0 1 0 0]);
figure(1)
clf
plot(sc1,scpf)
t=0.3:0.0001:4;
line(t,cpfitfunction3(t,cpf.pvals),'color','r','linewidth',2)
line([cpf.pvals([2 2])],[0 0.00045],'color','k')

scold=cut(sc1,[0 cpf.pvals(2)]);
shot=cut(sc1,[cpf.pvals(2) 2]);

sllc=-1*setfield(scold,'x',abs(getfield(scold,'x')/cpf.pvals(2)-1))+cpf.pvals(4)+cpf.pvals(6);
sllh=-1*setfield(shot,'x',abs(getfield(shot,'x')/cpf.pvals(2)-1))+cpf.pvals(4);
figure(2)
clf
plot(sllc,sllh,'loglog')
line(abs(t/cpf.pvals(2)-1),cpfitfunction3(t,cpf.pvals.*[-1 1 1 0 -1 0 0.1]'));
line(abs(t/cpf.pvals(2)-1),cpfitfunction(t,cpf.pvals.*[-1 1 1 0 -1 0 1]'));
