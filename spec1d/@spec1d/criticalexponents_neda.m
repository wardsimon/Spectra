
%========== Load c-data

feld=[0];
files=cellstr(['0   ']);
for n=1:length(files)
    d= dlmread(['LiErF4_c_',char(files(n)),'t.dat'],'\t',7,0);
    d(:,2)=212.192*d(:,2)/8.314472; %Einheiten umrechnen c/R
    spc(n)=spec1d(d(:,1),d(:,2),0*d(:,2)+0.0001);
end
% m=1:end are different fields. m=1 is zero field
m=1;
s=cut(spc(m),[0,0.7]); %cut away high T tail
plot(s)

%========== HMR

%=== fit to power law

% guess Tc
tcest=0.371;
%tcest=0.369;


% fit both sides simultaneously
sc=cut(s ,[1.002 0.998]*tcest);
sc=cut(sc,[0.7 1.05]*tcest);
[scpf,cpf]=fits(sc,'cpfitfunction',[-1 tcest 0.5 1.92 -1 1.92],[1 0 1 0 1 1]);
figure(1)
clf
plot(s,scpf) % plot on linear scale
line([cpf.pvals([2 2])],[0 2],'color','r')
%print -depsc cp_hmr1

% plot on reduced T, C-B and loglog scale
sllf=-1*(setfield(scpf,'x',abs(getfield(scpf,'x')/cpf.pvals(2)-1))-cpf.pvals(6));
sll=-1*(setfield(s,'x',abs(getfield(s,'x')/cpf.pvals(2)-1))-cpf.pvals(6));
figure(2)
clf
plot(sll,sllf,'loglog')
%axis([3e-3 1 0 1.6])
%print -depsc cp_hmr2


% redo fit for t<0, 0<t<t_crossover and t>t_crossover seperately
[sf(1),f(1)]=fits(cut(s,[0.7 0.998]*tcest),'cpfitfunction',[-5 tcest 0.28 1.92 -3.4 2],[0 0 1 0 1 1]);
sl(1)=-1*(setfield(sf(1),'x',abs(getfield(sf(1),'x')/f(1).pvals(2)-1))-f(1).pvals(6));

[sf(2),f(2)]=fits(cut(s,[1.002 1.05]*tcest),'cpfitfunction',f(1).pvals,[1 0 0 0 0 0]);
sl(2)=-1*(setfield(sf(2),'x',abs(getfield(sf(2),'x')/f(2).pvals(2)-1))-f(2).pvals(6));

[sf(3),f(3)]=fits(cut(s,[1.05 2]*tcest),'cpfitfunction',f(1).pvals,[1 0 1 0 0 0]);
sl(3)=-1*(setfield(sf(3),'x',abs(getfield(sf(3),'x')/f(3).pvals(2)-1))-f(3).pvals(6));

sll=-1*(setfield(s,'x',abs(getfield(s,'x')/f(1).pvals(2)-1))-f(1).pvals(6));
figure(3)
clf
plot(sl,sll,'loglog')

sll(1)=cut(s,[0 0.998]*tcest);
sll(1)=-1*(setfield(sll(1),'x',abs(getfield(sll(1),'x')/f(1).pvals(2)-1))-f(1).pvals(6));
sll(2)=cut(s,[1.002 2]*tcest);
sll(2)=-1*(setfield(sll(2),'x',abs(getfield(sll(2),'x')/f(1).pvals(2)-1))-f(1).pvals(6));
figure(4)
clf
line(1-[0 0.998],-cpfitfunction([0 0.998]*tcest,f(1).pvals)+f(1).pvals(6),'color','k','linestyle','-','linewidth',2)
line([1.002 2]-1,-cpfitfunction([1.002 2]*tcest,f(2).pvals)+f(1).pvals(6),'color','k','linestyle','-','linewidth',2)
line([1.01 2]-1,-cpfitfunction([1.01 2]*tcest,f(3).pvals)+f(1).pvals(6),'color','k','linestyle','-','linewidth',2)
hold on
h=plot(sll([2 1]));set(h(2),'color','r','linewidth',2);set(h(1),'markerfacecolor','none','linewidth',2)
set(gca,'yscale','log','xscale','log')
print -depsc cp_hmr

if 0 %======= logarithmic: cpfitfunction2
tcest=0.37;
sc=cut(s ,[1.0 0.999]*tcest);
sc=cut(sc,[0.7 1.]*tcest);
[scpf,cpf]=fits(sc,'cpfitfunction',[0 tcest -0.5 0 1 -0.1],[0 1 1 0 1 1]);
[scpf2,cpf]=fits(sc,'cpfitfunction2',[0 tcest -1 0 1 -0.1],[0 0 1 1 1 1]);
clf
plot(s,scpf2)
line([cpf.pvals([2 2])],[0 2],'color','r')
print -depsc cp_hmr1

%cpf.pvals(2)=0.38;
sllf=setfield(scpf2,'x',abs(getfield(scpf,'x')/cpf.pvals(2)-1));
sll=setfield(s,'x',abs(getfield(s,'x')/cpf.pvals(2)-1));
clf
plot(sll,sllf,'semilogx')
axis([3e-3 1 0 1.6])
print -depsc cp_hmr2

end % if 0
