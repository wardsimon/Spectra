function rc_lab(Bragg,phonw,R0,vi,vf,qx,qy)
%
% MATLAB function to put labels on the rescal figure specified by handle h
%
%
% DFM 20.7.95

if ~isempty(findobj('Tag','Rescal: Projections'))
   close(findobj('Tag','Rescal: Projections'))
end

h=figure('Name','Rescal: Projections',...
                'Tag','Rescal: Projections',...
		'NumberTitle','off',...
		'Menubar','None',...
		'Visible','on','position',[300 100 50 50 ]);

h1=uimenu(h,'Label','Print');
uimenu(h1,'Label','Print Figure','Callback','print');

%----- Set WYSIWYG characteristics

unis = get(h,'units');
ppos = get(h,'paperposition');
set(h,'units',get(h,'paperunits'));
pos = get(h,'position');
pos(3:4) = ppos(3:4);
set(h,'position',pos);
set(h,'units',unis);
drawnow

group1y=350;
group2y=275;
group3y=210;
group4y=180;
group5y=150;
group6y=105;
dgroupy=15;

method=get(findobj('tag','hrc_rescal_method'),'Userdata');
if strcmp(method,'rc_cnmat'); method=' Cooper-Nathan''s'; else; method=' Popovici''s'; end
h=text(0,0,'Method:','units','points','position',[-40 375]);
h=text(0,0,method,'units','points','position',[20 374],'Fontsize',10);

%----- Spectrometer

h=text(0,0,'Spectrometer:','units','points','position',[-40 group1y]);
set(gca,'visible','off')

%----- Get parameter values

p=rc_savp('respars');

%----- Create matrix of labels of instrument parameters 

str=['DM  =' sprintf('%7.3f',p(1))];
mat=[str]; 
str=['ETAM=' sprintf('%7.3f',p(3))];
mat=[mat; str2mat(str)];
str=['SM  =' sprintf('%7.3f',p(6))];
mat=[mat; str2mat(str)];
str=['KFIX=' sprintf('%7.3f',p(9))];
mat=[mat; str2mat(str)];
str=['ALF1=' sprintf('%7.3f',p(11))];
mat=[mat; str2mat(str)];
str=['BET1=' sprintf('%7.3f',p(15))];
mat=[mat; str2mat(str)];

str=['DA  =' sprintf('%7.3f',p(2))];
mat=[mat; str2mat(str)]; 
str=['ETAA=' sprintf('%7.3f',p(4))];
mat=[mat; str2mat(str)];
str=['SS  =' sprintf('%7.3f',p(7))];
mat=[mat; str2mat(str)];
str=['FX  =' sprintf('%7.3f',p(10))];
mat=[mat; str2mat(str)];
str=['ALF2=' sprintf('%7.3f',p(12))];
mat=[mat; str2mat(str)];
str=['BET2=' sprintf('%7.3f',p(16))];
mat=[mat; str2mat(str)];

str=['            '];
mat=[mat; str2mat(str)]; 
str=['ETAS=' sprintf('%7.3f',p(5))];
mat=[mat; str2mat(str)];
str=['SA  =' sprintf('%7.3f',p(8))];
mat=[mat; str2mat(str)];
str=['            '];
mat=[mat; str2mat(str)];
str=['ALF3=' sprintf('%7.3f',p(13))];
mat=[mat; str2mat(str)];
str=['BET3=' sprintf('%7.3f',p(17))];
mat=[mat; str2mat(str)];

str=['            '];
mat=[mat; str2mat(str)]; 
str=['            '];
mat=[mat; str2mat(str)];
str=['            '];
mat=[mat; str2mat(str)];
str=['            '];
mat=[mat; str2mat(str)];
str=['ALF4=' sprintf('%7.3f',p(14))];
mat=[mat; str2mat(str)];
str=['BET4=' sprintf('%7.3f',p(18))];
mat=[mat; str2mat(str)];

%----- Add instrument parameter labels

basepos=[10 group1y-dgroupy];
fontheight=10;

n=1;
while n<=6;

    line=mat(n,:);
    pos=basepos-[0 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=6;

    line=mat(n+6,:);
    pos=basepos-[-100 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=6;

    line=mat(n+12,:);
    pos=basepos-[-200 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=6;

    line=mat(n+18,:);
    pos=basepos-[-300 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

%----- Create matrix of labels of sample parameters 

h=text(0,0,'Lattice:','units','points','position',[-40 group2y]);

str=['AS  =' sprintf('%7.3f',p(19))];
mat=[str]; 
str=['AA  =' sprintf('%7.3f',p(22))];
mat=[mat; str2mat(str)];
str=['AX  =' sprintf('%7.3f',p(25))];
mat=[mat; str2mat(str)];
str=['BX  =' sprintf('%7.3f',p(28))];
mat=[mat; str2mat(str)];
str=['QH  =' sprintf('%7.3f',p(31))];
mat=[mat; str2mat(str)];

str=['BS  =' sprintf('%7.3f',p(20))];
mat=[mat; str2mat(str)]; 
str=['BB  =' sprintf('%7.3f',p(23))];
mat=[mat; str2mat(str)];
str=['AY  =' sprintf('%7.3f',p(26))];
mat=[mat; str2mat(str)];
str=['BY  =' sprintf('%7.3f',p(29))];
mat=[mat; str2mat(str)];
str=['QK  =' sprintf('%7.3f',p(32))];
mat=[mat; str2mat(str)];

str=['CS  =' sprintf('%7.3f',p(21))];
mat=[mat; str2mat(str)]; 
str=['CC  =' sprintf('%7.3f',p(24))];
mat=[mat; str2mat(str)];
str=['AZ  =' sprintf('%7.3f',p(27))];
mat=[mat; str2mat(str)];
str=['BZ  =' sprintf('%7.3f',p(30))];
mat=[mat; str2mat(str)];
str=['QL  =' sprintf('%7.3f',p(33))];
mat=[mat; str2mat(str)];

str=['            '];
mat=[mat; str2mat(str)]; 
str=['            '];
mat=[mat; str2mat(str)];
str=['            '];
mat=[mat; str2mat(str)];
str=['            '];
mat=[mat; str2mat(str)];
str=['EN  =' sprintf('%7.3f',p(34))];
mat=[mat; str2mat(str)];

%----- Add sample parameter labels

basepos=[10 group2y-dgroupy];
fontheight=10;

n=1;
while n<=5;

    line=mat(n,:);
    pos=basepos-[0 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=5;

    line=mat(n+5,:);
    pos=basepos-[-100 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=5;

    line=mat(n+10,:);
    pos=basepos-[-200 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=5;

    line=mat(n+15,:);
    pos=basepos-[-300 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

%----- Create matrix of Bragg width parameters 

h=text(0,0,'Bragg widths:','units','points','position',[-40 group3y]);

str=['Qx  =' sprintf('%7.3f',Bragg(1))];
mat=[str]; 
str=['Qy  =' sprintf('%7.3f',Bragg(2))];
mat=[mat; str2mat(str)];
str=['Qz  =' sprintf('%7.3f',Bragg(3))];
mat=[mat; str2mat(str)];
str=['V   =' sprintf('%7.3f',Bragg(4))];
mat=[mat; str2mat(str)];
str=['DEE =' sprintf('%7.3f',Bragg(5))];
mat=[mat; str2mat(str)];

%----- Add Bragg width parameter labels

basepos=[10 group3y-dgroupy];
fontheight=10;

n=1;
while n<=5;

    line=mat(n,:);
    pos=basepos-[-(n-1)*100 0];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

%----- Resolution normalisation labels

h=text(0,0,'Resolution normalisation:','units','points','position',[-40 group4y]);
str=['R0  =' sprintf('%7.3e',R0)];
mat=[str]; 
str=['Vi  =' sprintf('%7.3e',vi)];
mat=[mat; str2mat(str)];
str=['Vf  =' sprintf('%7.3e',vf)];
mat=[mat; str2mat(str)];

%----- Add resolution normalisation labels

basepos=[10 group4y-dgroupy];
fontheight=10;

n=1;
while n<=3;

    line=mat(n,:);
    pos=basepos-[-(n-1)*100 0];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

%----- Create matrix of labels of dispersion parameters 

h=text(0,0,'Dispersion surface:','units','points','position',[-40 group5y]);

str=['DH  =' sprintf('%7.3f',p(35))];
mat=[str]; 
str=['GH  =' sprintf('%7.3f',p(39))];
mat=[mat; str2mat(str)];
str=['DK  =' sprintf('%7.3f',p(36))];
mat=[mat; str2mat(str)];
str=['GK  =' sprintf('%7.3f',p(40))];
mat=[mat; str2mat(str)];

str=['DL  =' sprintf('%7.3f',p(37))];
mat=[mat; str2mat(str)];
str=['GL  =' sprintf('%7.3f',p(41))];
mat=[mat; str2mat(str)]; 
str=['DE  =' sprintf('%7.3f',p(38))];
mat=[mat; str2mat(str)];
str=['GMOD=' sprintf('%7.3f',p(42))];
mat=[mat; str2mat(str)];

%----- Add Dispersion labels

basepos=[10 group5y-dgroupy];
fontheight=10;

n=1;
while n<=2;

    line=mat(n,:);
    pos=basepos-[0 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=2;

    line=mat(n+2,:);
    pos=basepos-[-100 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=2;

    line=mat(n+4,:);
    pos=basepos-[-200 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=2;

    line=mat(n+6,:);
    pos=basepos-[-300 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

%----- Phonon Width label

h=text(0,0,'Phonon Width:','units','points','position',[-40 group6y]);


str=['        ' sprintf('%7.3f',phonw) '   (meV, 1/Angs or a mixture of both)'];
mat=[str]; 

%----- Add resolution normalisation labels

basepos=[10 group6y-dgroupy];
fontheight=10;

n=1;
while n<=1;

    line=mat(n,:);
    pos=basepos-[-(n-1)*100 0];
    h=text(0,0,line,'units','points','position',[30 105],'Fontsize',fontheight);
   set(h,'units','norm');
    n=n+1;

end

h=text(0,0,'Projection axes in terms of AS BS CS:','units','points','position',[-40 -40]);
h=text(0,0,str2mat(['Qx = ' sprintf('%7.2f%7.2f%7.2f',qx(1),qx(2),qx(3))]),...
'units','points','position',[200 -40],'Fontsize',fontheight); 
h=text(0,0,str2mat(['Qy = ' sprintf('%7.2f%7.2f%7.2f',qy(1),qy(2),qy(3))]),...
'units','points','position',[300 -40],'Fontsize',fontheight); 





