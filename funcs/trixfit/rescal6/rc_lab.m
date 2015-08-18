function rc_lab(Bragg,phonw,R0,vi,vf,qx,qy)
%
% MATLAB function to put labels on the rescal figure specified by handle h
%
%
% DFM 20.7.95 EF 09.97

fontheight=9;

if ~isempty(findobj('Tag','Rescal: Projections'))
	Fp=get(findobj('tag','Rescal: Projections'),'Position');
	Fp = [ Fp(1:2) 600 600 ];
   	close(findobj('Tag','Rescal: Projections'))
else
	Fp=[50 50 600 600];
end


h=figure('Name','Rescal: Projections',...
                'Tag','Rescal: Projections',...
		'NumberTitle','on',...
                'Units','pixels',...
                'Position',Fp,...
   		'Resize','off',...
		'Visible','on');

% 'Menubar','none',...

set(h,'units','points');
Fp=get(h,'position');

h1=uimenu(h,'Label','Print');
uimenu(h1,'Label','Print Figure','Callback','rc_prt');

%----- Set WYSIWYG characteristics

%unis = get(h,'units');
%ppos = get(h,'paperposition');
%set(h,'units',get(h,'paperunits'));
%pos = get(h,'position');
%pos(3:4) = ppos(3:4);
%set(h,'position',pos);
%set(h,'units',unis);

dgroupy=15;
group0y=Fp(4)-5*fontheight-2*dgroupy;
group1y=group0y-dgroupy;
group2y=group1y-5*fontheight-dgroupy;
group3y=group2y-4*fontheight-dgroupy;
group4y=group3y-0*fontheight-dgroupy;
group5y=group4y-0*fontheight-dgroupy;
group6y=group5y-1*fontheight-dgroupy;
group7y=group6y-0*fontheight-dgroupy;

xshift=-80;
xst=30;

method=get(findobj('tag','hrc_rescal_method'),'Userdata');
if strcmp(method,'rc_cnmat'); method=' Cooper-Nathan''s'; else; method=' Popovici''s'; end
h=text(0,0,'{\bf Method:}','units','points','position',[-40 group0y],'Fontsize',fontheight);
h=text(0,0,method,'units','points','position',[xst group0y],'Fontsize',fontheight);

%----- Spectrometer

h=text(0,0,'{\bf Spectrometer:}',...
           'Units','points',...
           'Fontsize',fontheight,...
           'Position',[-40 group1y]);
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

basepos=[xst group1y];

n=1;
while n<=6;

    line=mat(n,:);
    pos=basepos-[0 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
    n=n+1;

end

n=1;
while n<=6;

    line=mat(n+6,:);
    pos=basepos-[xshift (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
    n=n+1;

end

n=1;
while n<=6;

    line=mat(n+12,:);
    pos=basepos-[2*xshift (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=6;

    line=mat(n+18,:);
    pos=basepos-[3*xshift (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

%----- Create matrix of labels of sample parameters 

h=text(0,0,'{\bf Lattice:}',...
           'Units','points',...
           'Fontsize',fontheight,...
           'Position',[-40 group2y]);

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

basepos=[xst group2y];

n=1;
while n<=5;

    line=mat(n,:);
    pos=basepos-[0 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=5;

    line=mat(n+5,:);
    pos=basepos-[xshift (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=5;

    line=mat(n+10,:);
    pos=basepos-[2*xshift (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=5;

    line=mat(n+15,:);
    pos=basepos-[3*xshift (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

%----- Create matrix of Bragg width parameters 

h=text(0,0,'{\bf Bragg widths:}',...
           'Units','points',...
           'Fontsize',fontheight,...
           'position',[-40 group3y]);

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

basepos=[xst group3y];

n=1;
while n<=5;

    line=mat(n,:);
    pos=basepos-[(n-1)*xshift 0];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

%----- Resolution normalisation labels

h=text(0,0,'{\bf Normalisation:}',...
           'Units','points',...
           'Fontsize',fontheight,...
           'Position',[-40 group4y]);

str=['R0  =' sprintf('%7.3e',R0)];
mat=[str]; 
str=['Vi  =' sprintf('%7.3e',vi)];
mat=[mat; str2mat(str)];
str=['Vf  =' sprintf('%7.3e',vf)];
mat=[mat; str2mat(str)];

%----- Add resolution normalisation labels

basepos=[xst group4y];

n=1;
while n<=3;

    line=mat(n,:);
    pos=basepos-[(n-1)*xshift 0];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

h=text(0,0,'Unit vector projection axes in terms of AS BS CS:',...
'units','points','position',[-40 -42],'Fontsize',fontheight);
h=text(0,0,str2mat(['Qx = ' sprintf('%7.2f%7.2f%7.2f',qx(1),qx(2),qx(3))]),...
'units','points','position',[165 -42],'Fontsize',fontheight); 
h=text(0,0,str2mat(['Qy = ' sprintf('%7.2f%7.2f%7.2f',qy(1),qy(2),qy(3))]),...
'units','points','position',[275 -42],'Fontsize',fontheight); 

%----- Create matrix of labels of dispersion parameters 

h=text(0,0,'{\bf Dispersion:}',...
           'units','points',...
           'Fontsize',fontheight,...
           'position',[-40 group5y]);

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

basepos=[xst group5y];

n=1;
while n<=2;

    line=mat(n,:);
    pos=basepos-[0 (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=2;

    line=mat(n+2,:);
    pos=basepos-[xshift (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=2;

    line=mat(n+4,:);
    pos=basepos-[2*xshift (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

n=1;
while n<=2;

    line=mat(n+6,:);
    pos=basepos-[3*xshift (n-1)*fontheight/1.1];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

%----- Phonon Width label

h=text(0,0,'{\bf Phonon Width:}',...
           'units','points',...
           'Fontsize',fontheight,...
           'position',[-40 group6y]);

str=['        ' sprintf('%7.3f',phonw) '   (meV, 1/Angs or a mixture of both)'];
mat=[str]; 

%----- Add resolution normalisation labels

basepos=[xst group6y];

n=1;
while n<=1;

    line=mat(n,:);
    pos=basepos-[(n-1)*xshift 0];
    h=text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
%    set(h,'units','norm');
    n=n+1;

end

%----- Label Axes

%h=text(0,0,...
%   '{\bf Unit vector projection axes in terms of AS BS CS:}',...
%   'units','points',...
%   'position',[-40 group7y],'Fontsize',fontheight);
%h=text(0,0,str2mat(['Qx = ' sprintf('%7.2f%7.2f%7.2f',qx(1),qx(2),qx(3))]),...
%'units','points','position',[0   group7y-fontheight],'Fontsize',fontheight); 
%h=text(0,0,str2mat(['Qy = ' sprintf('%7.2f%7.2f%7.2f',qy(1),qy(2),qy(3))]),...
%'units','points','position',[100 group7y-fontheight],'Fontsize',fontheight); 

%----- Add the instrumental parameters

method=get(findobj('tag','hrc_rescal_method'),'Userdata');
if strcmp(method,'rc_cnmat')
   return
end

hrc_inst_paras=get(findobj('Tag','hrc_inst_paras'),'Userdata');

pinst(1) =get(hrc_inst_paras(1),'Value');
pinst(2) =str2num(get(hrc_inst_paras(2),'String'));
pinst(3) =str2num(get(hrc_inst_paras(3),'String'));
if pinst(1)  == 1, pinst(3)=pinst(2); end              % Circular source

pinst(4) =get(hrc_inst_paras(4),'Value');
pinst(5) =str2num(get(hrc_inst_paras(5),'String'));
pinst(6) =str2num(get(hrc_inst_paras(6),'String'));
if pinst(4)  == 1, pinst(5)=0; pinst(6)=0; end         % No guide

pinst(7) =get(hrc_inst_paras(7),'Value');
pinst(8) =str2num(get(hrc_inst_paras(8),'String'));
pinst(9) =str2num(get(hrc_inst_paras(9),'String'));
pinst(10)=str2num(get(hrc_inst_paras(10),'String'));
if pinst(7)  == 1, pinst(10)=pinst(9); pinst(9)=pinst(8); end  % Cylindrical sample

pinst(11)=get(hrc_inst_paras(11),'Value');
pinst(12) =str2num(get(hrc_inst_paras(12),'String'));
pinst(13) =str2num(get(hrc_inst_paras(13),'String'));
if pinst(11) == 1, pinst(13)=pinst(12); end            % Circular detector

pinst(1)=pinst(1)-1;   pinst(4)=pinst(4)-1;   pinst(7)=pinst(7)-1;   pinst(11)=pinst(11)-1;

for i=14:length(hrc_inst_paras)
        if length(str2num(get(hrc_inst_paras(i),'String'))) ~= 1
	   disp([ ' Cannot get parameter n.' num2str(i) ' = ' get(hrc_inst_paras(i),'String') ' from Instr window'])
	 else 
           pinst(i)=str2num(get(hrc_inst_paras(i),'String'));
	 end
end

%----- Positions for text

dgroupy=12;
i1y=group7y-20*fontheight-dgroupy+50; % added +50 for Matlab 6.5
i2y=i1y-0*fontheight-dgroupy;
i3y=i2y-0*fontheight-dgroupy;
i4y=i3y-0*fontheight-dgroupy;
i5y=i4y-0*fontheight-dgroupy;
i6y=i5y-0*fontheight-dgroupy;
i7y=i6y-0*fontheight-dgroupy;
i8y=i7y-0*fontheight-dgroupy;

%----- Source

h=text(0,0,'{\bf Source:}',...
           'Units','points',...
           'Fontsize',fontheight,...
           'Position',[-40 i1y]);

if pinst(1)==1

   str=['Diameter =' sprintf('%7.3f',pinst(2))];
   mat=[str]; 
   str=['                 '];
   mat=[mat; str2mat(str)];


else

   str=['Width    =' sprintf('%7.3f',pinst(2))];
   mat=[str]; 
   str=['Height   =' sprintf('%7.3f',pinst(3))];
   mat=[mat; str2mat(str)];

end

basepos=[xst i1y];

n=1;
while n<=2;

    line=mat(n,:);
    pos=basepos-[(n-1)*xshift 0];
    text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
    n=n+1;

end

%----- Guide

h=text(0,0,'{\bf Guide:}',...
           'Units','points',...
           'Fontsize',fontheight,...
           'Position',[-40 i2y]);

if pinst(4)==1

   str=['Hor. Div.=' sprintf('%7.3f',pinst(5))];
   mat=[str]; 
   str=['Ver. Div.=' sprintf('%7.3f',pinst(6))];
   mat=[mat; str2mat(str)];


else

   str=['No guide         '];
   mat=[str]; 
   str=['                 '];
   mat=[mat; str2mat(str)];

end

basepos=[xst i2y];

n=1;
while n<=2;

    line=mat(n,:);
    pos=basepos-[(n-1)*xshift 0];
    text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
    n=n+1;

end

%----- Sample

h=text(0,0,'{\bf Sample:}',...
           'Units','points',...
           'Fontsize',fontheight,...
           'Position',[-40 i3y]);

if pinst(7)==1

   str=['Depth    =' sprintf('%7.3f',pinst(8))];
   mat=[str]; 
   str=['Width    =' sprintf('%7.3f',pinst(9))];
   mat=[mat; str2mat(str)];
   str=['Height   =' sprintf('%7.3f',pinst(10))];
   mat=[mat; str2mat(str)];


else

   str=['Diameter =' sprintf('%7.3f',pinst(8))];
   mat=[str]; 
   str=['Height   =' sprintf('%7.3f',pinst(10))];
   mat=[mat; str2mat(str)];
   str=['                 '];
   mat=[mat; str2mat(str)];

 
end

basepos=[xst i3y];

n=1;
while n<=3;

    line=mat(n,:);
    pos=basepos-[(n-1)*xshift 0];
    text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
    n=n+1;

end

%----- Source

h=text(0,0,'{\bf Detector:}',...
           'Units','points',...
           'Fontsize',fontheight,...
           'Position',[-40 i4y]);

if pinst(11)==0

   str=['Diameter =' sprintf('%7.3f',pinst(12))];
   mat=[str]; 
   str=['                 '];
   mat=[mat; str2mat(str)];


else

   str=['Width    =' sprintf('%7.3f',pinst(12))];
   mat=[str]; 
   str=['Height   =' sprintf('%7.3f',pinst(13))];
   mat=[mat; str2mat(str)];

end

basepos=[xst i4y];

n=1;
while n<=2;

    line=mat(n,:);
    pos=basepos-[(n-1)*xshift 0];
    text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
    n=n+1;

end

%----- Monochromator

h=text(0,0,'{\bf Monochromator:}',...
           'Units','points',...
           'Fontsize',fontheight,...
           'Position',[-40 i5y]);

str=['Thickness=' sprintf('%7.3f',pinst(14))];
mat=[str]; 
str=['Width    =' sprintf('%7.3f',pinst(15))];
mat=[mat; str2mat(str)];
str=['Height   =' sprintf('%7.3f',pinst(16))];
mat=[mat; str2mat(str)];

basepos=[xst i5y];

n=1;
while n<=3;

    line=mat(n,:);
    pos=basepos-[(n-1)*xshift 0];
    text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
    n=n+1;

end

%----- Monochromator

h=text(0,0,'{\bf Analyser:}',...
           'Units','points',...
           'Fontsize',fontheight,...
           'Position',[-40 i6y]);

str=['Thickness=' sprintf('%7.3f',pinst(17))];
mat=[str]; 
str=['Width    =' sprintf('%7.3f',pinst(18))];
mat=[mat; str2mat(str)];
str=['Height   =' sprintf('%7.3f',pinst(19))];
mat=[mat; str2mat(str)];

basepos=[xst i6y];

n=1;
while n<=3;

    line=mat(n,:);
    pos=basepos-[(n-1)*xshift 0];
    text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
    n=n+1;

end

%----- Distances

h=text(0,0,'{\bf Distances:}',...
           'Units','points',...
           'Fontsize',fontheight,...
           'Position',[-40 i7y]);

str=['L0       =' sprintf('%10.3f',pinst(20))];
mat=[str]; 
str=['L1       =' sprintf('%10.3f',pinst(21))];
mat=[mat; str2mat(str)];
str=['L2       =' sprintf('%10.3f',pinst(22))];
mat=[mat; str2mat(str)];
str=['L3       =' sprintf('%10.3f',pinst(23))];
mat=[mat; str2mat(str)];

basepos=[xst i7y];

n=1;
while n<=4;

    line=mat(n,:);
    pos=basepos-[(n-1)*xshift 0];
    text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
    n=n+1;

end

%----- Distances

h=text(0,0,'{\bf Distances:}',...
           'Units','points',...
           'Fontsize',fontheight,...
           'Position',[-40 i8y]);

str=['ROMH     =' sprintf('%7.3f',pinst(24))];
mat=[str]; 
str=['ROMV     =' sprintf('%7.3f',pinst(25))];
mat=[mat; str2mat(str)];
str=['ROAH     =' sprintf('%7.3f',pinst(26))];
mat=[mat; str2mat(str)];
str=['ROAV     =' sprintf('%7.3f',pinst(27))];
mat=[mat; str2mat(str)];

basepos=[xst i8y];

n=1;
while n<=4;

    line=mat(n,:);
    pos=basepos-[(n-1)*xshift 0];
    text(0,0,line,'units','points','position',pos,'Fontsize',fontheight);
    n=n+1;

end






