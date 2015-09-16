function mf_print

[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;
figure(hmf_data);
oldu = get(gcf,'units');
oldpu = get(gcf,'paperunits');
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
p=get(gcf,'position');
width=p(3); height=p(4);

paper=get(gcf,'papersize');
ppos=[(paper(1)-width)/2 (paper(2)-height)/2 width height];
set(gcf,'paperposition',ppos);

if strcmp(computer,'PCWIN')
   print -dwin -v -painters
else
   print -painters
end
datafile = get(findobj('tag','mf_DataFile'),'string');
to = printopt;
fprintf(1,'Print figure %i (MFit Data:%s) with : %s\n',hmf_data,datafile,to);
set(gcf,'units',oldu);
set(gcf,'paperunits',oldpu);
