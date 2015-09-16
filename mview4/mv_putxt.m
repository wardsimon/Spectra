function mv_puttxt(cmd)

[hmv_ctrl, hmv_data]=mv_figs;
if hmv_data==0 
	return
end

return  % force return

figure(hmv_data);
mv_msg('Click on text insertion point');
h=gtext(' ');

set(h,'units','points');
basepos=get(h,'Position');

h=get(hmv_pars,'userdata');
s1=sprintf('%10s %10s %10s',...
				get(h(1,3),'String'),...
				get(h(1,1),'String'),...
				get(h(1,2),'String'));
for j=2:size(h,1)
	s2=sprintf('%10s %10s %10s',...
				get(h(j,3),'String'),...
				get(h(j,1),'String'),...
				get(h(j,2),'String'));
	s1=str2mat(s1,s2);
end

[lines cols]=size(s1);

n=1;
while n<=lines
    line=s1(n,:);
    pos=basepos-[0 (n-1)*10];
    h=text(0,0,line,...
    		'units','points',...
    		'Position',pos,...
    		'units','norm',...
    		'Fontsize',10,...
    		'Tag','mv_text_user');
    n=n+1;
end

mv_msg(' ');
