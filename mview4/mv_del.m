function mv_del(all)
%
% MATLAB function to delete one or more buffers
% in mview.
%
% DFM 12.12.95
%

%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');
[dummy,nbuffs]=size(buffers);

%----- Unload handles of radio buttons

if nargin==0

   tmv_radio=get(findobj('Tag','tmv_radio'),'Userdata');
   noper_buffs=length(tmv_radio);

   if noper_buffs==0

      mv_msg(' Must select one or more buffers to delete.');
      return

   end

   oper_buffs = zeros([noper_buffs,1]);
   for i=1:noper_buffs

      oper_buffs(i)=str2num(get(tmv_radio(i),'String'));
      if oper_buffs(i) > nbuffs; 

         mv_msg(' Not all selected buffers contain data.');
         return;

      end

   end

elseif nargin==1

%----- Delete all

   set(findobj('Tag','hmv_buffers'),'Userdata',[]);
   hmv_text=get(findobj('Tag','hmv_text'),'Userdata');
   for i=1:nbuffs; set(hmv_text(i),'String',[]); set(hmv_text(i),'ToolTipString',''); end
   mv_rtidy(0)	  
   return

end

%----- Make a string matrix of existing text windows

current_labels=[];
current_tips=[];
hmv_text=get(findobj('Tag','hmv_text'),'Userdata');
for i=1:nbuffs
   current_labels=strvcat(current_labels,get(hmv_text(i),'String'));
   if ~isempty(get(hmv_text(i),'ToolTipString'))
	current_tips = strvcat(current_tips,get(hmv_text(i),'ToolTipString'));
  else
	current_tips = strvcat(current_tips,' ');
  end
   set(hmv_text(i),'String','');
   set(hmv_text(i),'ToolTipString','');
end

%----- Delete specified buffers and labels

oper_buffs=sort(oper_buffs);

buffers(oper_buffs)=[];
current_labels(oper_buffs,:)=[];
current_tips(oper_buffs,:) = [];

set(findobj('Tag','hmv_buffers'),'Userdata',buffers);

%----- Write back remaining labels

[dummy,nbuffs]=size(buffers);

for i=1:nbuffs

   set(hmv_text(i),'String',current_labels(i,:),'ToolTipString', current_tips(i,:));

end

%----- Reset radio buttons

mv_rtidy(0)

return
