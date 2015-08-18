function mv_rtidy(flag)
%
% Mview function to turn off/on radio buttons 
% in mview.
%
% DFM 21.4.97
%

%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');
[dummy,nbuffs]=size(buffers);

if flag==0

   value=0;

   tmv_radio=get(findobj('Tag','tmv_radio'),'Userdata');
   for i=1:length(tmv_radio)
     set(tmv_radio(i),'Value',value);
   end
   set(findobj('Tag','tmv_radio'),'Userdata',[]);

elseif flag==1

   list=[];
   for i=1:nbuffs
      buf_tag=['mv_buf' num2str(i)];
      h=findobj('Tag',buf_tag);
      set(h,'Value',1);
      list=[list h];
   end
   set(findobj('Tag','tmv_radio'),'Userdata',list);

end
