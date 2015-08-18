function mv_ulabs

cb=get(findobj('Tag','hmv_CurrentBuffer'),'Userdata');

if ~isempty(cb)
	if iscell(cb), cb=cb{1}; end
end

if ~isempty(cb) & cb~=0

   mv_t=get(findobj('Tag','mv_text_title'),'Userdata');
   mv_x=get(findobj('Tag','mv_text_xlabel'),'Userdata');
   mv_y=get(findobj('Tag','mv_text_ylabel'),'Userdata');

   new_g_label=get(mv_t,'String');
   new_x_label=get(mv_x,'String');
   new_y_label=get(mv_y,'String');

   buffers=get(findobj('Tag','hmv_buffers'),'Userdata');

   buffers(cb).x_label=new_x_label;
   buffers(cb).y_label=new_y_label;
   buffers(cb).g_label=new_g_label;

   set(findobj('Tag','hmv_buffers'),'Userdata',buffers);

end
