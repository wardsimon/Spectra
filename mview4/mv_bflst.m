function [list]=mv_bflst(buftag);
%
% function mv_bflst(h)
%
% Update the list of selected buffers. The 'userdata' of the
% object with tag 'tmv_radio' is a list of handles of radio
% buttons (whose 'userdata' are the buffer data) 
% selected buffers in the order that they were selected.
% The parameter h is the handle of the radiobutton pushed.


list=get(findobj('Tag','tmv_radio'),'Userdata');

if nargin==1

	h=findobj('Tag',buftag);
	if ~isempty(list)
		p=find(list==h);
	else
		p=[];
	end
	if isempty(p)
		list=[list h];
	else
		list(p)=[];
	end

	set(findobj('Tag','tmv_radio'),'Userdata',list);

end

