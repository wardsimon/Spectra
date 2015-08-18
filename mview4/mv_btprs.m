function mv_btprs
%
% MFIT function  mv_btprs
% Process button press on data window
% MZ 14.3.95
%

objtype=get(gco,'type');

if strcmp(objtype,'text')==1
	mv_text(gco);
elseif strcmp(objtype,'axes')==1
	mv_rbbox('go');
end
