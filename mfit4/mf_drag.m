function mf_drag
% Drag text line, or group of text lines, around

hh=gco;
if strcmp(get(hh,'Type'),'text')
   h=get(hh,'userdata');

   if length(h)>1                         % Line liked to group
      p=get(gca,'CurrentPoint');          % Current point
      p=p(1,:);
      delta=p-get(hh,'position');         % Offset of current object
      for i=1:length(h)
         pos=get(h(i),'position');        
         set(h(i),'position',pos+delta);  % Add offset to all lines
      end
   else                                   % Single line
      p=get(gca,'CurrentPoint');
      p=p(1,:);
      set(hh,'position',p);
   end      
end

