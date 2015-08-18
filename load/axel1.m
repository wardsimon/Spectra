function axel1

hok = findobj('Tag','mfexpr_ok');
hcan= findobj('Tag','mfexpr_cancel');


while (isempty(gco))
    drawnow
 end

   
 while((hok ~= gco) & (hcan ~=gco ))
     drawnow
 end
    
 if(isempty(gco))
     axel1
 end
 
 return