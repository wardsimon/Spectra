function axel

hok = findobj('Tag','mfload_ok');
hcan= findobj('Tag','mfload_cancel');



while (isempty(gco))
    drawnow
 end

   
 while((hok ~= gco) & (hcan ~=gco ))
     drawnow
 end
    
 if(isempty(gco))
     axel
 end
 
 return