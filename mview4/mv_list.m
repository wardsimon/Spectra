function mv_list
%
% MATLAB function to get data into mview
%
% DFM 26.2.95
%

%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');
[dummy,nbuffs]=size(buffers);

%----- Perform operation on buffers: list

%----- First print out ranges 

fprintf('\n-----------------------------------------------------------------------------')
fprintf('\n                            Buffer Details             ')
fprintf('\n  #  Points  X min      X max     X range    Y min      Y max      Y range\n')
fprintf('-----------------------------------------------------------------------------\n')

for i=1:nbuffs
  
      x_b=buffers(i).xobs;
      y_b=buffers(i).yobs;
      title=buffers(i).g_label;
      npts=length(x_b); 
      minx=min(x_b);
      maxx=max(x_b);
      rangex=max(x_b)-min(x_b);    
      miny=min(y_b);
      maxy=max(y_b);
      rangey=max(y_b)-min(y_b);
      fprintf('%3i  %4i  %6.2e  %6.2e  %6.2e  %6.2e  %6.2e  %6.2e \n',...
               i,npts,minx,maxx,rangex,miny,maxy,rangey)

end

fprintf('\n\n')
fprintf('\n------------------------------------------------------')
fprintf('\n             Statistical Analysis                 ')
fprintf('\n  #   Y Sum      Y Mean     1st moment  2nd moment\n')
fprintf('------------------------------------------------------\n')

for i=1:nbuffs
   
      x_b=buffers(i).xobs;
      y_b=buffers(i).yobs;
      npts=length(x_b);
      title=buffers(i).g_label;
      minx=min(x_b);
      maxx=max(x_b);
      stepsz=(max(x_b)-min(x_b))/(npts-1);    
      maxy=max(y_b);
      sumy=sum(y_b);
      aver_y=sumy/npts;
      fmom=sum(x_b.*y_b)/sum(y_b);
      smom=sqrt(sum(x_b.^2.*y_b)/sum(y_b)-(sum(x_b.*y_b)/sum(y_b))^2);

      fprintf('%3i  %6.2e   %6.2e   %6.2e   %6.2e\n',...
               i,sumy,aver_y,fmom,smom)

end

return

