function ml_mskpt
%
% function ml_mskpt
%
% MLAUE function to choose points for the background
% ARW 04.08.07
%
% Last modified:

%===== Find the extracted data window ============================
hml_exwin = findobj(0,'Tag','ml_ExtDataWindow'); 
if isempty(hml_exwin) hml_exwin = 0; end
%===== Find the extracted data axes ==============================
hml_exdat = findobj(hml_exwin,'Tag','ml_ExtractedImage');
imclim=get(hml_exdat,'Clim');
data=get(hml_exdat,'UserData');
szdat=fliplr(size(data));
%===== Get background data that has been previously stored ===
hml_bkgrad=findobj('tag','ml_BackGroundRadio');
borders=get(hml_bkgrad,'UserData');if length(borders)==1, borders=[];end
hml_bkgdat=findobj('tag','ml_BkgTit');
xyi=get(hml_bkgdat,'UserData');if length(xyi)==1, xyi=[]; end

%===== Find the status of the button =============================
hml_mskrad=findobj(hml_exwin,'Tag','ml_BackGroundRadio');
butstat=[get(hml_mskrad,'Value')];
if butstat(1,1) == 1
    ml_msg('Click and drag to choose background points, right click to end');
    axes(hml_exdat);hold on
    k=waitforbuttonpress;
    button=1;
    while button ~= 3
        [x1,y1,button]=ginput(1);
        pt1=get(hml_exdat,'CurrentPoint');
        pt1=round(pt1(1,1:2));
        pt1(pt1<=0)=1;
        if pt1(1) > szdat(1), pt1(1)=szdat(1); end;
        if pt1(2) > szdat(2), pt1(2)=szdat(2); end;
        if button ~=3
            mskbox=rbbox;
            pt2=get(hml_exdat,'CurrentPoint');
            pt2=round(pt2(1,1:2));
            pt2(pt2<=0)=1;
            if pt2(1) > szdat(1), pt2(1)=szdat(1); end;
            if pt2(2) > szdat(2), pt2(2)=szdat(2); end;
            p1=min(pt1,pt2);
            msz=abs(pt1-pt2);
%----- Write the borders to a variable for later plotting -----
            borders=[borders;
                     p1(1) p1(2) ;...
                     p1(1)+msz(1) p1(2) ;...
                     p1(1)+msz(1) p1(2)+msz(2) ;...
                     p1(1) p1(2)+msz(2) ;...
                     p1(1) p1(2) ];
            xrange=p1(2):p1(2)+msz(2);
            yrange=p1(1):p1(1)+msz(1);
%----- Write the x and y pixels for the background to an n*2 array -----
            for i=1:length(xrange)
                for j=1:length(yrange)
                     xyi=[xyi;xrange(i) yrange(j)];
                end
            end
        end
%----- Check that there are no doubles -----
       if ~isempty(xyi)
           txyi=[];
           for i=1:length(xyi(:,1))
               colindex=find(xyi(:,1)==xyi(i,1) & xyi(:,2)==xyi(i,2));
               if length(colindex) ==1
                   txyi=[txyi;xyi(i,:)];
               else
                   xyi(i,:)=[-i -i];
               end
           end
       end
       xyi=txyi;
%----- Write the background information to the UserData of the background radio button ---
       set(hml_bkgrad,'UserData',borders);     
       set(hml_bkgdat,'UserData',xyi);     
%----- Update the Integrated Intensities -----------------------
       ml_updint;
%----- Plot data -----
       ml_extplot;
    end
    ml_msg('Background selection ended');
    set(hml_mskrad,'Value',0);
end
    
    
