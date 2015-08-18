function ml_mskpt
%
% function ml_mskpt
%
% MLAUE function to mask points in extracted data
% ARW 03.08.07
%
% Last modified:

%===== Find the extracted data window ============================
hml_exwin = findobj(0,'Tag','ml_ExtDataWindow'); 
if isempty(hml_exwin) hml_exwin = 0; end
%===== Find the extracted data axes ==============================
tagname='ml_ExtractedImage';
hml_exdat = findobj(hml_exwin,'Tag',tagname);
imclim=get(hml_exdat,'Clim');
data=get(hml_exdat,'UserData');
szdat=fliplr(size(data));

%===== Find the status of the button =============================
hml_mskrad=findobj(hml_exwin,'Tag','ml_MaskRadio');
butstat=[get(hml_mskrad,'Value')];
if butstat(1,1) == 1
    ml_msg('Click and drag to mask points, right click to end');
    axes(hml_exdat);
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
            data(p1(2):p1(2)+msz(2),p1(1):p1(1)+msz(1))=0;
            set(hml_exdat,'UserData',data);
%----- Calculate Integrated Intensities ------
            ml_updint;
%----- Update Plot ---------------------------
            ml_extplot;
        end
    end
    ml_msg('Mask selection ended');
%----- Set RadioButton to 'Off' --------------
    set(hml_mskrad,'Value',0);
end
    
    
