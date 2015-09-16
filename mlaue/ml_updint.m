function ml_updint
%  function ml_updint
%
%  MLAUE function to update the Integrated Counts fields
%
% ARW 6.8.07
%
% Last modified:

%===== Find the objects ============================
%----- Addresses to take information
hml_exdat=findobj('Tag','ml_ExtractedImage');   % Contains the (masked) data
hml_bkgtit=findobj('tag','ml_BkgTit');          % Contains the data for the background
%----- Addresses to write outputs -------------------
hml_TotIntBox=[findobj('Tag','ml_TotIntBox') findobj('Tag','ml_eTotIntBox')];
hml_BkgBox=[findobj('Tag','ml_BkgBox') findobj('Tag','ml_eBkgBox')];
hml_IminusBBox=[findobj('Tag','ml_IminusBBox') findobj('Tag','ml_eIminusBBox')];

%===== Get the information
data=get(hml_exdat,'UserData');
bkgdatapts=get(hml_bkgtit,'UserData');

%===== Calculate the integrated intensity.  Errors are sqrt(counts) =====
Itot=sum(sum(data));
Itot=[Itot sqrt(Itot)];

%===== Calculate the Background =========================================
if ~isempty(bkgdatapts)
    bkgdata=[];
    for i=1:length(bkgdatapts(:,1)) 
        tbkg=data(bkgdatapts(i,1),bkgdatapts(i,2));
        if tbkg ~=0, bkgdata=[bkgdata,tbkg]; end
    end
    if ~isempty(bkgdata)
        nbkg=length(bkgdata);
        bkg=[sum(bkgdata)/nbkg sqrt(sum(bkgdata))/nbkg];
    else
        bkg=[0 0];
    end
else
    bkg=[0 0];
end
    
%===== Background subtracted data =======================================
dmb=data(data~=0);
dmb=dmb-bkg(1);
edmb=sqrt(dmb-bkg(2));
ImB=[sum(sum(dmb)) sqrt(sum(sum(edmb.^2)))];

%===== Write the information to the output boxes
set(hml_TotIntBox(1),'String',sprintf('%g',Itot(1)));
set(hml_TotIntBox(2),'String',sprintf('%g',Itot(2)));
set(hml_BkgBox(1),'String',sprintf('%g',bkg(1)));
set(hml_BkgBox(2),'String',sprintf('%g',bkg(2)));
set(hml_IminusBBox(1),'String',sprintf('%g',ImB(1)));
set(hml_IminusBBox(2),'String',sprintf('%g',ImB(2)));
    
