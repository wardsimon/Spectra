function [y, name, pnames, pin,separate]=multifunc(x,p, flag)
% multifunc : Multi Function handler window
% function [y, name, pnames, pin]=multifunc((x,p, flag)
%
% MFIT Multi Function handler window
% * Select sub-function number to define.
% * You can use any expression in 'fn=' fields, including user variables ('user dummy')
%    and any evaluable string in 'y=' field.
% * The constrain field enables to enter any matlab expression, such as 'p(1) = 0;'
% * User can choose to use an external function (pre-compiled on disk) or an internal
%    evaluated function (slower).
% * User can enable separate plotting for each sub-function by clicking on 'fn' fields
%    and re-building function.
%
% Batch direct command:
% flag can be a command for function building
% 'n=..'  changes the number of sub-functions (for example 'n=3')
% 'fn=..' changes the sub-function n ('f2=gauss')
% 'y='    set the final expression ('y=f1.*f2+f3')
% 'c=..'  set the optional constrain field
% 'e'     use 'eval' mode (slow but efficient) instead of external file mode
% 'getbatch' returns in 'y' the equivalent mf_batch command to build multifunc
%
% a mfit batch example is :
% exec multifunc(x,p,'n=3');
% exec multifunc(x,p,'f1=gauss');
% exec multifunc(x,p,'f2=lorz ID');       % will add 'ID' comment in parameter names
% exec multifunc(x,p,'f3=user MyParam');  % create 1 user param
% exec multifunc(x,p,'y=f1+f2+f3');
%
% if function name is 'user', a user parameter is defined, name follows
% or is prompted if not given. ex : f1=user temperature

% Author:  EF <farhi@ill.fr> 07.97 KD 05.99
% Description:  Multi Function handler window

%=========== Make window... =======================




[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

if (nargin < 3)
  flag = 0;
end

name = 'multifunc';
y=[];
pin = [];
pnames = '';
separate={};
NL = sprintf('\n');

hmf_multif=findobj('Tag','mf_mfunc');

% handle window as long as needed : user nb_func change, no window,...

if isempty(hmf_multif)
   mf_msg('Building multifunc window');
   nfunc = 2;

%------- Create figure window -------------------------------------
   hmf_multif=figure('Position',[400 400 270 (4+nfunc*26+28+28+28+28+28)],...
    'Tag','mf_mfunc',...
      'MenuBar','none',...
      'Name','MFIT: Multi Functions',...
      'Color',get(0,'DefaultUicontrolBackgroundColor'), ...
      'Resize','off',...
          'Visible','on');

%-------- Nb funcs objects : text, label, edit -----------------------

   uicontrol(hmf_multif,...
      'Style','text',...
    'Tag','mff_text1',...
      'Position',[50 (4+nfunc*26+28+28+28+28) 130 16],...
      'ForegroundColor',[1 1 1],...
    'ToolTipString',[ 'Enter here the name of the' NL 'functions to be assembled' ], ...
      'String','Function names');


   uicontrol(hmf_multif,...     % edit label
      'Style','text',...
      'ForegroundColor',[1 1 1],...
      'Position',[10 32 150 16],...
    'ToolTipString',[ 'Modify the number of sub-functions to use' NL 'and click on "Make functions"' ], ...
      'String','Number of functions');

   uicontrol(hmf_multif,...                       % number of functions edit field
      'Tag','mf_mfn',...
      'Position',[160 32 50 18],...
    'Style','edit',...
      'String','2',...
      'BackgroundColor',[1 1 1],...
      'ForegroundColor',[0 0 0],...
    'Visible','on',...
    'Clipping','off');

   uicontrol(hmf_multif,...     % expression label
      'Style','text',...
      'ForegroundColor',[1 1 1],...
      'Position',[10 60 70 16],...
    'ToolTipString',[ 'This is the expresion to be evaluated' NL 'in order to compute the multi function value.' NL ' You may use f1-fn values' ], ...
      'String','Final y=');

   uicontrol(hmf_multif,...                       % expression field
      'Tag','mf_mfexpr',...
      'Position',[80 60 170 20],...
    'Style','edit',...
      'String','f1+f2',...
      'BackgroundColor',[1 1 1],...
      'ForegroundColor',[0 0 0],...
    'Visible','on',...
    'Clipping','off');

   uicontrol(hmf_multif,...     % expression label
      'Style','text',...
      'ForegroundColor',[1 1 1],...
      'Position',[10 88 70 16],...
    'ToolTipString',[ 'You may enter any command to be executed when evaluating the final function' NL 'This may e.g. be constrains on parameters such as "p(1) = p(3)"' ], ...
      'String','Constrain');

   uicontrol(hmf_multif,...                       % expression field
      'Tag','mf_mfconst',...
        'Position',[80 88 170 20],...
    'Style','edit',...
      'String','',...
      'BackgroundColor',[1 1 1],...
      'ForegroundColor',[0 0 0],...
    'Visible','on',...
    'Clipping','off');

%------- ok and cancel buttons -------------------------------------------
   uicontrol(hmf_multif,...
      'Style','PushButton',...
      'String','Make Function',...
    'ToolTipString',[ 'Click here to create the multi-function' ], ...
      'Tag','mfmf_ok',...
      'Position',[4 4 90 20]);
   uicontrol(hmf_multif,...
      'Style','PushButton',...
      'String','Cancel',...
      'Tag','mfmf_cancel',...
      'CallBack','multifunc([],[],''hide'');', ...
      'Position',[100 4 50 20]);

   uicontrol(hmf_multif,...
    'Style','checkbox',....
    'Tag','mfmf_extf',...
    'Value',0,...
        'ToolTipString',[ 'If checked, will save the function into a file' NL 'else work in memory' ], ...
    'String','Use Extn. File',...
    'Position',[155 4 110 20],...
        'visible','on');




% eval str save field ------------------------------------------
   uicontrol(hmf_multif,...
      'Style','text',...
      'String','y = 0*x;',...
      'Tag','mfmf_fcn',...
    'Visible','off');

elseif (flag) | (~flag & strcmp(get(hmf_multif,'Visible'),'off'))
  set(hmf_multif,'Visible','on');
end

%------ batch direct command in flag param as a string ------------------

waitflag = 1;
if ischar(flag)
  flag = deblank(flag);
  if flag(1) == 'n'
    j = findstr(flag,'=');
    nfunc = flag((j(1)+1):length(flag));
    if str2num(nfunc) > 1
      set(findobj('Tag','mf_mfn'),'String',nfunc);
    end
  end
  if flag(1) == 'f'
    i = sscanf(flag(2:length(flag)),'%i');
    j = findstr(flag,'=');
    fn = flag((j(1)+1):length(flag));
    hedit = get(hmf_multif,'UserData'); % func edit fields
    nhedit = length(hedit)/3;
    if i<= nhedit, set(hedit(3*i-2),'String',fn); end
  end
  if flag(1) == 'y'
    j = findstr(flag,'=');
    fn = flag((j(1)+1):length(flag));
    set(findobj('Tag','mf_mfexpr'),'String',fn);
  end
  if flag(1) == 'c'
    j = findstr(flag,'=');
    fn = flag((j(1)+1):length(flag));
    set(findobj('Tag','mf_mfconst'),'String',fn);
  end
  if flag(1) == 'e'
    set(findobj('Tag','mfmf_extf'),'Value',0);
  end
  if ~isempty(findstr(flag,'getbatch'))
    hedit = get(hmf_multif,'UserData'); % func edit fields
    nhedit = length(hedit)/3;
    if nhedit > 0
      toreturn = sprintf('exec multifunc([],[],''n=%i'');\n',nhedit);
      for i=1:nhedit
        fname = get(hedit(3*i-2),'String');
        if ~isempty(fname)
          toreturn = sprintf('%sexec multifunc([],[],''f%i=%s'');\n',toreturn,i,fname);
        end
      end
      evalconst=get(findobj('Tag','mf_mfconst'),'String');
      if ~isempty(evalconst)
        toreturn = sprintf('%sexec multifunc([],[],''c=%s'');\n',toreturn,evalconst);
      end
      evalfunc = get(findobj('Tag','mf_mfexpr'),'String');
      if ~isempty(evalfunc)
        toreturn = sprintf('%sexec multifunc([],[],''y=%s'');\n',toreturn,evalfunc);
      end
      if  get(findobj('Tag','mfmf_extf'),'Value') == 0
        toreturn = sprintf('%sexec multifunc([],[],''e'');\n',toreturn);
      end
      y = toreturn;
    else
      y = '% no multifunction';
    end
    return;
  end
  if ~isempty(findstr(flag,'hide'))
    set(hmf_multif,'Visible','off')
  end

  waitflag = 0; % no user 'ok/cancel' wait.
end

%------ Make/update funcs name edit fields -------------------------------

hedit = get(hmf_multif,'UserData'); % func edit fields
nhedit = length(hedit)/3;
nfunc=str2num(get(findobj('Tag','mf_mfn'),'String')); % nedit field (user req)

while (nhedit ~= nfunc)
   flag = 'user';
   mf_msg('Update func fields');
   hedit = get(hmf_multif,'UserData'); % func edit fields
   nhedit = length(hedit)/3;
   nfunc=str2num(get(findobj('Tag','mf_mfn'),'String')); % nedit field (user req)
   if ~nfunc, return; end

% resize window
  figpos = get(findobj('Tag','mf_mfunc'),'Position');
  figpos(4) = (4+nfunc*26+28+28+28+28+28);
   set(findobj('Tag','mf_mfunc'),'Position',figpos);
   set(findobj('Tag','mff_text1'),'Position',[50 (4+nfunc*26+28+28+28+28) 130 16]);

% adds func fields if needed
   for i=(nhedit+1):nfunc

     hedit(3*i-2)=uicontrol(hmf_multif,...
    'Style','edit',...
      'BackgroundColor',[1 1 1],...
      'ForegroundColor',[0 0 0],...
      'HorizontalAlignment','left',...
    'Clipping','off',...
    'Position',[40 (4+i*26+28+28+28) 210 20]);



     hedit(3*i-1)  =uicontrol(hmf_multif,...
    'Style','checkbox',...
      'value',1,... %'ForegroundColor',[1 1 1],...
    'String',[ 'f' num2str(i) '='],...
    'ToolTipString',[ 'f' num2str(i) ' expression' ], ...
    'Tag', 'mfmf_click',...
    'Position',[5 (4+i*26+28+28+28) 35 20]);
% 'ToolTipString',[ 'f' num2str(i) ' definition' NL 'if ckecked, will show separate contribution to the final value'],...

hedit(3*i)  = uicontrol(hmf_multif,...
         'Style','pushbutton',...
         'Position',[250 (4+i*26+28+28+28) 20 20],'String','>', 'FontWeight','bold', ...
         'ForegroundColor','blue',...
         'ToolTipString',['Click here to select function ' num2str(i) ' from the Fit Functions menu'],...
         'UserData', {}, ...
         'callback',sprintf('selection(%d);',3*i-2));


    end
% removes func fields if needed
   for i=(nfunc+1):nhedit
  delete(hedit(3*i-2));
  delete(hedit(3*i-1));
    delete(hedit(3*i));
   end
   hedit = hedit(1:(3*nfunc)); % reshape
   nhedit = length(hedit)/3;
   nfunc=str2num(get(findobj('Tag','mf_mfn'),'String')); % nedit field (user req)

   set(hmf_multif,'UserData',hedit);
   delete(gca);
   set(findobj('Tag','mf_mfuncs'),'Visible','on');          % Switch window on


%============ Now wait until ok or cancel pressed =========================

   hok = findobj('Tag','mfmf_ok');
   set(hok,'callback','multifunc([],[],1);');
   if ~flag
  flag = 1;
   end
   nfunc=str2num(get(findobj('Tag','mf_mfn'),'String')); % nedit field (user req)
end % while nfields

if ischar(flag)
  return;
end

if (flag) % make eval string, call fcns...
   [pp, dp, fixed,pnames]=mf_rpars;
   figure(hmf_ctrl);
   mf_msg('OK. Init Multifunc.');
   name = 'MF';
   pin = [];
   pnames = [];
   separate = {};
   hedit = get(hmf_multif,'UserData'); % func edit fields
   nhedit = length(hedit)/3;
   evalstr = ''; evalstr2 = '';
   evalfunc = get(findobj('Tag','mf_mfexpr'),'String');
   evalconst=get(findobj('Tag','mf_mfconst'),'String');
   index = 1;
   if ~isempty(hmf_data) & hmf_data
    figure(hmf_data);
  end

   for i=1:nhedit
  fname = get(hedit(3*i-2),'String');
  fprintf(1,'MF : Expression %i : %s - ',i,fname);
  if isempty(fname)
    fname = '0*x';
  end
  [word1 opt]=strtok(fname);
  if strcmp(word1,'user')
    exprtype = 2;
  else
    exprtype = exist(word1);
  end
  if ~isempty(find(exprtype == [ 0 1 3 4 7 ]))
    isg = evalin('base',[ 'isglobal(' fname ')' ], '-1');
    fprintf(1,'user variable ');
    if (isg  == 1)
      fprintf(1,'(global)\n');
      toadd = [ 'eval(''global ' fname '''); f' num2str(i) '= ' fname ';' ];
      evalstr = [ evalstr toadd ];
      evalstr2 = [ evalstr2 toadd ];
      name = [ name '+' fname ];
    elseif (isg == 0)
      fprintf(1,'(base)\n');
      toadd = [ 'f' num2str(i) '= evalin(''base'',''' fname ''');' ];
      evalstr = [ evalstr toadd ];
      evalstr2 = [ evalstr2 toadd ];
      name = [ name '+' fname ];
    else
      fprintf(1,'expression (might produce error)\n');
      if ~isempty(fname)
        toadd = [ 'f' num2str(i) '= ' fname ');' ];
        evalstr = [ evalstr toadd ];
        evalstr2 = [ evalstr2 toadd ];
        name = [ name '+user'  ];
      end
    end
  elseif ~isempty(find(exprtype == [ 2 5 6 ]))
    if index <= length(p)
      pforfunc = p(index:length(p));
    else
      pforfunc = [];
    end

    fy = strcmp(word1,'user');
    if ~fy
      [fy, ffname, fpnames, fpin] = feval(word1,x,pforfunc,flag);
    else
      fy = fliplr(deblank(fliplr(deblank(opt))));
      if isempty(fy)
        [fy, ffname, fpnames, fpin] = feval(word1,x,pforfunc,'user_ask');
        if ~isempty(fpnames)
          set(hedit(3*i-2),'String',[ word1 ' ' fpnames ]);
        end
      else
        word1 = 'user';
        [fy, ffname, fpnames, fpin] = feval(word1,x,pforfunc,fy);
      end
      opt = '';
    end
    opt = fliplr(deblank(fliplr(deblank(opt))));
    if ~isempty(opt)
      opt = [ ':' opt '.' ];
    else
      opt = ':';
    end

    fpin = fpin(:);
    [lpnames,dummy] = size(fpnames);
    if length(fpin) < lpnames
      fpin = [ fpin ; zeros(lpnames - length(fpin) ,1) ];
    elseif length(fpin) > lpnames
      fpin = fpin(1:lpnames);
    end
    fprintf(1,'fit function (%s) %i parameters\n',ffname, lpnames);
    name = [ name '+' ffname ];
    endindex = index + lpnames -1;
    evalstr = [ evalstr 'f' num2str(i) '=' word1 ];
    evalstr2 = [ evalstr2 '  [f' num2str(i) ',nop,fpn,pf] = ' word1 ];
    if lpnames > 0
      evalstr = [ evalstr '(x,p(' num2str(index) ':' ...
      num2str(endindex) ')); ' ];
      if strcmp(word1,'user') & ~isempty(fpnames)
        evalstr2 = sprintf('%s(x,p(%i:%i),flag);\n  pin = [ pin(:) ; pf(:) ];\n  pnames = strvcat(pnames, strcat(''f%i:%s'',fpn));\n',...
        evalstr2, index, endindex,i,fpnames);
      else
        evalstr2 = sprintf('%s(x,p(%i:%i),flag);\n  pin = [ pin(:) ; pf(:) ];\n  pnames = strvcat(pnames, strcat(''f%i%s'',fpn));\n',...
        evalstr2, index, endindex,i,opt);
      end
    else
      evalstr = [ evalstr '(x);' ];
      evalstr2 = [ evalstr2 '(x,[],flag);' ];
    end

%================= separate structure===========================

    hclick=findobj('tag','mfmf_click');

    j=length(hclick);
    for k=length(hclick):-1:1
      clickk=get(hclick(k),'Value');
      rec=j+1-k;
            click(rec)=clickk;
    end


    if click(i)==1
      evalstr2 =sprintf('%s\n separate{%i}= f%i;\n', evalstr2, i, i);
      %evalstr2 =[evalstr2  'separate{' num2str(i) '}=' ' f' num2str(i) '; ' ];
      %evalstr = [evalstr   'separate{'  num2str(i) '}=' 'f' num2str(i) '; ' ];
      evalstr =sprintf('%s\n separate{%i}= f%i;\n', evalstr, i, i);

    else
      evalstr2 =sprintf('%s\n separate{%i}= {};\n', evalstr2, i);
      evalstr =sprintf('%s\n separate{%i}={};\n', evalstr, i);
    end

%================================================================
    index = endindex+1;
    for j=1:lpnames
      if ((i == 1) & (j==1) | isempty(pnames))
        pnames = [ 'f' num2str(i) opt deblank(mat2str(fpnames(j,:))) ];
      else
        pnames  = strvcat(pnames,[ 'f' num2str(i) opt deblank(mat2str(fpnames(j,:))) ]);
      end
    end
    pin = [ pin(:) ; fpin(:) ];
  else

    fprintf(1,'Unknown expression (might produce error)\n');
    evalstr = [ evalstr 'f' num2str(i) '= ' fname ';' ];
    evalstr2 = [ evalstr2 'f' num2str(i) '= ' fname ';' ];
    name = [ name '+user'  ];
  end

   end

   if isempty(evalfunc)
  evalfunc = 'y=0*x';
   end
   evalstr = [ evalstr  ' y=' evalfunc ';' ];
   evalstr2 = sprintf('%s \n  if flag == 1\n    y=0*x;\n  else\n    y=%s;\n  end', evalstr2, evalfunc);
   if ~isempty(evalconst)
  disp('WARNING: the constrain line will apply, but won''t be visible into parameter window.')
  evalstr = [ evalconst '; ' evalstr ];

   end
   set(findobj('Tag','mfmf_fcn'),'String',evalstr);

   pfix = 0*pin;
   dpmf = get(findobj('Tag','mf_mfn'),'Userdata');
   if flag < 2
  if (length(pp) == length(pin))
      pin = pp;
    pfix = fixed;
  elseif (length(dpmf) == length(pin))
    pin = dpmf(:,1);
    dp = dpmf(:,2);
    if size(dpmf,2) > 2
      pfix = dpmf(:,3);
    else
      pfix = [];
    end
    disp('Restore multifunc old parameters')
  elseif ~isempty(pp)
    disp('Old parameter values are :')
    fprintf(1,'%.3g ',pp(:)); fprintf(1,'\n');
    dp = [];
    end
   end
   p = pin;

   y=[];
   mf_pwin(name,pnames,p);
   if ~isempty(pnames)
  mf_upars([],dp,pfix);
   end
   if (flag == 1)
for i=1:nhedit
  if sum(click)==0
    equifunc = sprintf('function [y, name, pnames, pin] = mf_mftmp(x,p, flag)\n');
  else
    equifunc = sprintf('function [y, name, pnames, pin,separate] = mf_mftmp(x,p, flag)\n');
  end
end
  equifunc = sprintf('%s%% Multifunction Auto-create for : \n', equifunc);
  equifunc = sprintf('%s%% %s\n', equifunc, name);
  t = fix(clock);
  equifunc = sprintf('%s%% date %d.%d.%d   %d:%d:%d\n', equifunc, t(3:-1:1), t(4:6));

%=============================================================
  if click(i)==1
  equifunc=sprintf('%s%% Separate ploting implemented\n', equifunc);
  end

  equifunc = sprintf('%s  name = ''%s'';\n', equifunc, name);
  equifunc = sprintf('%s  pnames = '''';\n', equifunc);
  dps = sprintf('%f ',p(:)');
  equifunc = sprintf('%s  pin = [ %s ];\n', equifunc, dps);
  equifunc = sprintf('%s  separate = {};\n', equifunc);
  equifunc = sprintf('%sif nargin == 2\n',equifunc);
  equifunc = sprintf('%s  %s\n', equifunc, evalstr);
  equifunc = sprintf('%selse\n', equifunc);
  equifunc = sprintf('%s  pin = [];\n', equifunc);
  equifunc = sprintf('%s%s\n', equifunc, evalstr2);
  equifunc = sprintf('%send\n\n%% end of auto multi-function\n', equifunc);

  if get(findobj('Tag','mfmf_extf'),'Value')
    fprintf(1,'Stored mf_mftmp.m temporary function in %s\n',pwd);
    fid = fopen('mf_mftmp.m','w');
    if fid == -1
      disp('Warn : couldn''t write/create file. Using eval mode.')
      set(findobj('Tag','mfmf_extf'),'Value',0);
    else
      fprintf(fid,'%s\n',equifunc);
      fclose(fid);
      clear mf_mftmp
      path(path);
    end
  end
  if get(findobj('Tag','mfmf_extf'),'Value') == 0
    disp('Equivalent function is : ');
    disp(' ');
    if ~isempty(dir('mf_mftmp.m')), delete('mf_mftmp.m'); end
    clear mf_mftmp
    disp(equifunc);
  end
  h = findobj('Tag','mf_FitFuncName');
  if ~isempty(h)
    set(h,'String',name);
  end
  set(findobj('Tag','mf_mfn'),'Userdata', [ p(:) dp(:) pfix(:) ]);
    % store associated starting parameters
   end
   mf_msg('done');
else
   if ~isempty(dir('mf_mftmp.m'))
  if flag
    [y, name, pnames, pin, separate] = feval('mf_mftmp',x,p,flag);
  else
    hclick=findobj('tag','mfmf_click');

    j=length(hclick);
    for k=length(hclick):-1:1
      clickk=get(hclick(k),'Value');
      rec=j+1-k;
            click(rec)=clickk;
    end
    for i=1:nhedit
      if sum(click)==0
        [y, name, pnames, pin] = feval('mf_mftmp',x,p);
      else
        [y, name, pnames, pin,separate] = feval('mf_mftmp',x,p);
      end
    end
  end
   else
  evalstr = get(findobj('Tag','mfmf_fcn'),'String');

    for i=1:size(evalstr,1)
eval(evalstr(i,:));
end

   end
   if length(y) ~= length(x)
  disp('Warn : fix length problem in multifunc. Perhaps wrong')
  if ~isempty(y) & length(x)
    y = y(1:length(x));
  else
    y = 0*x;
  end
   end
end


% function selection(ind)
%
% sel = listdlg('ListString', 'axel', 'ListSize', [300 160], ...
%       'Name', 'Select a directory from the history', 'SelectionMode', 'single');


