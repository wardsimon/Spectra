function mf_batch(file)
%
% MFIT function mf_batch(file or command string)
%    Do batch fitting.
% commands :
% par <parindex> <parvalue|expr> set parameter
% fix <parindex> | all           fix parameter
% free <parindex> | all          free parameter
% fit                            fit current data with current function
% load {filename with funcname}  load data using MFit config
% save                           save params, data and fit
% plot                           plot data and fit
% axes                           changes axes limits
% guess (for autoguess)          autoguess before fit
% <mf_var> value|expr            set mf_var value
% exec <command>                 execute Matlab command
% varlist | vars                 show Mfit variables and values
% help                           list of commands
% rem                            comment (same as %)
% stop                           stop batch
%
% For the exec, mf_variables and par command :
% 1- you can use or modify variables x,y,err,selected,p,dp,fixed
% 2- in 'exec' any matlab command is possible such as : set(findobj('tag','mf_dp'),'string',0.1)
%
% M.Zinkin 29.11.94 EF 14.11.97

[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;
if isempty(hmf_ctrl) | ~hmf_ctrl
  mfit;
end

%========= Dialog to choose file if needed ============
if nargin==0

  batchdir=get(findobj('tag','mf_BatchDir'),'string');
  [batchfile batchdir]=uigetfiles(batchdir, 'Select batch file','MultiSelect','off');
  if batchfile==0
    return;
  end
  file=[batchdir batchfile];
  set(findobj('tag','mf_BatchDir'),'string',batchdir);
end

%--------- Initialize some variables ------------------
endoffile=-1;
endofpars=1;
cont=0;
eol = [ 10 13]; % end of line chars

name=get(findobj('tag','mf_VarList'),'userdata');

if iscellstr(file)
  file = char(file);
end

%------------- Open file -------------------------------
fid = -1;
if length(file < 32)
    fid=fopen(file,'r');
end
flag=cont;
ln=0;
if fid ~= -1
  mf_msg('Running batch file...');
  disp(['*Opening batch file ' file]);
  line=fgetl(fid);
elseif isstr(file)
  fprintf(1,'Executing command string (%i chars)\n',length(file))
  if size(file,1) > 1
    fid = '';
    for i=1:size(file,1)
      fid = sprintf('%s\n%s',fid,file(i,:));
    end
  else
    fid = file;
  end
  [line, fid] = strtok(fid,eol);
else
  disp('unknown batch command format')
  return
end

% get MFit data -------------------------------

if ~isempty(hmf_data) & hmf_data
  userdata = get(hmf_data,'userdata');
  x=userdata(:,1);
  y=userdata(:,2);
  err=userdata(:,3);
  selected=userdata(:,4);
end
[p,dp,fixed] = mf_rpars;


%======= Execute file, line by line =====================

ln=1;

while (~any(line==-1))     % While not end-of-file
  if (line(1) == '%')
    line = [ 'rem ' line ];
  end

%------- Break line into 'words'-space delimited strings -------------
  linesav = line;
% f=(line=='=');                               % Change '=' into ' '
% if any(f~=0)
%     line(f)=setstr(' '*ones(size(line(f))));
% end
  hpars=get(findobj('tag','mf_ParWindow'),'userdata');
  [word1 val]=strtok(line,' '); wordo = word1;
  [word2 line]=strtok(val,' ');
  [word3 line]=strtok(line,' ');
  [word4 rest]=strtok(line,' ');
  word1=lower(word1);                          % Make commands case insensitive

%======= Set a parameter =============================================
  if strcmp(word1,'par')
    pindex=str2num(['0' word2]);           % get parameter number and value...
    pval=eval( word3 ,'[]');
    if isempty(pval)                       % (parameter name specified, so value is word 4)
    pval=eval( word4 ,'0');
  end
  p(pindex)=pval;                        % set the parameter
  mf_upars(p,[]);                        % update parameter list
  disp(sprintf('%4d: set par %d to %.4e',ln,pindex,pval));


%======== Process commands ============================================
  else

%----- Load data file and fit functions --------------------------------
    if strcmp(word1,'load') & hmf_ctrl
      if ~isempty(word2)
        set(findobj('Tag','mf_DataFile'),'String',word2)
        if strcmp(word3,'with')
          if ~isempty(word4)
            set(findobj('Tag','mf_LoadRoutineFile'),'String',word4)
          end
        end
      end
      fprintf(1,'%4d: loading data file %s with %s\n',ln,get(findobj('Tag','mf_DataFile'),'String'),get(findobj('Tag','mf_LoadRoutineFile'),'String'));

      mf_gdata('noprompt');
      mf_newfn('noprompt');
                  userdata = get(hmf_data,'userdata');
%---I inserted this next bit, ARW 1-12-98
                  x=userdata(:,1);
                  y=userdata(:,2);
                  err=userdata(:,3);
                  selected=userdata(:,4);

%----- Do a fit -------------------------------------------------------
    elseif strcmp(word1,'fit') & hmf_data
      disp(sprintf('%4d: do fit',ln));
      mf_dofit;
      [p,dp,fixed]= mf_rpars;

%----- Fix a parameter ------------------------------------------------
    elseif strcmp(word1,'fix') & ~isempty(hpars)
      if strcmp(word2,'all')
        pindex = 1:size(hpars,1);
      else
        pindex=str2num(['0' word2]);
      end
      fprintf(1,'%4d: fix parameter ',ln);
      fprintf(1,'%d ', pindex);
      fprintf(1,'\n');
      fixed(pindex) = 1;
      set(hpars(pindex,3),'Value',1);

%----- Free a parameter -----------------------
    elseif strcmp(word1,'free') & ~isempty(hpars)
      if strcmp(word2,'all')
        pindex = 1:size(hpars,1);
      else
        pindex=str2num(['0' word2]);
      end
      fprintf(1,'%4d: free parameter ',ln);
      fprintf(1,'%d ', pindex);
      fprintf(1,'\n');
      fixed(pindex) = 0;
      set(hpars(pindex,3),'Value',0);

%----- Save parameters -----------------------
    elseif strcmp(word1,'save')==1
      if ~isempty(findstr('data',word2))
        disp(sprintf('%4d: saving data',ln));
        mf_save('data');
      elseif ~isempty(findstr('curve',word2))
      disp(sprintf('%4d: saving curve',ln));
      mf_save('curve');
    else
      disp(sprintf('%4d: saving parameters',ln));
      mf_save('parameters');
    end

%----- Pause -----------------------
    elseif strcmp(word1,'pause')
      secs=str2num(['0' word2]);
      disp(sprintf('%4d: pause - hit a key to continue',ln));
      if secs>0
        pause(secs);
      else
        pause;
      end

%----- Update plot -----------------------
    elseif strcmp(word1,'plot')==1 & hmf_data
      disp(sprintf('%4d: Update plot',ln));
      mf_uplot('all');

%----- Change axes limits -----------------------
    elseif strcmp(word1,'axes')==1 & hmf_data
      disp(sprintf('%4d: Setting axes limits',ln));
      mf_xylim(word2);

%----- AutoGuess command -----------------------
    elseif strcmp(word1,'guess')==1 & hmf_data
      disp(sprintf('%4d: Auto guess command',ln));
      mf_guess(1);

%----- Execute command -----------------------
    elseif strcmp(word1,'exec')==1
      [word1 val]=strtok(linesav,' ');
      disp(sprintf('%4d: execute command %s',ln,val));

      eval(val);

    elseif strcmp(word1,'help')==1
      help mf_batch

    elseif strcmp(word1,'news')==1
      edit('Install.txt')

    elseif strcmp(word1,'rem')==1
      disp(sprintf('%4d: comment %s',ln,linesav));

    elseif strcmp(word1,'stop')==1
      disp(sprintf('%4d: stop',ln));
      line = -1;
      fid = -1;

%----- See if it's a command to change a file or dir or mf_value...
    else
      for i=1:size(name,1);
        var=deblank(name(i,:));

        if strcmp(wordo,var)

          f=(val=='=');                               % Change '=' into ' '
          if any(f~=0)
            val(f)=setstr(' '*ones(size(val(f))));
          end
          val=fliplr(deblank(fliplr(deblank(val))));
          h=findobj('tag',['mf_' var]);
          par=get(h,'userdata');
          if ~isempty(val)
            disp(sprintf('%4d: setting %s to %s',ln, wordo,val));
            if isempty(par) | size(par,1)>1
              set(h,'string',val);
            else
              set(h,'userdata',val);
            end
          else
            if isempty(par) | size(par,1)>1
              par = get(h,'string');
            else
              par = get(h,'userdata');
            end
            disp(sprintf('%4d: getting %s is : %s',ln, wordo,par));
          end
        end
      end
      if ~isempty(findstr(word1,'varlist')) | ~isempty(findstr(word1,'vars'))
        disp(sprintf('%4d: getting %s: \n',ln, word1));
        for i=1:size(name,1)
          h=findobj('tag',['mf_' deblank(name(i,:))]);
          par=get(h,'userdata');
          if isempty(par) | size(par,1)>1
            par=get(h,'string');
          end
          fprintf(1,'%s \t = %s\n',name(i,:),par);
        end
      end
    end

% update MFit data -------------------------------

    if ~isempty(hmf_data) & hmf_data
      userdata(:,1) = x;
      userdata(:,2) = y;
      userdata(:,3) = err;
      userdata(:,4) = selected;
      set(hmf_data,'userdata',userdata);
      mf_upars(p,dp,fixed);
      mf_gdata('noload+noplot');
    end
  end

  if isstr(fid)
    [line, fid] = strtok(fid,eol);
    if ~isempty(line)
      ln=ln+1;
    else
      fid = -1;
      line = -1;
    end
  elseif fid ~= -1
    line=fgetl(fid);
    ln=ln+1;
  end
end

%------- Close the file --------------------------------
if fid ~= -1
  disp('*Batch file complete.')
  mf_msg('Batch file complete');
  fclose(fid);
end
%mf_upars(p,dp,fixed);
mf_gdata('noload');
