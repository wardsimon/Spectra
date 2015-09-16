function [pval, perr, pnames]=mftconv(file, outfile)
%
%MFTCONV convert .mft file to multicolumn format
%
%[pval perr pnames]=mftconv(mftfile, [outfile])
%  reads .mft file mftfile (as produced by MFIT), and produces matrix of
%  parameters: pval is an nfits*npars matrix containing the parameter values;
%  perr contains the corresponding uncertainties. If an output file is
%  specified, the results will be saved - to a .mat file if the filename has
%  the .mat extension, and to a multicolumn text file with a header otherwise.
%  NB. All entries in the .mft file must be from the same fit function
%  or the output will be meaningless.
%
%  MZ 5.6.95

%------------- Open data file-----------------------------------------
fid=fopen(file,'r');
if (fid<0)
  error(['Couldn''t find ' file]);
end

%------------- Now read file line by line ------------------------
line=fgetl(fid);
row=1; col=1;
pnames='junk';

while line~=-1

   while (~strcmp(strtok(line,' '),'par') & line~=-1)       % read to next 'par'
    line=fgetl(fid);
   end

   col=1;
   while (strcmp(strtok(line,' '),'par') & line~=-1)        % while lines start with 'par'
      if row==1
         pnames=str2mat(pnames, sscanf(line,'%*s%*s%s',1));  % extract name
      end
      pval(row,col)=eval(sscanf(line,'%*s%*s%*s%s',1));            % extract value
      perr(row,col)=eval(sscanf(line,'%*s%*s%*s%*s%s',1));         % extract uncertainty
    line=fgetl(fid);
    col=col+1;
   end
   while (isempty(findstr(line,'Chi')) & line~=-1)           % read to line with 'Chi' in
      line=fgetl(fid);
   end
   if line~=-1
      i=findstr(line,'Chi');
      chisq(row)=sscanf(line(i:length(line)),'%*s%f');
   end

   row=row+1;
end

pnames(1,:)=[];                                             % Cut out junk pname
fclose(fid);

%------------- Write output file if requested -----------------------
if nargin==2                                                % Output file specified

   if ~isempty(dir(outfile))
      error('File exists.')
   end
   if findstr(outfile,'.mat')                               % Save .mat file
      eval(['save ' outfile ' pnames pval perr chisq']);
   else                                                     % Save multicolumn text file

      fid=fopen(outfile,'w');                               % Open output file
      if fid==-1
         error('Bad output file')
      end

      fprintf(fid,['Data from: ' file '\n']);
      if size(pnames,2)>11                                  % Check parameter names <=11 chars...
         pnam=pnames(:,1:11);                               % ...if not trim to 11
      else
         pnam=pnames;
      end

      for i=1:size(pnam,1)                                % Save column headings
%         nam=[pnam(i,:) '           '];
         fprintf(fid,'  %11s    e:%11s',pnam(i,:),pnam(i,:));
      end
      fprintf(fid,'%11s','Chi^2');
      for i=1:size(pval,1)                                  % Save parameter/error pairs
         fprintf(fid,'\n');
         p=[pval(i,:); perr(i,:)];
         p=[reshape(p,1,2*size(pval,2)) chisq(i)];
         fprintf(fid,'%15.6e',p);
      end

      fclose(fid);
   end
end


