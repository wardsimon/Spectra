function new_file = geturl(filename)
% geturl: Fetch an URL link and create a local copy
%    filename=geturl(URLname)
%    Opens the specified URL, handles it and create a local copy of its content
%    The function returns the name of this copy.
%    Supported URLs include 'http://','ftp://','file://', as well as
%    ZIP and GZIP compressed data extensions (.zip, .gz, .tgz, .gtar)
%
% Part of: iFiles utilities (ILL library)
% Author:  E. Farhi <farhi@ill.fr>. Feb 4th, 2003.

% calls: java.net.URL, java.util.zip, InputStreamByteWrapper


  % open Java URL
  new_file = [];
  if ~exist('usejava')
    disp( 'iFiles/geturl: No Java support available. Can not access network.');
    return;
  end
  disp(['iFiles/geturl: Retrieving URL="' filename '"' ]);
  if ~usejava('jvm')
    disp( 'iFiles/geturl: Error Can not retrieving URL data. Java Virtual Machine is inactive.');
    return
  end

  Jurl = java.net.URL(filename);
  % create buffered java.io.InputStream from java.net.URL.openStream java method
  try
    Jis = openStream(Jurl);
  catch
    disp(['iFiles/geturl: Error: Could not open URL.'])
    disp(lasterr);
    return
  end
  % test if we got a zip or gzip file
  [path,name,ext,versn]=fileparts(filename);
  if ~isempty(strmatch(lower(ext),{'.zip'}))
    Jzip_is = java.util.zip.ZipInputStream(Jis);
    Jzip_is.getNextEntry;
    Jfile = Jzip_is;
    ext = '';
  elseif ~isempty(strmatch(lower(ext),{'.gz','.gzip','.tgz','.gtar'}))
    Jzip_is = java.util.zip.GZIPInputStream(Jis);
    Jfile = Jzip_is;
    ext = '';
  else
    Jfile = Jis;
  end
  try
    % get the length with method java.lang.io.inputStream.available
    Jn_is = available(Jfile);
  catch
    disp(['iFiles/geturl: Error: Could not get size of URL.'])
    disp(lasterr);
    return
  end
  disp(['                ' num2str(Jn_is) ' bytes. Please be patient...' ]);
  if Jn_is > 0
    % open an InputStreamByteWrapper to read file faster
    ReaderSize = 100*1024;  % 100 ko Java buffer
    try
      reader = InputStreamByteWrapper(uint32(ReaderSize));
    catch
      disp(['iFiles/geturl: Error: Can not retrieve data (InputStreamByteWrapper).'])
      return
    end
    % initialize URL content array
    Jbyte_array = uint8([]);
    % display wait bar
    try
      waitbar(0,h);
    catch
      h=[];
    end
    wb_title = ['load: URL=' filename ];
    wb_exists = 0;
    if isempty(h) % create or re-use existing waitbar if any
      h = waitbar(0, wb_title,'CreateCancelBtn','delete(gcf)');
    else
      waitbar(0, h, wb_title);
      wb_exists = 1;
    end
    ht = findall(h,'string', wb_title);
    set(ht,'interpreter','none');
    % create local file
    new_file = [ name ext versn ];
    % remove any existing instances of new_file
    if ~isempty(dir(new_file)), delete(new_file); end
    [fid, message] = fopen(new_file, 'a');
    if fid < 0
      disp(['iFiles/geturl: Error: Could not create local copy ' new_file ' of URL.']);
      disp(message);
      return
    end
    % read full URL content, until read == -1
    WBcounter  = 0;
    len = 1;
    index = 1; barindex=0;
    while (len > 0)
      if index > WBcounter
        WBcounter=WBcounter+Jn_is/10;
        drawnow;
        if ~isempty(ht),
          try
            waitbar(barindex/Jn_is, h);
          catch
            h=[];
          end
        end
        if isempty(h) & ~isempty(ht)
          disp(['iFiles/geturl: aborting URL ' filename ' (' num2str(index) '/' num2str(Jn_is) ').']);
          return
        end
        set(ht,'string', [ 'URL ' filename ' (' num2str(index) '/' num2str(Jn_is) ').']);
        set(h,'name', wb_title);
        if barindex>Jn_is, barindex=0; end
      end
      len = reader.read(Jfile, uint32(ReaderSize)); % Read 4096 bytes at a time (offset=0).
      if (len > 0)
        Jbyte_array = uint8(reader.bfr(1:len)); % Transfert bytes to Matlab buffer.
        index   = index+len;
        barindex=barindex+len;
        % copy data localy
        count = fwrite(fid, Jbyte_array, 'uint8');
        if Jn_is < ReaderSize, Jn_is = ReaderSize; end
        if count ~= len
          disp(['iFiles/geturl: Error: Could not write local copy ' new_file ' of URL.']);
          new_file = [];
          return
        end
      end
    end
    % close inputStream
    if ~isempty(h) & ~wb_exists, close(h); end
    close(Jfile);
    clear Jis Jurl Jfile Jbyte_array reader
    disp(['                wrote local copy "' new_file '" (' num2str(index) ' bytes).']);
    fclose(fid);
  else
    disp('iFiles/geturl: Error: URL contains no data (0 length).');
    return
  end
