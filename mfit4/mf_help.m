function mf_help(app)

  if nargin == 0, app = []; end
  if isempty(app), app='mfit'; end

  m_location    = fileparts(which(mfilename));
  this_location = '';

  % ./doc/mfit.html
  if isempty(this_location)
    this_location = [ m_location filesep ...
    'doc' filesep app '.html' ];
  end

  % ./doc/mfit4/mfit.html
  if isempty(dir(this_location))
    this_location = [ m_location filesep ...
    'doc' filesep 'mfit4' filesep app '.html' ];
  end

  % ../../doc/mfit.html
  if isempty(dir(this_location))
    this_location = [ m_location filesep ...
    '..' filesep '..' filesep 'doc' filesep app '.html' ];
  end

  % ../../doc/mfit4/mfit.html
  if isempty(dir(this_location))
    this_location = [ m_location filesep ...
    '..' filesep '..' filesep 'doc' filesep 'mfit4' filesep app '.html' ];
  end

  % ../doc/mfit.html
  if isempty(dir(this_location))
    this_location = [ m_location filesep ...
    '..' filesep 'doc' filesep app '.html' ];
  end

  % ../doc/mfit4/mfit.html
  if isempty(dir(this_location))
    this_location = [ m_location filesep ...
    '..' filesep 'doc' filesep 'mfit4' filesep app '.html' ];
  end

  if ~isempty(dir(this_location))
    disp([ app ': URL=file:' this_location ]);
  else
    this_location = 'http://www.ill.fr/tas/matlab/doc';
    disp([ app ': URL=file:' this_location ]);
  end

web(this_location);