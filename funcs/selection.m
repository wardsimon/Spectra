function selection(ind)
% this function fills multifunc f_n lines with function names

hmf_multif=findobj('Tag','mf_mfunc');
hedit     = get(hmf_multif,'UserData');

if isempty(hedit), return; end
if hedit < ind, return; end

[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

func_names = [];
mat_names  = [];

if isempty(hmf_ctrl) | hmf_ctrl == 0

  func_names = str2mat('Gaussian', 'Gaussian squared','Gaussian (sloping background)','Lorentzian','Lorentzian squared','Airy (for Fabry Perot)','Cusp (double power law)','Damped Harmonic Oscillator','Triangular peak (width fwhm)','Voigt function');
  func_names = str2mat(func_names,'Dirac peak','Heaviside peak','Two Gaussians','Two Lorentzians','Lorentzian + Gaussian','Green function','N Gaussians','N Lorentzians',' Power law (below x0)','Power law (above x0)','Ellipse','Constant (Background)','Straight line (ax+b)');
  func_names = str2mat(func_names,'Quadratic (ax^2+bx+c)','Polynomial','Exponential decay','Reflectivity (N layers)');

  mat_names = str2mat('gauss', 'Gauss2','sgauss','Lorz','lorz2','airyfp','cusp','dho','triangl','voigt','dirac','heaviside','gaussx2','lorzx2','lorgss','green','ngauss','nlorz','pow','pow1','ellipse','background','strline','quadrat','polynomial','expon','refl');

else

  SecTag  = 'Fit Functions';
  NameTag = 'mf_FitFuncName';
  FileTag = 'mf_FitFuncFile';
  DirTag  = 'mf_FitFuncDir';
  MenuTag = 'mf_FitFuncMenu';

  Name=get(findobj('tag',deblank(NameTag)),'userdata');
  Dir =get(findobj('tag',deblank(DirTag)), 'userdata');
  File=get(findobj('tag',deblank(FileTag)),'userdata');

  hmenu   = findobj(hmf_ctrl,'tag',deblank(MenuTag));
  hmchild = get(hmenu,'children');

  for j=1:size(Name,1)
    if ~isempty(deblank(File(j,:)))
      mat_names = strvcat(mat_names,  deblank(File(j,:)));
      func_names= strvcat(func_names, [ Name(j,:) '[' deblank(File(j,:)) ']' ]);
    end
  end

end

% add the 'User' item
mat_names = strvcat(mat_names, 'user parameter1');
func_names= strvcat(func_names, 'User Parameter (change parameter name) [user parameter1]');

% display the dialog list

sel = listdlg('ListString', func_names, 'ListSize', [300 160], ...
      'Name', 'Select a function', 'SelectionMode', 'single');

if sel == 0, return; end

hedit     = get(hmf_multif,'UserData');
set(hedit(ind),'string',deblank(mat_names(sel,:)))
