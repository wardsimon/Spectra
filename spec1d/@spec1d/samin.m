% Calculate the optimum point using SA method similar to samin routine from
% Octave
%   func: function to minimize
%   args: a row vector contain all the function arguments
%   controls: a struct contain all the SA parameters

function [varargout] = samin(varargin)

global pt

j=1;
controls=[];
for i=1:nargin
    temp=varargin{i};
    if isa(temp,'spec1d')
        for k=1:length(temp)
            s(j)=temp(k);
            j=j+1;
        end
    end
    if isa(temp,'struct')
        controls=temp;
    end
end

if isempty(pt)
    pt=0;
else
    if pt>1
        warning('Global pt is not 0 or 1, ignoring command')
        pt=0;
    end
end

func=varargin{2};

pn=1;
if strcmp(func(1),'@')
    func=eval(func);
    pn=0;
end


args=varargin{3};
if nargin>3
    fix=cast(varargin{4},'logical');
else
    fix=ones(size(args));
end

if isempty(controls)
    warning('Upper and lower bounds are set to pin+-1.5*pin')
    controls = struct(...
        'lb',args - (1.5*args),...
        'ub',args + (1.5*args),...
        'nt',320,...
        'ns',3,...
        'rt',0.9,...
        'maxEval',750000*sum(fix),...
        'neps',20,...
        'functol',1.0e-5,...
        'paramtol',0.001,...
        'mfunc','Emc',...
        'gamma',0.1);
end

if pt
    figh=figure;
end

% main settings
lb = controls.lb;
ub = controls.ub;
nt = controls.nt;
ns = controls.ns;
rt = controls.rt;
maxEval = controls.maxEval;
neps = controls.neps;
functol = controls.functol;
paramtol = controls.paramtol;
mfunc = controls.mfunc;
mct=find(strcmp(mfunc,{'Emc','Estun','Estun2','Esinh','Etanh'})==1);
gam = controls.gamma;
for b=1:length(s)
    
    x=s(b).x;
    y=s(b).y;
    e=s(b).e;
    
    p = args;
    bounds = ub - lb;
    n = size(p, 2);
    popt = p;
    
    % initial values
    nacc = 0;           % total accepted trial
    t = 1000;           % initial temperature
    converge = 0;       % convergence indicator
    coverage_ok = 0;    % search space coverage indicator. When this turns 1 the temperature will start to fall
    fstar = 1e9*ones(neps, 1);
    
    % initial obj_values
    f = feval(func,x,args);
    resid=y-f;
    r2=sum(resid.^2)/sum((y-mean(y)).^2);
    ropt = r2;
    feval_f = 1;
    yopt=f;
    fopt_r=NaN(1,maxEval);  fopt_r(1)=1;
    ropt_r=NaN(1,maxEval);  ropt_r(1)=r2;
    iter=1;
    if pt
        hba = waitbar(0,'Initializing waitbar...');
        c=onCleanup(@()delete(hba));
        figure(figh)
        subplot(2,2,1)
        p1_2=bar((1:n),popt(:),'FaceColor','r');
        set(p1_2,'YdataSource','popt')
        xlabel('Pin','Interpreter','tex','Fontname','Helvetica','Fontsize',14)
        ylabel('Pin value','Interpreter','tex','Fontname','Helvetica','Fontsize',14)
        set(gca,'LineWidth',1.5,'Fontname','Helvetica','Fontsize',14,'XaxisLocation','Top')
        subplot(2,2,2)
        p2_1=plot(x,y,'Color','b');
        hold on
        p2_2=plot(x,yopt,'Color','r');
        set(p2_2,'YdataSource','yopt')
        xlabel('X','Interpreter','tex','Fontname','Helvetica','Fontsize',14)
        ylabel('Y value','Interpreter','tex','Fontname','Helvetica','Fontsize',14)
        set(gca,'LineWidth',1.5,'Fontname','Helvetica','Fontsize',14,'XaxisLocation','Top')
        subplot(2,2,3)
        p4_1=plot(fopt_r,ropt_r,'LineStyle','None','Marker','.','Color','b');
        set(p4_1,'XdataSource','fopt_r','YdataSource','ropt_r');
        xlabel('Evaluations','Interpreter','tex','Fontname','Helvetica','Fontsize',14)
        ylabel('Convergence','Interpreter','tex','Fontname','Helvetica','Fontsize',14)
        subplot(2,2,4)       
        p3_1=plot(x,resid);
        set(p3_1,'YdataSource','resid')
        xlabel('X','Interpreter','tex','Fontname','Helvetica','Fontsize',14)
        ylabel('Y-Y_{best value}','Interpreter','tex','Fontname','Helvetica','Fontsize',14)
        set(gca,'LineWidth',1.5,'Fontname','Helvetica','Fontsize',14)
    end
    % main loop, first increase temperature until parameter space covered,
    % then reduce until convergence
    while (converge == 0)
        if sum(fix)==0
            yopt = feval(func,x,p);
            resid=y-yopt;
            r2opt=sum(resid.^2)/sum((y-mean(y)).^2);
            break
        end
        % statistics to report at each temp change, set back to zeros
        nup = 0;
        ndown = 0;
        nrej = 0;
        nnew = 0;
        lnobds = 0;
        % repeat nt times then adjust temperature
        for m = 1:nt
            nacp = zeros(n, 1);
            % repeat ns times then adjust bounds
            for j = 1:ns
                % generate a new points by taking last and adding a random
                % value to each of elements, in turn
                for h = 1:n
                    if fix(h)
                        pp = popt;
                        % This is the parameter sampling area.
                        % ------ NEEDS to be a more adaptive area. --------
                        pp(h) = p(h) + (2*rand-1).*bounds(h);
                        if (pp(h)<lb(h)) || (pp(h)>ub(h))
                            pp(h) = lb(h) + (ub(h) - lb(h))*rand;
                            lnobds = lnobds + 1;
                        end
                        
                        % evaluate function at new point
                        fp = feval(func,x,pp);
                        resid=y-fp;
                        r2p=sum(resid.^2)/sum((y-mean(y)).^2);
                        feval_f = feval_f + 1;
                        % accept the point if the function value decreases
                        if r2p<=r2
                            p = pp;
                            r2 = r2p;
                            nacc = nacc + 1;
                            nacp(h) = nacp(h) + 1;
                            nup = nup + 1;
                            % if lower than the current fopt then record
                            if r2p < ropt
                                ropt = r2p;
                                popt = pp;
                                yopt=fp;                              
                                iter=iter+1;
                                fopt_r(iter)=iter;
                                ropt_r(iter)=ropt;
                                if pt
                                    waitbar(feval_f/maxEval,hba,sprintf('%0.2f%% of maximum time',100*feval_f/maxEval))
                                    refreshdata(figh,'caller')
                                    drawnow
                                end
                                nnew = nnew + 1;
                            end
                        else
                            % the point is higher, use Metropolis criteria to
                            % decide (Emc), or proposals for tunneling in http://dx.doi.org/10.1209%2Fepl%2Fi2006-10058-0
                            switch mct
                                case 1
                                    mc  = @(r2p,r2,gam,t) exp(-(r2p-r2)/t); % Metropolis
                                case 2
                                    mc  = @(r2p,r2,gam,t) (1-exp(-gam*(r2p-r2))); % Stun
                                case 3
                                    mc  = @(r2p,r2,gam,t) log(gam*(r2p-r2) + sqrt(1 + (gam^2 * (r2p-r2)^2))); %Stun2
                                case 4
                                    mc  = @(r2p,r2,gam,t) sinh(gam*(r2p-r2)); % sinh
                                case 5
                                    mc  = @(r2p,r2,gam,t) tanh(gam*(r2p-r2)); % tanh
                            end
                            if rand < mc(r2p,r2,gam,t)
                                p = pp;
                                r2 = r2p;
                                nacc = nacc + 1;
                                nacp(h) = nacp(h) + 1;
                                ndown = ndown + 1;
                            else
                                nrej = nrej + 1;
                            end
                        end
                        
                        % if maxEval exceeded then terminate
                        if feval_f >= maxEval
                            sloop.x=x;
                            sloop.y=y;
                            sloop.e=e;
                            sloop.x_label=s(b).x_label;
                            sloop.y_label=s(b).y_label;
                            sloop.datafile=s(b).datafile;
                            sloop.yfit=yopt;
                            
                            sout(b)=spec1d(sloop);
                            
                            if pn
                                [dummy,dummy,pnames]=feval(func,x,p,1);
                            else
                                for i=1:length(p)
                                    pnames{i}=['p' num2str(i)];
                                end
                            end
                            v=length(y)-length(find(fix));
                            ChiSq = sum(((y-yopt)./e).^2 )/v;
                            fitdata(b).pvals=popt;
                            fitdata(b).evals=(bounds.*fix)/2;
                            fitdata(b).pnames=pnames;
                            fitdata(b).chisq=ChiSq;
                            fitdata(b).rsq=1-ropt;
                            varargout{1}=sout;
                            varargout{2}=fitdata;
                            warning('Convergence criteria are not met. maxEval exceeded.');
                            return
                        end
                    end
                end
            end
            
            % adjust bounds so that approximately half of all evaluations
            % are accepted
            test = 0;
            for i = 1:n
                if fix(i)
                    ratio = nacp(i)/ns;
                    if ratio>0.6
                        bounds(i) = bounds(i) * (1.0 + 2.0 * (ratio - 0.6) / 0.4);
                    elseif ratio<0.4
                        bounds(i) = bounds(i) / (1.0 + 2.0 * ((0.4 - ratio) / 0.4));
                    end
                    % keep within initial bounds
                    if bounds(i) > (ub(i) - lb(i))
                        bounds(i) = ub(i) - lb(i);
                        test = test + 1;
                    end
                end
            end
            if ~coverage_ok
                coverage_ok = (test == sum(fix));
            end
        end
        % check for convergence, if we have cover the params space
        if (coverage_ok)
            % last value close enough to neps value?
            fstar(1) = r2;
            test = (sum(abs(r2-fstar)>functol) > 0);  % conv. failed if test > 0
            % last value close enough to overall best
            if (abs(ropt - r2) <= functol) && (~test)
                for i = 1:n
                    if fix(i)
                        if bounds(i) > paramtol
                            converge = 0;   % no conv. if bounds too wide
                            break;
                        else
                            converge = 1;
                        end
                    end
                end
            end
            
            if (converge > 0)
                break
            end
            
            % reduce the temp, record function value in the neps list and
            % loop again t=rt&8 is the cauc
            t = rt*t;
            for i = neps:-1:2
                fstar(i) = fstar(i-1);
            end
            r2 = ropt;
            p = popt;
        else
            % increase the temp to quickly cover the params space
            t = t^2;
            for i = neps:-1:2
                fstar(i) = fstar(i-1);
            end
            r2 = ropt;
            p = popt;
        end
        
    end
    sloop.x=x;
    sloop.y=y;
    sloop.e=e;
    sloop.x_label=s(b).x_label;
    sloop.y_label=s(b).y_label;
    sloop.datafile=s(b).datafile;
    sloop.yfit=yopt;
    
    sout(b)=spec1d(sloop);
    
    if pn
        [dummy,dummy,pnames]=feval(func,x,p,1);
    else
        for i=1:length(p)
            pnames{i}=['p' num2str(i)];
        end
    end
    v=length(y)-length(find(fix));
    ChiSq = sum(((y-yopt)./e).^2 )/v;
    fitdata(b).pvals=p;
    fitdata(b).evals=(bounds.*fix)/2;
    fitdata(b).pnames=pnames;
    fitdata(b).chisq=ChiSq;
    fitdata(b).rsq=1-ropt;
end

varargout{1}=sout;
varargout{2}=fitdata;
