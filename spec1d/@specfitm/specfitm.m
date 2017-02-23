classdef specfitm < specfit
    %SPECFITM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent)
       pvals
       notfixed
       fit_ind = NaN;
    end
    
    properties
        x_per_spec
        param_keep
        multifit_ind
    end
    
    properties (Hidden)
        notfixed_original
    end
    
    methods
        function obj = specfitm(s,varargin)
            if nargin == 0
                return
            elseif isa(s,'specfit')
                obj.notfixed_original = obj.notfixed;
                obj = obj.copy(s);
            else
                s = specfit(s,varargin{:});
                obj.notfixed_original = obj.notfixed;
                obj = obj.copy(s);
            end
        end
        
        function [obj, s_out] = multifit_ini(obj,s_in)
            %% Multifit initialisation procedure
            % s = arrray of spec1d objects
            % pin = parameters for fitting
            % flag = 0, 1, 2
            % sep = separation of spec1d objects.
            % Very loosly based on the work of MM
            
            
            %----- Create a big spec1d file that contains, in order, all the smaller
            %      ones. Also, memorize the number of points per spec1d file so that
            %      this operation stays reversible!

            obj.x_per_spec = zeros(1,length(s_in));
            
            x = []; y = []; e = [];
            
            for il = 1:length(s_in)
                s_in(il) = clean(s_in(il));
                [xil, yil, eil] = extract(s_in(il));
%                 %----- Remove zeros from e
%                 xil(eil==0) = [];
%                 yil(eil==0) = [];
%                 eil(eil==0) = [];
                obj.x_per_spec(il) = length(xil);
                %----- Fold all data in the same array
                x((length(x)+1):length(x)+length(xil)) = xil(:);
                y((length(y)+1):length(y)+length(yil)) = yil(:);
                e((length(e)+1):length(e)+length(eil)) = eil(:);
            end
            
            % Handle the basic data.
            s_out = s_in(1).copy;
            s_out.x = x;
            s_out.y = y;
            s_out.e = e;
            s_out.yfit = [];
            
            % This deals with additional data in each spectra.
            add_dat = {s_in.userdata};
            for i = 1:length(add_dat)
                add_dat{i} = rmfield(add_dat{i},'rind');
                if isfield(add_dat{i},'combine')
                    if i == 1
                        ns = mergestruct(rmfield(add_dat{i},'combine'),add_dat{i}.combine);
                    else
                        f = fieldnames(add_dat{i}.combine);
                        for j = 1:length(f)
                            if isfield(ns,f{j})
                                ns.(f{j}) = vertcat(ns.(f{j}),add_dat{i}.combine.(f{j}));
                            else
                                ns.(f{j}) = add_dat{i}.combine.(f{j});
                            end
                        end
                    end
                end
            end
            ns.rind = s(1).userdata.rind;
            s_out.userdata = ns;
            
            s_out = feval(class(s(1)),s_out);
            obj.specID = s_out.ident;
        end
        
        function pout = get.pvals(obj)
            pin = @specfit/get.p(obj);
            flag = f_in.notfixed;
            pout = zeros(1,(length(pin)+length(find(flag>1))*(length(obj.x_per_spec)-1)));
%             dpin = pout;
%             obj.param_keep = dpin;
%             
            ll = 1;
            for j = 1:length(obj.notfixed_original)
                if obj.notfixed_original(j) == 2
                    pout(ll:(ll+length(s_in)-1)) = pin_old{j};
%                     dpin(ll:(ll+length(s_in)-1)) = 1;
                    obj.param_keep(ll:(ll+length(s_in)-1)) = j;
                    ll = ll+length(s_in);
                elseif obj.notfixed_original(j) == 3
%                     pout(ll:(ll+length(s_in)-1)) = pin_old{j};
%                     dpin(ll:(ll+length(s_in)-1)) = 0;
%                     obj.(ll:(ll+length(s_in)-1)) = j;
%                     ll = ll+length(s_in);
%                 elseif (dpin_old(j) == 1) || (dpin_old(j) == 0)
%                     pout(ll) = pin_old{j};
%                     dpin(ll) = dpin_old(j);
%                     obj.(ll) = j;
%                     ll = ll+1;
%                 else
%                     disp('Please use following notation: 0 = fixed, 1 = same for all objects, 2 = independently fitted in all objects, 3 = fixed and different for all objects')
%                 end
%             end
    end
