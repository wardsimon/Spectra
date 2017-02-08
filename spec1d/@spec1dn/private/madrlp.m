function obj = madrlp(obj)
    
sample = obj.resolution.sample;
sample = makeG(sample);
sample = makeVolume(sample);
sample = makeReciprical(sample);
sample = makeB(sample);
obj.crystal = sample;
EXP = obj.resolution.EXP;

     
    % Calculates the matrix S to convert from HKL to Scatting Plane    
    aspv=[EXP.orient1(:)';EXP.orient2(:)'];
    aspv(3,:)=0;
    vv = (sample.bmat*aspv')';
    vv(3,:)=cross(vv(1,:)',vv(2,:)')';
    vv(2,:)=-cross(vv(1,:)',vv(3,:)')';
    c = sum(vv'.^2);
    if sum(abs(c)<eps)>0;
        error('Failed to calculate C');
    end
    c = sqrt(c);
    vv = vv./repmat(c',1,3);
    obj.resolution.rl2sp = vv*sample.bmat;
end

        
function obj = makeReciprical(obj)
obj.astar = obj.b * obj.c * sind(obj.alpha) / obj.V;
obj.bstar = obj.a * obj.c * sind(obj.beta) / obj.V;
obj.cstar = obj.a * obj.b * sind(obj.gamma) / obj.V;
obj.alphastar = acosd((cosd(obj.beta) * cosd(obj.gamma) - cosd(obj.alpha)) / (sind(obj.beta) * sind(obj.gamma)));
obj.betastar  = acosd((cosd(obj.alpha) * cosd(obj.gamma) - cosd(obj.beta)) / (sind(obj.alpha) * sind(obj.gamma)));
obj.gammastar = acosd((cosd(obj.alpha) * cosd(obj.beta) - cosd(obj.gamma)) / (sind(obj.alpha) * sind(obj.beta)));
end

function obj = makeB(obj)
%             a=[obj.as obj.bs obj.cs];
%             alpha=[obj.aa obj.bb obj.cc];
%
%             a = a/(2*pi);
%
%             cosa=cosd(alpha);
%             sina=sind(alpha);
%
%             c=sqrt(1 - sum(cosa.^2) + 2*prod(cosa));
%             b=sina./(a*c);
%
%             cosb=zeros(1,3);
%             for i=1:3
%                 j=mod(i,3)+1;
%                 k=mod(j,3)+1;
%                 cosb(i)=(cosa(j)*cosa(k) - cosa(i))/(sina(j)*sina(k));
%             end
%             sinb=sqrt(1-cosb.^2);
%
%             obj.bmat=[ b(1) b(2)*cosb(3) b(3)*cosb(2);...
%                 0    b(2)*sinb(3) -b(3)*sinb(2)*cosa(1);...
%                 0    0            1/a(3)];
obj.bmat = 2*pi*[obj.astar, obj.bstar*cosd(obj.gammastar), obj.cstar * cosd(obj.betastar);
    0, obj.bstar*sind(obj.gammastar), -obj.cstar*sind(obj.betastar)*cosd(obj.alpha);
    0, 0, 1/obj.c];
end

        function obj = makeVolume(obj)
%             obj.V = obj.as*obj.bs*obj.cs*sqrt(1 - cosd(obj.aa)^2 - cosd(obj.bb)^2 - cosd(obj.cc)^2 +2*cosd(obj.aa)*cosd(obj.bb)*cosd(obj.cc));
              obj.V = sqrt(det(obj.G));
        end
        
        function obj = makeG(obj)
              obj.G = [obj.a^2, obj.a*obj.b*cosd(obj.gamma), obj.a*obj.c*cosd(obj.beta);
                       obj.a*obj.b*cosd(obj.gamma), obj.b^2, obj.b*obj.c*cosd(obj.alpha);
                       obj.a*obj.c*cosd(obj.beta), obj.b*obj.c*cosd(obj.alpha), obj.c^2];
        end