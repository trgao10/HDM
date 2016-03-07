function [Transplan12,Kernel12,AugKernel12,V1onV2] = MapSoftenKernel(Coords1,Coords2,F2,V1,V2,epsilon,augParam)
%MAPSOFTENKERNEL: Produce a transport plan between Coords1 and Coords2
%   Coords1, Coords2: matrices of dimensions 2-x-nV1 and 2-x-nV2
%   F2:               connectivity list of faces on Coords2
%   V1, V2:           vertices parametrized by Coords1, Coords2
%   epsilon:          diffusion bandwidth; 'auto' uses the mean edge length
%                     of the triangulation (V2,F2)
%   Transplan12:      matrix of dimension nV1-x-nV2, which means that each
%                     row specifies how to transport a unit mass on
%                     Coords1 to Coords2
%
%   Tingran Gao, Duke University
%   trgao10@math.duke.edu
%   last modified: Oct 15, 2014
%

if nargin<7
    augParam = 1.5;
end

flat = @(x) x(:);

nV1 = size(V1,2);
nV2 = size(V2,2);

TR = triangulation(F2',Coords2');
[ti,BC] = pointLocation(TR,Coords1');
% [ti,BC] = myPointLocation(TR,Coords1');
NaNInds = find(isnan(ti));
BaseInds = 1:nV1;
BaseInds(NaNInds) = [];
BC(NaNInds,:) = [];
ti(NaNInds) = [];
BCTargets = TR.ConnectivityList(ti,:);
BCPlan = sparse(repmat(BaseInds',3,1),flat(BCTargets),BC(:),nV1,nV2);
V1onV2 = V2*BCPlan';

V1onV2(:,NaNInds) = [];

BCTargets_x = reshape(V2(1,BCTargets),size(BCTargets));
BCTargets_y = reshape(V2(2,BCTargets),size(BCTargets));
BCTargets_z = reshape(V2(3,BCTargets),size(BCTargets));

Dists2BCTargets = sqrt((repmat(V1onV2(1,:),3,1)-BCTargets_x').^2 +...
                       (repmat(V1onV2(2,:),3,1)-BCTargets_y').^2 +...
                       (repmat(V1onV2(3,:),3,1)-BCTargets_z').^2);

Threshold = zeros(1,nV1);
Threshold(BaseInds) = augParam*mean(Dists2BCTargets);
V1onV2 = V2*BCPlan';
clear BCPlan

if strcmpi(epsilon, 'auto')
    Edges = TR.edges;
    epsilon = median(sum((V2(:,Edges(:,1))-V2(:,Edges(:,2))).^2))/2;
end

Kernel_Dists2BCTargets = exp(-Dists2BCTargets.^2/epsilon);
Kernel12 = sparse(repmat(BaseInds',3,1),flat(BCTargets),flat(Kernel_Dists2BCTargets'),nV1,nV2);
Kernel_Dists2BCTargets = repmat(1./sum(Kernel_Dists2BCTargets),3,1).*Kernel_Dists2BCTargets;
Transplan12 = sparse(repmat(BaseInds',3,1),flat(BCTargets),flat(Kernel_Dists2BCTargets'),nV1,nV2);

DistMatrix = pdist2(V1onV2',V2');
% keyboard
DistMatrix(DistMatrix>repmat(Threshold',1,size(DistMatrix,2))) = Inf;
AugKernel12 = sparse(exp(-DistMatrix.^2/epsilon));

end

