function [Transplan12,V1onV2] = MapSoftenBC(Coords1,Coords2,F2,V1,V2)
%MAPSOFTENBC: Produce a transport plan between Coords1 and Coords2
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

flat = @(x) x(:);

nV1 = size(V1,2);
nV2 = size(V2,2);

TR = triangulation(F2',Coords2');
[ti,BC] = pointLocation(TR,Coords1');
NaNInds = find(isnan(ti));
BaseInds = 1:nV1;
BaseInds(NaNInds) = [];
BC(NaNInds,:) = [];
ti(NaNInds) = [];
BCTargets = TR.ConnectivityList(ti,:);
BCPlan = sparse(repmat(BaseInds',3,1),flat(BCTargets),BC(:),nV1,nV2);
V1onV2 = V2*BCPlan';
Transplan12 = BCPlan;

end

