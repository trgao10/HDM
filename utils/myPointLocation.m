function [ti,BCs] = myPointLocation(TR,dataPoints)
%MYPOINTLOCATION: Counterpart of pointLocation for (non-Delauney)
%                 triangulations in MATLAB version < 2014b
%   TR:           MATLAB triangulation structure
%   dataPoints:   N-by-D matrx; N is the number of points, D the dimension
%   ti:           index of triangle that contains each dataPoint
%   BC:           barycentric coordinates of each dataPoint in the
%                 enclosing triangle
%
% Tingran Gao, Duke University
% trgao10@math.duke.edu
% last modified: Oct 18, 2014
%

nV = size(dataPoints,1);
nF = size(TR.ConnectivityList,1);
ti = zeros(nV,1);
BCs = zeros(nV,3);
for j = 1:nV
    BC = TR.cartesianToBarycentric((1:nF)',repmat(dataPoints(j,:),nF,1));
    tind = find(all(BC>-1e-10,2));
    if(numel(tind)>0)
        ti(j) = tind(1);
    else
        ti(j) = NaN;
        continue;
    end
    BCs(j,:) = BC(tind(1),:);
end


end

