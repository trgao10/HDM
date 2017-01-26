%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HDM for a pilot morphological dataset of 50 lemur teeth (5 groups, each
%%% containing 10 specimens).
%%% This is the dataset featured in "Diffusion Geometry of Fibre Bundles".
%%% last modified: August 13, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
BaseEps = 0.04;
BNN = 4;
FibrEps = 1e-3;
MapType = 'cPMST';
FeatureFix = '';
GroupLevel = 'Genus';
GroupNames = {'Alouatta','Ateles','Brachyteles','Callicebus','Saimiri'};

%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/HDM/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
sample_path = '../DATA/HDM/samples/';
result_path = '../ArchivedResults/HDM/cPMST/FeatureFixOff/';
soften_path = [result_path 'soften/'];

%% load taxa codes
taxa_file = [data_path 'hdm_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);
ChunkSize = 25;

%% options that control the diffusion eigenvector visualization
options.sample_path = sample_path;
options.DisplayLayout = [5,10];
options.DisplayOrient = 'Horizontal';
options.boundary = 'on';
options.names = 'off';
options.linkCamera = 'on';

%% useful inline functions
ChunkIdx = @(TAXAind1,TAXAind2) ceil(((TAXAind1-1)*GroupSize+TAXAind2)/ChunkSize);

%% parse GroupNames
[~,ClTable,~] = xlsread(spreadsheet_path);
Names = {};
NamesByGroup = cell(1,length(GroupNames));
TaxaByGroup = cell(1,length(GroupNames));
for j=1:length(GroupNames)
    NamesJ = ClTable(strcmpi(ClTable(1:end,strcmpi(ClTable(1,:),GroupLevel)),GroupNames{j}),1);
    Names = [Names,NamesJ{:}];
    NamesByGroup{j} = NamesJ;
    TaxaByGroup{j} = cellfun(@(x) find(strcmpi(taxa_code, x)), NamesJ);
end

GroupSize = length(Names);
DiffMatrixSizeList = zeros(GroupSize,1);
TAXAinds = zeros(GroupSize,1);
NamesDelimit = zeros(GroupSize+1,2);
meshList = cell(1,GroupSize);
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,Names{j}));
    load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']);
    meshList{j} = G;
    DiffMatrixSizeList(j) = G.nV;
    NamesDelimit(j+1,1) = NamesDelimit(j,2)+1;
    NamesDelimit(j+1,2) = NamesDelimit(j+1,1)+G.nV-1;
end
Names = taxa_code(TAXAinds); % match upper/lower cases
NamesDelimit(1,:) = [];
nVList = DiffMatrixSizeList;
nVListCumsum = cumsum(nVList);

PerGroupSize = zeros(1,length(GroupNames));
for j=1:length(NamesByGroup)
    for k=1:length(NamesByGroup{j})
        NamesByGroup{j}{k} = taxa_code{strcmpi(taxa_code,NamesByGroup{j}{k})};
    end
    PerGroupSize(j) = length(NamesByGroup{j});
end
CumsumPerGroupSize = cumsum(PerGroupSize);

PerGroupDelimit = [[1,CumsumPerGroupSize(1:end-1)+1]', CumsumPerGroupSize'];
colorsList = [228,26,28;0,126,204;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191]/255;

%% collection rigid motions
rigid_motions = load([data_path 'cPMSTinitRms.mat']);
R = rigid_motions.R(TAXAinds,TAXAinds);

options.R = cell(1,GroupSize);
for j=1:length(options.R)
    options.R{j} = R{1,j};
end
% options.R = sdpSyncRoundOff(cell2mat(R), 3, GroupSize, struct('isreal',1,'issym',1,'maxit',1000,'disp',0));
% options.R = specSyncRoundOff((cell2mat(R)+cell2mat(R)')/2, 3, GroupSize, struct('isreal',1,'issym',1,'maxit',1000,'disp',0));
% options.R = rigid_motions.R(TAXAinds,TAXAinds);

%% process base diffusion
load([result_path MapType 'DistMatrix.mat']);
if strcmpi(MapType,'cP')
    eval('BaseDistMatrix = cPDistMatrix(TAXAinds,TAXAinds);');
else
    eval('BaseDistMatrix = ImprDistMatrix(TAXAinds,TAXAinds);');
end
BaseDistMatrix = BaseDistMatrix-diag(diag(BaseDistMatrix));

%%% only connect BNN-nearest-neighbors
[sDists,rowNNs] = sort(BaseDistMatrix,2);
sDists = sDists(:,2:(1+BNN));
rowNNs = rowNNs(:,2:(1+BNN));
BaseWeights = sparse(repmat((1:GroupSize)',1,BNN),rowNNs,sDists);
BaseWeights = min(BaseWeights, BaseWeights');
for j=1:GroupSize
    sDists(j,:) = BaseWeights(j,rowNNs(j,:));
end
sDists = exp(-sDists.^2/BaseEps);

%% build diffusion kernel matrix
DiffMatrixSize = sum(DiffMatrixSizeList);
DiffMatrixSizeList = cumsum(DiffMatrixSizeList);
DiffMatrixSizeList = [0; DiffMatrixSizeList];
GroupDelimit = zeros(length(GroupNames)+1,2);
for j=2:(length(GroupNames)+1)
    GroupDelimit(j,1) = GroupDelimit(j-1,2)+1;
    GroupDelimit(j,2) = DiffMatrixSizeList(CumsumPerGroupSize(j-1)+1);
end
GroupDelimit(1,:) = [];
DiffMatrixSizeList(end) = []; % treated as block shifts
DiffMatrixRowIdx = [];
DiffMatrixColIdx = [];
DiffMatrixVal = [];

cback = 0;
for j=1:GroupSize
    G1 = load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']); G1 = G1.G;
    for nns = 1:BNN
        if (sDists(j,nns) == 0)
            continue;
        end
        k = rowNNs(j,nns);
        G2 = load([sample_path taxa_code{strcmpi(taxa_code,Names{k})} '.mat']); G2 = G2.G;
        
        %%% load texture coordinates
        TAXAind1 = TAXAinds(j);
        TAXAind2 = TAXAinds(k);
        load([soften_path 'soften_mat_' num2str(ChunkIdx(TAXAind1, TAXAind2)) '.mat']);
        AugKernel12 = cPSoftMapsMatrix{TAXAind1, TAXAind2};
        
        [rowIdx, colIdx, val] = find(AugKernel12);
        DiffMatrixRowIdx = [DiffMatrixRowIdx; rowIdx+DiffMatrixSizeList(j)];
        DiffMatrixColIdx = [DiffMatrixColIdx; colIdx+DiffMatrixSizeList(k)];
        DiffMatrixVal = [DiffMatrixVal; sDists(j,nns)*val];

        [rowIdx, colIdx, val] = find(AugKernel12');
        DiffMatrixRowIdx = [DiffMatrixRowIdx; rowIdx+DiffMatrixSizeList(k)];
        DiffMatrixColIdx = [DiffMatrixColIdx; colIdx+DiffMatrixSizeList(j)];
        DiffMatrixVal = [DiffMatrixVal; sDists(j,nns)*val];
    end
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf(['%4d/' num2str(GroupSize) ' done.\n'],j);
end

H = sparse(DiffMatrixRowIdx,DiffMatrixColIdx,DiffMatrixVal,DiffMatrixSize,DiffMatrixSize);
clear DiffMatrixColIdx DiffMatrixRowIdx DiffMatrixVal rowIdx colIdx val
clear TextureCoords1Matrix TextureCoords2Matrix

%% eigen-decomposition
sqrtInvD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,1./sqrt(sum(H)));
H = sqrtInvD*H*sqrtInvD;
H = (H+H')/2;

eigopt = struct('isreal',1,'issym',1,'maxit',5000,'disp',0);
tic;
[U, lambda] = eigs(H, 101, 'LM', eigopt);
lambda = diag(lambda);
disp(['Eigs completed in ' num2str(toc) ' seconds']);

%%
%==========================================================================
%%% HBDM (Horizontal Base Diffusion Maps)
%==========================================================================
sqrtInvD(isinf(sqrtInvD)) = 0;
BundleHDM = sqrtInvD*U(:,2:end);
HBDM = zeros(GroupSize, nchoosek(size(BundleHDM,2),2));
for j=1:GroupSize
    BundleHDM_Block = normc(BundleHDM(NamesDelimit(j,1):NamesDelimit(j,2),:));
    BundleHDM_Block = BundleHDM_Block*sparse(1:(size(U,2)-1), 1:(size(U,2)-1), sqrt(lambda(2:end)));
    HBDM(j,:) = pdist(BundleHDM_Block', @(x,y) y*x');
end
% Yhbdm = tsne(HBDM, [], 3, 50, 30);
HBDM_dist = pdist(HBDM);
% save('HBDM', 'HBDM');
[Yhbdm,stress] = mdscale(HBDM_dist,4,'criterion','metricstress');

figure('Name','HBDM');
% axis equal
% hold on
for j=1:length(GroupNames)
    plot3(-Yhbdm(PerGroupDelimit(j,1):PerGroupDelimit(j,2),1),...
        Yhbdm(PerGroupDelimit(j,1):PerGroupDelimit(j,2),3),...
        Yhbdm(PerGroupDelimit(j,1):PerGroupDelimit(j,2),2),...
        'Color', colorsList(j,:), 'Marker', '.', 'MarkerSize', 15, 'LineStyle', 'none');
    if (j == 1)
        grid on
        axis equal
        hold on;
    end
end
% view([-20,10]);
% set(gcf,'color','w');
legend(GroupNames);
%%%% ATTENTION! %%%%
%%%% The rows/columns of "BaseDistMatrix" do not follow "taxa_code"!
%%%% Instead, they follow variable "Names", which is build on group orders.
% nameMode = datacursormode(gcf);
% set(nameMode,'DisplayStyle','window');
% set(nameMode, 'UpdateFcn', {@showMeshName, Yhbdm, Names});

%%
%==========================================================================
%%% MDS on BaseDistMatrix
%==========================================================================
[Yb,stress] = mdscale(BaseDistMatrix,3,'criterion','metricstress');
figure('Name','BaseDistMatrix');
for j=1:length(GroupNames)
    plot3(Yb(PerGroupDelimit(j,1):PerGroupDelimit(j,2),1),...
        Yb(PerGroupDelimit(j,1):PerGroupDelimit(j,2),2),...
        Yb(PerGroupDelimit(j,1):PerGroupDelimit(j,2),3),...
        'Color', colorsList(j,:), 'Marker', '.', 'MarkerSize', 15, 'LineStyle', 'none');
    if (j == 1)
        grid on
        axis equal
        axis tight
        hold on;
    end
end
legend(GroupNames);

%%
%==========================================================================
%%% Diffusion Maps on BaseDistMatrix
%==========================================================================
tD = BaseDistMatrix;
tD = tD+diag(Inf(GroupSize,1));
epsilon = mean(min(tD, [] ,2));
W = exp(-tD.^2/epsilon^2);
D = sum(W,2);
L = diag(1./sqrt(D))*W*diag(1./sqrt(D));
L = (L+L')/2;

[Udm, Ldm] = eigs(L, GroupSize-1, 'LM', eigopt);
dims = sum(diag(Ldm)>0);
Udm = Udm(:,1:dims);
Ldm = Ldm(1:dims, 1:dims);
Ydm = diag(1./sqrt(D))*Udm*sqrt(Ldm);

% Ydm = tsne(Ydm, [], 3, 10, 30);
DM_dist = pdist(Ydm);
[Ydm,stress] = mdscale(DM_dist,3,'criterion','metricstress');

figure('Name','Diffusion Maps');
for j=1:length(GroupNames)
    plot3(Ydm(PerGroupDelimit(j,1):PerGroupDelimit(j,2),1),...
        Ydm(PerGroupDelimit(j,1):PerGroupDelimit(j,2),2),...
        Ydm(PerGroupDelimit(j,1):PerGroupDelimit(j,2),3),...
        'Color', colorsList(j,:), 'Marker', '.', 'MarkerSize', 15, 'LineStyle', 'none');
    if (j == 1)
        grid on
        axis equal
        hold on;
    end
end
legend(GroupNames);

%%
%==========================================================================
%%% consistent spectral clustering on each surface
%==========================================================================
% numClusters = 8; %%% interesting observation: one group missing segment 6, the other missing segment 4
numClusters = 12;
% numClusters = 12; %%% number of segments used in the HDM paper
% SignVectors = sqrtInvD*U(:,2:25); %%% number of coordinates used in the HDM paper
SignVectors = sqrtInvD*U(:,2:13);
sumd = Inf;
cback = 0;
for j=1:20
    [t_idx, ~, t_sumd] = kmeans(SignVectors,numClusters,'MaxIter',200);
    if sum(t_sumd) < sum(sumd)
        idx = t_idx;
        sumd = t_sumd;
    end
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('kmeans clustering iteration %4d/20 done.\n',j);
end
%%% Some idx might be +/-Inf, since sqrtInvD might contain +/-Inf
%%% Better to insert a piece of code here assigning a non-nan label to
%%% +/-Inf points in idx
[InfIdx,~] = find(isinf(SignVectors));
InfIdx = unique(InfIdx);
for j=1:length(InfIdx)
    IdxJ = find(nVListCumsum>=InfIdx(j),1);
    NamesJ = Names{IdxJ};
    load([sample_path NamesJ '.mat']);
    ValidVList = 1:G.nV;
    IdxOnG = idx(NamesDelimit(IdxJ,1):NamesDelimit(IdxJ,2));
    ValidVList(IdxOnG == idx(InfIdx(j))) = [];
    tmpDistMatrix = pdist2(G.V(:,InfIdx(j)-NamesDelimit(IdxJ,1)+1)',G.V(:,ValidVList)');
    [~,minInd] = min(tmpDistMatrix);
    idx(InfIdx(j)) = idx(ValidVList(minInd)+NamesDelimit(IdxJ,1)-1);
end
ViewBundleFunc(Names,idx,options);

keyboard

%% align every mesh to the first one, and recompute normals
for j = 1:GroupSize
    meshList{j}.V = options.R{j} * meshList{j}.V;
    [meshList{j}.Nv, meshList{j}.Nf] = meshList{j}.ComputeNormal();
end

save('HBDM.mat','HBDM','idx','numClusters','meshList','taxa_code','GroupSize','GroupDelimit','PerGroupDelimit','NamesDelimit','PerGroupSize','Names','NamesByGroup')

%% export segments on each tooth
numSegments = zeros(1,GroupSize);
segNumVertices = zeros(GroupSize,numClusters);
segAreas = zeros(GroupSize,numClusters);
segPath = './segments/';
segCentroids = cell(GroupSize,numClusters);
segCentroidsFlat = cell(GroupSize,numClusters);
segColorList = [83,123,113;...
                31,120,180;...
                178,223,138;...
                51,160,44;...
                251,154,153;...
                227,26,28;...
                253,191,111;...
                255,127,0;...
                202,178,214;...
                106,61,154;...
                51,51,51;...
                177,89,40]/255;
close(gcf);
for j = 1:GroupSize
    individualSegPath = [segPath Names{j} '/'];
    labelsOnMesh = idx(NamesDelimit(j,1):NamesDelimit(j,2));
    uniqueLabelsOnMesh = unique(labelsOnMesh);
    numSegments(j) = length(uniqueLabelsOnMesh);
    
    disp(uniqueLabelsOnMesh');
    
    touch(individualSegPath);
    avgSegVertNormal = zeros(numClusters,3); %%% normal always taken in R^3
    avgSegFaceNormal = zeros(numClusters,3); %%% normal always taken in R^3
    for k = 1:numClusters
        Gseg = Mesh(meshList{j});
        dInds = find(labelsOnMesh ~= uniqueLabelsOnMesh(k));
        Gseg.DeleteVertex(dInds);
        segNumVertices(j,k) = Gseg.nV;
        segAreas(j,k) = Gseg.ComputeSurfaceArea();
        segCentroids{j,k} = mean(Gseg.V,2);
        [~,closestVIdxOnSegment] = min(pdist2(segCentroids{j,k}', meshList{j}.V'));
        segCentroids{j,k} = meshList{j}.V(:, closestVIdxOnSegment);
        segCentroidsFlat{j,k} = meshList{j}.Aux.UniformizationV(:, closestVIdxOnSegment);
        Gseg.Write([individualSegPath sprintf('%02d', k) '.off'], 'off', []);
        
        [Gseg.Nv, Gseg.Nf] = Gseg.ComputeNormal();
        avgSegVertNormal(k,:) = mean(Gseg.Nv,2)';
        avgSegFaceNormal(k,:) = mean(Gseg.Nf,2)';
    end
    avgVertNormal = mean(meshList{j}.Nv,2)';
    avgFaceNormal = mean(meshList{j}.Nf,2)';
    
    writetable(table([avgVertNormal; avgSegVertNormal], 'RowNames', [{'total'}; cellfun(@(x) num2str(x), num2cell((1:numClusters)'), 'UniformOutput', false)]), [individualSegPath 'segVertNormals.csv'], 'WriteRowNames', true, 'WriteVariableNames', false);
    writetable(table([avgFaceNormal; avgSegFaceNormal], 'RowNames', [{'total'}; cellfun(@(x) num2str(x), num2cell((1:numClusters)'), 'UniformOutput', false)]), [individualSegPath 'segFaceNormals.csv'], 'WriteRowNames', true, 'WriteVariableNames', false);
    
    curFig = figure;
    subplot(1,2,1);
    
    %%%% draw points representing each segment and add them to legend
    slice_j = cell2mat(segCentroids(j,:));
    for kk=1:numClusters
        if kk==1
            hold on
        end
        plot3(slice_j(1,kk), slice_j(2,kk), slice_j(3,kk), 'Color', segColorList(kk,:), 'Marker', '.', 'MarkerSize', 20, 'LineStyle', 'none');
    end
    legend(cellfun(@(x) num2str(x), num2cell(1:numClusters), 'UniformOutput', false), 'Location', 'westoutside');
    
    meshList{j}.draw(struct('FaceColor', 'interp',...
        'FaceVertexCData', segColorList(labelsOnMesh,:),...
        'CDataMapping','scaled',...
        'EdgeColor', 'none',...
        'FaceAlpha', 1,...
        'AmbientStrength', 0.3,...
        'SpecularStrength', 0.0));
    hold on
    camlight('headlight');
    camlight(180,0);
    
    subplot(1,2,2);
    %%%% draw points representing each segment and add them to legend
    slice_j = cell2mat(segCentroidsFlat(j,:));
    for kk=1:numClusters
        if kk==1
            hold on
        end
        plot3(slice_j(1,kk), slice_j(2,kk), slice_j(3,kk), 'Color', segColorList(kk,:), 'Marker', '.', 'MarkerSize', 20, 'LineStyle', 'none');
    end
    legend(cellfun(@(x) num2str(x), num2cell(1:numClusters), 'UniformOutput', false), 'Location', 'eastoutside');
    
    meshList{j}.draw('flat',...
        struct('FaceColor', 'interp',...
        'FaceVertexCData', segColorList(labelsOnMesh,:),...
        'CDataMapping','scaled',...
        'EdgeColor', 'none',...
        'FaceAlpha', 1,...
        'AmbientStrength', 0.3,...
        'SpecularStrength', 0.0));
    hold on
    camlight('headlight');
    camlight(180,0);
    
    savefig(curFig, [individualSegPath 'segmentIDs.fig']);
    close(curFig);
end
normalizedSegAreas = segAreas./repmat(sum(segAreas, 2), 1, numClusters);

writetable(table([1:numClusters; segNumVertices], 'RowNames', [{'Seg ID'} Names]), [segPath 'segNumVertices.csv'], 'WriteRowNames', true, 'WriteVariableNames', false);
writetable(table([1:numClusters; segAreas], 'RowNames', [{'Seg ID'} Names]), [segPath 'segAreas.csv'], 'WriteRowNames', true, 'WriteVariableNames', false);
writetable(table([1:numClusters; normalizedSegAreas], 'RowNames', [{'Seg ID'} Names]), [segPath 'normalizedSegAreas.csv'], 'WriteRowNames', true, 'WriteVariableNames', false);
