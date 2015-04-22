%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
BaseEps = 0.03;
BNN = 4;
FibrEps = 1e-3;
MapType = 'cP';
FeatureFix = '';
GroupLevel = 'Genus';
% GroupNames = {'Alouatta','Ateles'};
GroupNames = {'Alouatta','Ateles','Brachyteles','Callicebus','Saimiri'};

%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/HDM/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
sample_path = '../cPdist/samples/HDM/';
result_path = '/media/trgao10/Work/MATLAB/ArchivedResults/HDM/cPdist/';
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
options.DisplayOrient = 'Vertical';
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
VList = cell(GroupSize,1);
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,Names{j}));
    load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']);
    DiffMatrixSizeList(j) = G.nV;
    VList{j} = G.V;
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
colorsList = [228,26,28;0,126,204;77,175,74;152,78,163;255,127,0]/255;

%% collection rigid motions
rigid_motions = load([data_path 'cPMSTinitRms.mat']);
options.R = rigid_motions.R(TAXAinds,TAXAinds);

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

sqrtInvD(isinf(sqrtInvD)) = 0;
BundleHDM = sqrtInvD*U(:,2:51);

CoordProj2Template = zeros(GroupSize, size(BundleHDM,2), 3);
for j=1:GroupSize
    CoordProj2Template(j,:,:) = (VList{j}*normc(BundleHDM(NamesDelimit(j,1):NamesDelimit(j,2), :)))';
end

figure('Name','CoordProj2Template');
d = 1;
for j=1:length(GroupNames)
    scatter3(CoordProj2Template(PerGroupDelimit(j,1):PerGroupDelimit(j,2),d,1),...
        CoordProj2Template(PerGroupDelimit(j,1):PerGroupDelimit(j,2),d,2),...
        CoordProj2Template(PerGroupDelimit(j,1):PerGroupDelimit(j,2),d,3),...
        30, colorsList(j,:), 'filled');
    if (j == 1)
        axis equal
        axis tight
        hold on;
    end
end
legend(GroupNames);


