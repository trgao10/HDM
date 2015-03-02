%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
BaseEps = 0.03;
BNN = 7;
FibrEps = 1e-3;
MapType = 'cPMST';
FeatureFix = 'Off';

%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
sample_path = '../cPdist/samples/Teeth/';
result_path = ['/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/' MapType '/' 'FeatureFix' FeatureFix '/'];
soften_path = [result_path 'soften/'];

%% load taxa codes
taxa_file = [data_path 'teeth_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);
ChunkSize = 55;
Names = taxa_code;

%% options that control the diffusion eigenvector visualization
options.sample_path = sample_path;
options.DisplayLayout = [5,6];
options.DisplayOrient = 'Horizontal';
options.boundary = 'on';
options.names = 'off';

%% useful inline functions
ChunkIdx = @(TAXAind1,TAXAind2) ceil(((TAXAind1-1)*GroupSize+TAXAind2)/ChunkSize);

%% parse names
GroupSize = length(Names);
DiffMatrixSizeList = zeros(GroupSize,1);
TAXAinds = zeros(GroupSize,1);
NamesDelimit = zeros(GroupSize+1,2);
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,Names{j}));
    load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']);
    DiffMatrixSizeList(j) = G.nV;
    NamesDelimit(j+1,1) = NamesDelimit(j,2)+1;
    NamesDelimit(j+1,2) = NamesDelimit(j+1,1)+G.nV-1;
end
Names = taxa_code(TAXAinds); % match upper/lower cases
NamesDelimit(1,:) = [];
nVList = DiffMatrixSizeList;
nVListCumsum = cumsum(nVList);

%% collection rigid motions
rigid_motions = load([data_path 'cPMSTinitRms.mat']);
options.R = reshape(rigid_motions.R(TAXAinds,TAXAinds),GroupSize,GroupSize);

%% process base diffusion
load([result_path MapType 'DistMatrix.mat']);
if strcmpi(MapType,'cP')
    eval(['BaseDistMatrix = ' MapType 'DistMatrix(TAXAinds,TAXAinds);']);
else
    eval(['BaseDistMatrix = ImprDistMatrix(TAXAinds,TAXAinds);']);
end

%%% only connect BNN-nearest-neighbors
[sDists,rowNNs] = sort(BaseDistMatrix,2);
sDists = sDists(:,2:(1+BNN));
rowNNs = rowNNs(:,2:(1+BNN));
BaseWeights = sparse(repmat((1:GroupSize)',1,BNN),rowNNs,sDists,GroupSize,GroupSize);
BaseWeights = min(BaseWeights, BaseWeights');
for j=1:GroupSize
    sDists(j,:) = BaseWeights(j,rowNNs(j,:));
end
sDists = exp(-sDists.^2/BaseEps);

%% build diffusion kernel matrix
DiffMatrixSize = sum(DiffMatrixSizeList);
DiffMatrixSizeList = cumsum(DiffMatrixSizeList);
DiffMatrixSizeList = [0; DiffMatrixSizeList];
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
sqrtD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,sqrt(sum(H)));
invD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,1./sum(H));
sqrtInvD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,1./sqrt(sum(H)));
% K = invD*H;
H = sqrtInvD*H*sqrtInvD;
H = (H+H')/2;

%%% this loop is slow but much less memory consuming
% for k=1:size(DiffMatrix,1)
%     DiffMatrix(k,:) = sqrtInvD(k)*DiffMatrix(k,:);
%     DiffMatrix(:,k) = sqrtInvD(k)*DiffMatrix(:,k);
% end

eigopt.isreal = 1;
eigopt.issym = 1;
eigopt.maxit = 5000;
eigopt.disp = 0;
tic;
[U, lambda] = eigs(H, 101, 'LM', eigopt);
lambda = diag(lambda);
disp(['Eigs completed in ' num2str(toc) ' seconds']);

%==========================================================================
%%% HDBM (Hypoelliptic Diffusion Base Maps)
%==========================================================================
sqrtInvD(isinf(sqrtInvD)) = 0;
BundleHDM = sqrtInvD*U(:,2:end);
% BundleHDM = sqrtInvD*U(:,2:end)*sparse(1:(size(U,2)-1), 1:(size(U,2)-1), sqrt(lambda(2:end)));
HDBM = zeros(GroupSize, nchoosek(size(BundleHDM,2),2));
for j=1:GroupSize
    BundleHDM_Block = normc(BundleHDM(NamesDelimit(j,1):NamesDelimit(j,2),:));
    BundleHDM_Block = BundleHDM_Block*sparse(1:(size(U,2)-1), 1:(size(U,2)-1), sqrt(lambda(2:end)));
    HDBM(j,:) = pdist(BundleHDM_Block', @(x,y) y*x');
%     HDBM(j,:) = pdist(BundleHDM(NamesDelimit(j,1):NamesDelimit(j,2),:)',@(x,y) y*x');
end
HDBM_dist = pdist(HDBM);
[Y,stress] = mdscale(HDBM_dist,3,'criterion','metricstress');

save('PNAS_HDBM.mat', 'lambda', 'U', 'HDBM', 'Names', 'Y');

