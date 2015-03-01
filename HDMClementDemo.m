%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
BaseEps = 0.03;
FibrEps = 1e-3;
MapType = 'cPMST';

%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/Clement/';
sample_path = '../cPdist/samples/Clement/';
result_path = '../DATA/Clement/cPMST/FeatureFixOff/';
TextureCoords1Path = [result_path 'TextureCoords1/'];
TextureCoords2Path = [result_path 'TextureCoords2/'];

%% load taxa codes
taxa_file = [data_path 'clement_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);
ChunkSize = 20; % Clement Data Set

%% options that control the diffusion eigenvector visualization
options.sample_path = sample_path;
options.DisplayLayout = [5,8];
options.DisplayOrient = 'Horizontal';
options.boundary = 'on';
options.names = 'off';

%% useful inline functions
ChunkIdx = @(TAXAind1,TAXAind2) ceil(((TAXAind1-1)*GroupSize+TAXAind2)/ChunkSize);

%% parse GroupNames
Names = getFileNames(sample_path);

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
BaseWeights = exp(-BaseDistMatrix.^2/BaseEps);
BaseWeights = BaseWeights - diag(diag(BaseWeights));
% BaseWeights = diag(1./sum(BaseWeights,2))*BaseWeights;
% keyboard

%% build diffusion kernel matrix
DiffMatrixSize = sum(DiffMatrixSizeList);
DiffMatrixSizeList = cumsum(DiffMatrixSizeList);
DiffMatrixSizeList = [0; DiffMatrixSizeList];
DiffMatrixSizeList(end) = []; % treated as block shifts
DiffMatrixRowIdx = [];
DiffMatrixColIdx = [];
DiffMatrixVal = [];
for j=1:GroupSize
    G1 = load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']); G1 = G1.G;
    for k=(j+1):GroupSize
        G2 = load([sample_path taxa_code{strcmpi(taxa_code,Names{k})} '.mat']); G2 = G2.G;
        
        %%% load texture coordinates
        TAXAind1 = TAXAinds(j);
        TAXAind2 = TAXAinds(k);
        load([TextureCoords1Path 'TextureCoords1_mat_' num2str(ChunkIdx(TAXAind1,TAXAind2)) '.mat']);
        load([TextureCoords2Path 'TextureCoords2_mat_' num2str(ChunkIdx(TAXAind1,TAXAind2)) '.mat']);
        TextureCoords1 = TextureCoords1Matrix{TAXAind1,TAXAind2};
        TextureCoords2 = TextureCoords2Matrix{TAXAind1,TAXAind2};
        [~,~,AugKernel12,~] = MapSoftenKernel(TextureCoords1,TextureCoords2,G2.F,G1.V,G2.V,FibrEps);
        [~,~,AugKernel21,~] = MapSoftenKernel(TextureCoords2,TextureCoords1,G1.F,G2.V,G1.V,FibrEps);
        AugKernel12 = max(AugKernel12,AugKernel21');
        
        [rowIdx, colIdx, val] = find(AugKernel12);
        DiffMatrixRowIdx = [DiffMatrixRowIdx; rowIdx+DiffMatrixSizeList(j)];
        DiffMatrixColIdx = [DiffMatrixColIdx; colIdx+DiffMatrixSizeList(k)];
        DiffMatrixVal = [DiffMatrixVal; BaseWeights(j,k)*val];
        [rowIdx, colIdx, val] = find(AugKernel12');
        DiffMatrixRowIdx = [DiffMatrixRowIdx; rowIdx+DiffMatrixSizeList(k)];
        DiffMatrixColIdx = [DiffMatrixColIdx; colIdx+DiffMatrixSizeList(j)];
        DiffMatrixVal = [DiffMatrixVal; BaseWeights(k,j)*val];
    end
    disp([num2str(j) '/' num2str(GroupSize) ' done.']);
end

H = sparse(DiffMatrixRowIdx,DiffMatrixColIdx,DiffMatrixVal,DiffMatrixSize,DiffMatrixSize);
clear DiffMatrixColIdx DiffMatrixRowIdx DiffMatrixVal rowIdx colIdx val
clear TextureCoords1Matrix TextureCoords2Matrix

%% eigen-decomposition
sqrtD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,sqrt(sum(H)));
invD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,1./sum(H));
sqrtInvD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,1./sqrt(sum(H)));
K = invD*H;
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
[U, lambda] = eigs(H, 100, 'LM', eigopt);
lambda = diag(lambda);
disp(['Eigs completed in ' num2str(toc) ' seconds']);
% clear H

%==========================================================================
%%% convergence (to statoinary distribution) rate won't work
%%% because K has too many eivenvaluves close to 1
%==========================================================================
% stdDistribution = U(:,1)'*sqrtD;
% x = ones(1,DiffMatrixSize);
% for j = 1:10000
%     x = x*K;
% end
% max(abs(x-stdDistribution))
% ViewBundleFunc(Names,abs(x-stdDistribution)',options);

%==========================================================================
%%% Visualize Template
%==========================================================================
% Template = sqrtInvD*U(:,2:50);
% colors = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;1,1,1];
% for j=1:length(GroupNames)
%     LocalTemplate = Template(GroupDelimit(j,1):GroupDelimit(j,2),:);
%     [r,~] = find(isinf(LocalTemplate));
%     LocalTemplate(r,:) = [];
%     figure;scatter3(LocalTemplate(:,1),LocalTemplate(:,2),LocalTemplate(:,3),1,colors(j,:),'filled');
% end
% pause();
% % atria = nn_prepare(Template);
% % [count, neighbors] = range_search(Template, atria, 1:size(Template,1),0.01,0);
% figure;scatter3(Template(:,1),Template(:,2),Template(:,3),1,'k','filled');
% T3 = Template(:,1:3);
% [r,c] = find(isinf(T3));
% T3(r,:) = [];
% GM = Mesh('VF',T3',[1;1;1]);
% options.pointCloud = 1;
% GM.Write('Template3.off','off',options);
% pause();

% hold on;
% [~,Inds] = sort(count);
% MaxInds = Inds(end-10:end);
% scatter3(Template(MaxInds,1),Template(MaxInds,2),Template(MaxInds,3),50,'r','filled');
% [mu,Sigma,z] = GMM(Template,5);
% scatter3(mu(:,1),mu(:,2),mu(:,3),50,'r','filled');
% ViewBundleFunc(Names,z,options);

%==========================================================================
%%% consistent spectral clustering on each surface
%==========================================================================
SignVectors = sqrtInvD*U(:,2:15);
% SignVectors(abs(SignVectors)<1e-10) = 0;
% SignVectors = sign(SignVectors);
idx = kmeans(SignVectors,15,'MaxIter',1000);
%%% TODO: some idx might be +/-Inf, since sqrtInvD might contain +/-Inf
%%% better insert a piece of code here assigning a non-nan label to +/-Inf
%%% points in idx
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
pause();

%==========================================================================
%%% estimate intrinsic dimensionality using multiscale SVD
%==========================================================================
% EstDimOpts = struct('NumberOfTrials',1,'verbose',0,'MAXDIM',100,...
%                     'MAXAMBDIM',100,'Ptwise',true,'NetsOpts',[],...
%                     'UseSmoothedS',true, 'EnlargeScales',true );
% Template = SignVectors;
% OriginalIdx = 1:size(Template,1);
% Template(InfIdx,:) = [];
% OriginalIdx(InfIdx) = [];
% EstDimOpts.RandomizedSVDThres = min([size(Template,1),size(Template,2),100]);
% [EstDim,EstDimStats,Stats] = EstDim_MSVD(Template', EstDimOpts);
% %%% TODO: write a routine to extract submesh according to specified
% %%% vertices/faces/both
% pause();

%==========================================================================
%%% view eigenvectors
%==========================================================================
eigen_ind = 0;
while (1)
    eigen_ind = eigen_ind+1;
    tic;
    ViewBundleFunc(Names,sqrtInvD*U(:,eigen_ind),options);
    toc;
    pause;
end


