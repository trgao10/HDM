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
% GroupNames = {'Alouatta'};
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
options.DisplayLayout = [2,5];
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
[U, lambda] = eigs(H, 4, 'LM', eigopt);
lambda = diag(lambda);
disp(['Eigs completed in ' num2str(toc) ' seconds']);

%%
%==========================================================================
%%% HDM (Hypoelliptic Diffusion Maps)
%==========================================================================
% sqrtInvD(isinf(sqrtInvD)) = 0;
% BundleHDM = sqrtInvD*U(:,2:end);
%==========================================================================
%%% Visualize Template
%==========================================================================
Template = sqrtInvD*U(:,2:4);
% colors = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;1,1,1];
% figure;
for j=1:length(GroupNames)
    LocalTemplate = Template(GroupDelimit(j,1):GroupDelimit(j,2),:);
    [r,~] = find(isinf(LocalTemplate));
    LocalTemplate(r,:) = [];
%     scatter3(LocalTemplate(:,1),LocalTemplate(:,2),LocalTemplate(:,3),1,colors(j,:),'filled');
%     if (j == 1)
%         hold on;
%     end
    T = LocalTemplate;
    [r,c] = find(isinf(T));
    T(r,:) = [];
    GM = Mesh('VF',T',[1;1;1]);
    options.pointCloud = 1;
    GM.Write([GroupNames{j} '.off'],'off',options);    
%     figure;scatter3(LocalTemplate(:,1),LocalTemplate(:,2),LocalTemplate(:,3),1,colors(j,:),'filled');
end


% T = Template(NamesDelimit(1,1):NamesDelimit(1,2),:);
% [r,~] = find(isinf(T));
% T(r,:) = [];
% GM = Mesh('VF',T',[1;1;1]);
% options.pointCloud = 1;
% GM.Write([Names{1} '.off'],'off',options);

% pause();
% atria = nn_prepare(Template);
% [count, neighbors] = range_search(Template, atria, 1:size(Template,1),0.01,0);
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
% SignVectors = sqrtInvD*U(:,2:25);
% % SignVectors(abs(SignVectors)<1e-10) = 0;
% % SignVectors = sign(SignVectors);
% idx = kmeans(SignVectors,10,'MaxIter',1000);
% %%% TODO: some idx might be +/-Inf, since sqrtInvD might contain +/-Inf
% %%% better insert a piece of code here assigning a non-nan label to +/-Inf
% %%% points in idx
% [InfIdx,~] = find(isinf(SignVectors));
% InfIdx = unique(InfIdx);
% for j=1:length(InfIdx)
%     IdxJ = find(nVListCumsum>=InfIdx(j),1);
%     NamesJ = Names{IdxJ};
%     load([sample_path NamesJ '.mat']);
%     ValidVList = 1:G.nV;
%     IdxOnG = idx(NamesDelimit(IdxJ,1):NamesDelimit(IdxJ,2));
%     ValidVList(IdxOnG == idx(InfIdx(j))) = [];
%     tmpDistMatrix = pdist2(G.V(:,InfIdx(j)-NamesDelimit(IdxJ,1)+1)',G.V(:,ValidVList)');
%     [~,minInd] = min(tmpDistMatrix);
%     idx(InfIdx(j)) = idx(ValidVList(minInd)+NamesDelimit(IdxJ,1)-1);
% end
% ViewBundleFunc(Names,idx,options);
% % pause();
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
% eigen_ind = 0;
% while (1)
%     eigen_ind = eigen_ind+1;
%     tic;
%     ViewBundleFunc(Names,sqrtInvD*U(:,eigen_ind),options);
%     toc;
%     pause;
% end


%==========================================================================
%%% sparse eigenvectors
%==========================================================================
% rho = 1.0;
% alpha = 1.0;
% 
% P = H^2;
% P = [P;ones(1,size(P,2))];
% b = [zeros(size(P,1),1);1];
% 
% check_size = size(P,1)-1;
% [x, ~] = basis_pursuit(P, b, rho, alpha);
% x_original = x;
% res = norm(P(1:1:check_size,:)*x_original);
% disp(['before thresholding, ||Ax-x||_2 = ' num2str(res)]);
% 
% % 
% % x = threshold_speigv(x, GroupSize);
% % x = x./sum(x);
% % 
% % res = norm(C(1:check_size,:)*x);
% % disp(['after thresholding, ||Ax-x||_2 = ' num2str(res)]);
