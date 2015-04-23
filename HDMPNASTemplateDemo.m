%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
BaseEps = 0.03;
BNN = 5;
FibrEps = 1e-3;
MapType = 'cPMST';
FeatureFix = 'Off';
GroupLevel = 'Genus';
% GroupNames = {'Purgatorius','Pronothodectes','Tupaia','Lemur',...
%     'Microcebus','Cantius','Arctocebus','Adapis','Lepilemur',...
%     'Eosimias','Cynocephalus','Leptacodon','Nycticebus'};
% GroupNames = {'Euprimates','Primates','Dermoptera','Scandentia','Incertae sedis'};
% GroupNames = {'Purgatorius'};
% GroupNames = {'Purgatorius','Pronothodectes'};
% GroupNames = {'Purgatorius','Pronothodectes','Tupaia','Lemur'};
% GroupNames = {'Purgatorius','Pronothodectes','Tupaia','Lemur',...
%     'Microcebus','Cantius','Arctocebus','Adapis','Lepilemur',...
%     'Eosimias','Cynocephalus'};
GroupNames = {'Purgatorius','Pronothodectes','Tupaia','Lemur',...
    'Microcebus','Cantius','Arctocebus','Adapis','Lepilemur',...
    'Eosimias','Cynocephalus','Varecia'};
% GroupNames = {'Donrussellia','Cheirogaleus','Avahi','Eulemur',...
%     'Hapalemur','Loris','Nycticebus','Leptacodon'};
% GroupNames = {'Tupaia','Galago'};
% GroupNames = {'Purgatorius','Pronothodectes','Tupaia','Microcebus','Lemur','Varecia'};
% GroupNames = {'Purgatorius','Pronothodectes','Tupaia','Microcebus'};

%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
sample_path = '../cPdist/samples/PNAS/';
result_path = ['/media/trgao10/Work/MATLAB/ArchivedResults/PNAS/' MapType '/' 'FeatureFix' FeatureFix '/'];
soften_path = [result_path 'soften/'];

viz_path = './results/TwelveLemurs/';
touch(viz_path);
% TextureCoords1Path = [result_path 'TextureCoords1/'];
% TextureCoords2Path = [result_path 'TextureCoords2/'];

%% load taxa codes
taxa_file = [data_path 'teeth_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);
ChunkSize = 55;

%% options that control the diffusion eigenvector visualization
options.sample_path = sample_path;
options.DisplayLayout = [4,12];
options.DisplayOrient = 'Vertical';
options.boundary = 'on';
options.names = 'off';

%% useful inline functions
ChunkIdx = @(TAXAind1,TAXAind2) ceil(((TAXAind1-1)*GroupSize+TAXAind2)/ChunkSize);

%% parse GroupNames
[~,ClTable,~] = xlsread(spreadsheet_path);
Names = {};
NamesByGroup = cell(1,length(GroupNames));
for j=1:length(GroupNames)
    NamesJ = ClTable(strcmpi(ClTable(1:end,strcmpi(ClTable(1,:),GroupLevel)),GroupNames{j}),1);
    Names = [Names,NamesJ{:}];
    NamesByGroup{j} = NamesJ;
end

GroupSize = length(Names);
DiffMatrixSizeList = zeros(GroupSize,1);
TAXAinds = zeros(GroupSize,1);
NamesDelimit = zeros(GroupSize+1,2);
MeshList = cell(GroupSize,1);
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,Names{j}));
    load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']);
    G.nV = size(G.V,2);
    G.nF = size(G.F,2);
    MeshList{j} = G;
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
BaseDistMatrix = BaseDistMatrix-diag(diag(BaseDistMatrix));

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
% BaseWeights = exp(-BaseWeights.^2/BaseEps);
% BaseWeights = BaseWeights - diag(diag(BaseWeights));
% % BaseWeights = diag(1./sum(BaseWeights,2))*BaseWeights;
% % keyboard

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
sqrtD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,sqrt(sum(H)));
invD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,1./sum(H));
sqrtInvD = sparse(1:DiffMatrixSize,1:DiffMatrixSize,1./sqrt(sum(H)));
% K = invD*H;
H = sqrtInvD*H*sqrtInvD;
H = (H+H')/2;

eigopt = struct('isreal',1,'issym',1,'maxit',5000,'disp',0);
tic;
[U, lambda] = eigs(H, 4, 'LM', eigopt);
lambda = diag(lambda);
disp(['Eigs completed in ' num2str(toc) ' seconds']);

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
    GM.Write([viz_path GroupNames{j} '.off'],'off',options);    
%     figure;scatter3(LocalTemplate(:,1),LocalTemplate(:,2),LocalTemplate(:,3),1,colors(j,:),'filled');
end

T = Template;
[r,~] = find(isinf(T));
T(r,:) = [];
GM = Mesh('VF',T',[1;1;1]);
options.pointCloud = 1;
GM.Write([viz_path 'preTemplate.off'],'off',options);

disp('RIMLS in MeshLab!');
keyboard

%%
%%% manually MLS template in MeshLab
TemplatePCloud = textread([viz_path 'Template.off']);
TemplatePCloud = TemplatePCloud(:,1:3);

ShapeVar = [];
for j=1:GroupSize
    LocalShapeVar = zeros(MeshList{j}.nV,1);
    LocalTemplate = Template(NamesDelimit(j,1):NamesDelimit(j,2),:)';
    cback = 0;
    for k=1:MeshList{j}.nV
        LocalShapeVar(k) = norm(LocalTemplate(:,k)-vanillaMLS(LocalTemplate(:,k), TemplatePCloud'));
        
        for cc=1:cback
            fprintf('\b');
        end
        cback = fprintf(['%4d/' num2str(MeshList{j}.nV) ' done.\n'], k);
    end
    ShapeVar = [ShapeVar;LocalShapeVar];    
%     invalidValues = isinf(LocalShapeVar) & isnan(LocalShapeVar);
%     LocalShapeVar(invalidValues) = mean(LocalShapeVar(~invalidValues));
%     Color = LocalShapeVar;
    % Color = (MeshList{1}.A*ShapeVar')';
    % Color = log(Color);
    % Color = log(Color)-min(log(Color));
    % Color = log(ShapeVar)-min(log(ShapeVar));
    % Color = (Color-mean(Color))/std(Color);
%     Color = clamp(Color, mean(Color)-2*std(Color), mean(Color)+2*std(Color));
%     ShapeVar = [ShapeVar;LocalShapeVar];
%     MeshList{1}.ViewFunctionOnMesh(Color', struct('mode','native'));
end

invalidValues = (isinf(ShapeVar) | isnan(ShapeVar));
ShapeVar(invalidValues) = mean(ShapeVar(~invalidValues));
% Color = ShapeVar;
Color = clamp(ShapeVar, mean(ShapeVar)-2*std(ShapeVar), mean(ShapeVar)+2*std(ShapeVar));

ViewBundleFunc(Names, Color, options);

save('TemplateVar.mat','Names','Color','options','ShapeVar','MeshList');

