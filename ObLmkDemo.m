%% preparation
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
GroupLevel = 'Genus';
GroupNames = {'Lemur','Adapis'};
% GroupNames = {'Arctocebus'};
% GroupNames = {'Lemur','Adapis','Arctocebus','Microcebus',...
%     'Cantius','Purgatorius','Pronothodectes','Tupaia',};

%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
sample_path = '../cPdist/samples/PNAS/';

%% load taxa codes
taxa_file = [data_path 'teeth_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);

%% options that control the diffusion eigenvector visualization
options.sample_path = sample_path;
% options.DisplayLayout = [1,4];
% options.DisplayLayout = [4,9];
options.DisplayLayout = [2,4];
options.DisplayOrient = 'Vertical';
options.boundary = 'on';
options.names = 'off';

options.type = 'full';
options.Landmarks = 'on';
options.LandmarksPath = [data_path 'landmarks_teeth.mat'];
options.MeshesPath = [data_path 'meshes/'];
options.MeshSuffix = '_sas.off';

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

%% collection rigid motions
rigid_motions = load([data_path 'cPMSTinitRms.mat']);
options.R = reshape(rigid_motions.R(TAXAinds,TAXAinds),GroupSize,GroupSize);

%% view observer landmarks
ViewObLmk(Names,options);
