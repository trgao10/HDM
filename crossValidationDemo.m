%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
NumCV = 10000;
mdsDim = 10;
GroupLevel = 'Genus';
GroupNames = {'Alouatta','Ateles','Brachyteles','Callicebus','Saimiri'};

%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/HDM/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
sample_path = '../cPdist/samples/HDM/';
result_path = '/media/trgao10/Work/MATLAB/ArchivedResults/HDM/cPdist/';
soften_path = [result_path 'soften/'];
% TextureCoords1Path = [result_path 'TextureCoords1/'];
% TextureCoords2Path = [result_path 'TextureCoords2/'];

%% load taxa codes
taxa_file = [data_path 'hdm_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);
ChunkSize = 25;

%% useful inline functions
ChunkIdx = @(TAXAind1,TAXAind2) ceil(((TAXAind1-1)*GroupSize+TAXAind2)/ChunkSize);

%% parse GroupNames
[~,ClTable,~] = xlsread(spreadsheet_path);
Names = {};
Labels = {};
NamesByGroup = cell(1,length(GroupNames));
TaxaByGroup = cell(1,length(GroupNames));
for j=1:length(GroupNames)
    NamesJ = ClTable(strcmpi(ClTable(1:end,strcmpi(ClTable(1,:),GroupLevel)),GroupNames{j}),1);
    Names = [Names,NamesJ{:}];
    NamesByGroup{j} = NamesJ;
    TaxaByGroup{j} = cellfun(@(x) find(strcmpi(taxa_code, x)), NamesJ);
    LabelJ = cell(size(NamesJ));
    LabelJ(:) = {GroupNames{j}};
    Labels = [Labels,LabelJ{:}];
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
colorsList = [1,0,0;0,1,0;0,0,1;0,1,1;1,1,0;1,0,1;0,0,0];
colorsList = [colorsList;colorsList/2];

%% load pre-computed MDS embeddings
load('results/Compare_DM_HBDMFiveGroupsOfTenResult.mat');
% Ydm = cmdscale(DM_dist);
% Ydm = Ydm(:,1:mdsDim);
Ydm = mdscale(DM_dist,mdsDim,'criterion','metricstress');
% Yhbdm = cmdscale(HBDM_dist);
% Yhbdm = Yhbdm(:,1:mdsDim);
Yhbdm = mdscale(HBDM_dist,mdsDim,'criterion','metricstress');

%% cross-validate DM and HBDM
DMCVs = zeros(1,NumCV);
HBDMCVs = zeros(1,NumCV);

cback = 0;
for j=1:NumCV
    if (mod(j,100) == 0)
        for cc=1:cback
            fprintf('\b');
        end
        cback = fprintf('%4d\n',j);
    end
    
    %%% use the same partition for DM and HBDM
    CValid = cvpartition(Labels,'k',5);
    
    err = zeros(CValid.NumTestSets,1);
    for i = 1:CValid.NumTestSets
        trIdx = CValid.training(i);
        teIdx = CValid.test(i);
        ytest = classify(Ydm(teIdx,:), Ydm(trIdx,:), Labels(trIdx));
        err(i) = sum(~strcmp(ytest,Labels(teIdx)'));
    end
    DMCVs(j) = sum(err)/sum(CValid.TestSize); 
    % cvErr = sum(err)/sum(CValid.TestSize);
    % disp(['Diffusion Map cvErr: ' num2str(cvErr)]);
    
    err = zeros(CValid.NumTestSets,1);
    for i = 1:CValid.NumTestSets
        trIdx = CValid.training(i);
        teIdx = CValid.test(i);
        ytest = classify(Yhbdm(teIdx,:), Yhbdm(trIdx,:), Labels(trIdx));
        err(i) = sum(~strcmp(ytest,Labels(teIdx)'));
    end
    HBDMCVs(j) = sum(err)/sum(CValid.TestSize); 
    % cvErr = sum(err)/sum(CValid.TestSize);
    % disp(['Hypoelliptic Diffusion Map cvErr: ' num2str(cvErr)]);
end

disp(['mean(DMCVs) = ' num2str(mean(DMCVs)) ', std(DMCVs) = ' num2str(std(DMCVs))]);
disp(['mean(HBDMCVs) = ' num2str(mean(HBDMCVs)) ', std(HBDMCVs) = ' num2str(std(HBDMCVs))]);

