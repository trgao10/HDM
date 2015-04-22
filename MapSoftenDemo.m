%%% preparation
% clear all;
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%%% setup names
Names = {'B03', 'D09'};

%%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
sample_path = '../cPdist/samples/PNAS/';
result_path = '/media/trgao10/Work/MATLAB/ArchivedResults/PNAS/cPDist/';
TextureCoords1Path = [result_path 'TextureCoords1/'];
TextureCoords2Path = [result_path 'TextureCoords2/'];

%%% load taxa codes
taxa_file = [data_path 'teeth_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);
ChunkSize = 55;

%%% useful inline functions
ChunkIdx = @(TAXAind1,TAXAind2) ceil(((TAXAind1-1)*GroupSize+TAXAind2)/ChunkSize);

%%% translate names to TAXAinds
TAXAind1 = find(strcmpi(taxa_code, Names{1}));
TAXAind2 = find(strcmpi(taxa_code, Names{2}));

%%% load texture coordinates
load([TextureCoords1Path 'TextureCoords1_mat_' num2str(ChunkIdx(TAXAind1,TAXAind2)) '.mat']);
load([TextureCoords2Path 'TextureCoords2_mat_' num2str(ChunkIdx(TAXAind1,TAXAind2)) '.mat']);
TextureCoords1 = TextureCoords1Matrix{TAXAind1,TAXAind2};
TextureCoords2 = TextureCoords2Matrix{TAXAind1,TAXAind2};

%%% load mesh faces and vertices
G1 = load([sample_path taxa_code{TAXAind1} '.mat']); G1 = G1.G;
G2 = load([sample_path taxa_code{TAXAind2} '.mat']); G2 = G2.G;

%%% MapSoften: BC (BaryCentric)
[Transplan12,V1onV2] = MapSoftenBC(TextureCoords1,TextureCoords2,G2.F,G1.V,G2.V);
% CheckInd = 1000;
% figure;G2.draw();hold on;
% scatter3(V1onV2(1,CheckInd),V1onV2(2,CheckInd),V1onV2(3,CheckInd),10,'g','filled');
% [SortedWeights,SortInds] = sort(Transplan12(CheckInd,:));
% scatter3(G2.V(1,SortInds(end)),G2.V(2,SortInds(end)),G2.V(3,SortInds(end)),20,'r','filled');
% disp(['Red point weights ' num2str(SortedWeights(end))]);
% scatter3(G2.V(1,SortInds(end-1)),G2.V(2,SortInds(end-1)),G2.V(3,SortInds(end-1)),20,'y','filled');
% disp(['Yellow point weights ' num2str(SortedWeights(end-1))]);
% scatter3(G2.V(1,SortInds(end-2)),G2.V(2,SortInds(end-2)),G2.V(3,SortInds(end-2)),20,'b','filled');
% disp(['Blue point weights ' num2str(SortedWeights(end-2))]);

[Transplan21,V2onV1] = MapSoftenBC(TextureCoords2,TextureCoords1,G1.F,G2.V,G1.V);
% CheckInd = 100;
% figure;G1.draw();hold on;
% scatter3(V2onV1(1,CheckInd),V2onV1(2,CheckInd),V2onV1(3,CheckInd),10,'g','filled');
% [SortedWeights,SortInds] = sort(Transplan21(CheckInd,:));
% scatter3(G1.V(1,SortInds(end)),G1.V(2,SortInds(end)),G1.V(3,SortInds(end)),20,'r','filled');
% disp(['Red point weights ' num2str(SortedWeights(end))]);
% scatter3(G1.V(1,SortInds(end-1)),G1.V(2,SortInds(end-1)),G1.V(3,SortInds(end-1)),20,'y','filled');
% disp(['Yellow point weights ' num2str(SortedWeights(end-1))]);
% scatter3(G1.V(1,SortInds(end-2)),G1.V(2,SortInds(end-2)),G1.V(3,SortInds(end-2)),20,'b','filled');
% disp(['Blue point weights ' num2str(SortedWeights(end-2))]);

%%% MapSoften: kernel
augParam = 1.5;
[Transplan12,Kernel12,AugKernel12,V1onV2] = MapSoftenKernel(TextureCoords1,TextureCoords2,G2.F,G1.V,G2.V,'auto',augParam);
CheckInd = 466;
figure;G2.draw();hold on;
scatter3(V1onV2(1,CheckInd),V1onV2(2,CheckInd),V1onV2(3,CheckInd),10,'g','filled');
SupportInds = find(AugKernel12(CheckInd,:));
scatter3(G2.V(1,SupportInds),G2.V(2,SupportInds),G2.V(3,SupportInds),20,'b','filled');

[Transplan21,Kernel21,AugKernel21,V2onV1] = MapSoftenKernel(TextureCoords2,TextureCoords1,G1.F,G2.V,G1.V,'auto',augParam);
% CheckInd = 235;
% figure;G1.draw();hold on;
% scatter3(V2onV1(1,CheckInd),V2onV1(2,CheckInd),V2onV1(3,CheckInd),10,'g','filled');
% SupportInds = find(AugKernel21(CheckInd,:));
% scatter3(G1.V(1,SupportInds),G1.V(2,SupportInds),G1.V(3,SupportInds),20,'b','filled');

%%%% if augParam is large, there are more points included, so one may
%%%% prefer to use SymAugKernel12 = min(AugKernel12,AugKernel21');
%%%% if augParam is small, there are fewer points included, so one may
%%%% prefer to use SymAugKernel12 = max(AugKernel12,AugKernel21');
SymAugKernel12 = max(AugKernel12,AugKernel21');
SymAugKernel21 = SymAugKernel12';
% CheckInd = 235;
figure;G2.draw();hold on;
scatter3(V1onV2(1,CheckInd),V1onV2(2,CheckInd),V1onV2(3,CheckInd),50,'g','filled');
SupportInds = find(SymAugKernel12(CheckInd,:));
scatter3(G2.V(1,SupportInds),G2.V(2,SupportInds),G2.V(3,SupportInds),50,'b','filled');



