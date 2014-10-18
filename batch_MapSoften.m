%%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%%% setup parameter
FibrEps = 1e-3;

%%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
sample_path = '../cPdist/samples/Teeth/';
result_path = '/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/cPDist/';
TextureCoords1Path = [result_path 'TextureCoords1/'];
TextureCoords2Path = [result_path 'TextureCoords2/'];
SoftMapsPath = [result_path 'cPSoftMapsMatrix.mat'];

%%% load taxa codes
taxa_file = [data_path 'teeth_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);
ChunkSize = 55;

%%% useful inline functions
ChunkIdx = @(TAXAind1,TAXAind2) ceil(((TAXAind1-1)*GroupSize+TAXAind2)/ChunkSize);

cPSoftMapsMatrix = cell(GroupSize,GroupSize);
%%% translate names to TAXAinds
for TAXAind1 = 1:GroupSize
    fprintf([num2str(TAXAind1) '/' num2str(GroupSize) '\n']);
    for TAXAind2 = (TAXAind1+1):GroupSize
        progressbar(TAXAind2,GroupSize,20);
        
        %%% load texture coordinates
        load([TextureCoords1Path 'TextureCoords1_mat_' num2str(ChunkIdx(TAXAind1,TAXAind2)) '.mat']);
        load([TextureCoords2Path 'TextureCoords2_mat_' num2str(ChunkIdx(TAXAind1,TAXAind2)) '.mat']);
        TextureCoords1 = TextureCoords1Matrix{TAXAind1,TAXAind2};
        TextureCoords2 = TextureCoords2Matrix{TAXAind1,TAXAind2};
        
        %%% load mesh faces and vertices
        G1 = load([sample_path taxa_code{TAXAind1} '.mat']); G1 = G1.G;
        G2 = load([sample_path taxa_code{TAXAind2} '.mat']); G2 = G2.G;        
        [~,~,AugKernel12,~] = MapSoftenKernel(TextureCoords1,TextureCoords2,G2.F,G1.V,G2.V,FibrEps);
        [~,~,AugKernel21,~] = MapSoftenKernel(TextureCoords2,TextureCoords1,G1.F,G2.V,G1.V,FibrEps);
        cPSoftMapsMatrix{TAXAind1,TAXAind2} = max(AugKernel12,AugKernel21');
        cPSoftMapsMatrix{TAXAind2,TAXAind1} = cPSoftMapsMatrix{TAXAind1,TAXAind2}';
    end
end

save(SoftMapsPath,'cPSoftMapsMatrix');

