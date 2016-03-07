function soften_ongrid(G1,G2,soften_mat,TAXAind1,TAXAind2,options)
%SOFTEN_ONGRID Summary of this function goes here
%   Detailed explanation goes here

G1 = load(G1); G1 = G1.G;
G2 = load(G2); G2 = G2.G;

taxa_code = options.TaxaCode;
GroupSize = length(taxa_code);
ChunkSize = options.ChunkSize;
TextureCoords1Path = options.TextureCoords1Path;
TextureCoords2Path = options.TextureCoords2Path;
FibrEps = str2num(options.FibrEps);

load(soften_mat);

%%% useful inline functions
ChunkIdx = @(TAXAind1,TAXAind2) ceil(((TAXAind1-1)*GroupSize+TAXAind2)/ChunkSize);

%%% load texture coordinates
load([TextureCoords1Path 'TextureCoords1_mat_' num2str(ChunkIdx(TAXAind1,TAXAind2)) '.mat']);
load([TextureCoords2Path 'TextureCoords2_mat_' num2str(ChunkIdx(TAXAind1,TAXAind2)) '.mat']);
TextureCoords1 = TextureCoords1Matrix{TAXAind1,TAXAind2};
TextureCoords2 = TextureCoords2Matrix{TAXAind1,TAXAind2};

%%%
tic;
disp(['Comparing ' G1.Aux.name ' vs ' G2.Aux.name '...']);
[~,~,AugKernel12,~] = MapSoftenKernel(TextureCoords1,TextureCoords2,G2.F,G1.V,G2.V,FibrEps);
[~,~,AugKernel21,~] = MapSoftenKernel(TextureCoords2,TextureCoords1,G1.F,G2.V,G1.V,FibrEps);
cPSoftMapsMatrix{TAXAind1,TAXAind2} = max(AugKernel12,AugKernel21');
%cPSoftMapsMatrix{TAXAind2,TAXAind1} = cPSoftMapsMatrix{TAXAind1,TAXAind2}';
save(soften_mat,'cPSoftMapsMatrix');
disp([G1.Aux.name ' vs ' G2.Aux.name ' done.']);
toc;

end

