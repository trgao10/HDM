function [rslt] = HypoellipticDiffusionMap(hdm)
%
%  Hypoelliptic Diffusion Map v0.1
%
% INPUT:
%   hdm.data	    : pxn matrix which represents n data points in R^p
%   hdm.BEpsilon	: bandwidth on the base
%   hdm.FEpsilon    : bandwidth on the fibre
%   hdm.BNN         : number of nearest neighbors on the base
%   hdm.FNN         : number of nearest neighbors on each fibre
%   hdm.T	        : diffusion time
%   hdm.delta	    : parameter for truncation
%
% (OPTION)
%   hdm.debug       : debug mode
%   hdm.embedmaxdim : largest allowable dimension for embedding (embedding
%                     is mandatory, according to truncation; this option
%                     sets the upper bound for allowable embeddings)
%
% OUTPUT:
%   rslt.embeddim   : the truncated HDM
%
% DEPENDENCE:
%   LocalPCA.m, NearestNeighbors/
%
% last modified by Tingran Gao (trgao10@math.duke.edu) 2014-12-12
%

%%% indispensible input fields
if ~isfield(hdm, 'BEpsilon')
    error('missing field "BEpsilon" from input');
end
if ~isfield(hdm, 'FEpsilon')
    error('missing field "FEpsilon" from input');
end
if ~isfield(hdm, 'BNN')
    error('missing field "BNN" from input');
end
if ~isfield(hdm, 'FNN')
    error('missing field "FNN" from input');
end

hdm.alpha         = getoptions(hdm, 'alpha', 1);
hdm.delta         = getoptions(hdm, 'delta', 0.9);
hdm.T             = getoptions(hdm, 'T', 2);
hdm.FSampleType   = getoptions(hdm, 'FSampleType', 'uniform');
hdm.debug         = getoptions(hdm, 'debug', 0);
hdm.embedmaxdim   = getoptions(hdm, 'embedmaxdim', 100);

eigopt.isreal = 1;
eigopt.issym  = 1;
eigopt.maxit  = 3000;
eigopt.disp   = 0;

%%% number of samples on base and fibre
NB = size(hdm.data,2);
NF = hdm.NF;
% NF = hdm.BNN; %%% subject to change
FNN = hdm.FNN; %%% BNN will be corrected

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Step 1: Nearest Neighbor Search
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if hdm.debug==1
    fprintf('(DEBUG:HDM) Step 1: NN search and Data preparation.\n');
end

atria = nn_prepare(hdm.data');
[BaseIndex,BaseWeight] = nn_search(hdm.data',atria,(1:NB)',hdm.BNN+1,-1,0.0);
%%% BaseWeight was BaseDistance before.
%%% We'll fill in exp(-BaseWeight.^2/hdm.BEpsilon) later in space

if hdm.debug==1
    fprintf(['(DEBUG:HDM) BNN=',num2str(hdm.BNN),'\n']);
    fprintf(['(DEBUG:HDM) minimal farthest distance=',num2str(min(BaseWeight(:,end))),'\n']);
    fprintf(['(DEBUG:HDM) maximal farthest distance=',num2str(max(BaseWeight(:,end))),'\n']);
    fprintf(['(DEBUG:HDM) median farthest distance=',num2str(median(BaseWeight(:,end))),'\n']);
    fprintf(['(DEBUG:HDM) 1.5*sqrt(min farthest distance)=',num2str(1.5*sqrt(min(BaseWeight(:,end)))),'.\n']);
    fprintf(['(DEBUG:HDM) 1.5*sqrt(max farthest distance)=',num2str(1.5*sqrt(max(BaseWeight(:,end)))),'.\n']);
    fprintf(['(DEBUG:HDM) 1.5*sqrt(median farthest distance)=',num2str(1.5*sqrt(median(BaseWeight(:,end)))),'.\n']);
end

%%% patchno is set to convert the NN information to the \sqrt{h} information
patchno = hdm.BNN*ones(1,NB);

if hdm.debug==1
    fprintf('(DEBUG:HDM) neighbor points with kernel value less than exp(-5*1.5^2) are trimmed.\n');
end

for ii=1:NB
    patchno(ii) = sum(BaseWeight(ii,:) <= 1.5*sqrt(hdm.BEpsilon));
    BaseWeight(ii, BaseWeight(ii,:) > 1.5*sqrt(hdm.BEpsilon)) = Inf;
end

%%% truncate distance and index by patchno
BNN = max(patchno)-1;
BaseWeight = BaseWeight(:, 1:(BNN+1));
BaseIndex = BaseIndex(:, 1:(BNN+1));

if hdm.debug==1
    if (quantile(patchno,0.9) == hdm.BNN)
        %%% it means the NN is not big enough so that the decay of the
        %%% kernel won't be fast enough for the error to be small
        fprintf('(DEBUG:HDM:WARNING) BNN should be chosen larger\n');
    end
    
    fprintf(['(DEBUG:HDM) the number of points with distance less then\n']);
    fprintf(['           1.5*sqrt(BEpsilon)=' num2str(1.5*sqrt(hdm.BEpsilon)) ...
        ' is (min,max,median) = (' num2str(min(patchno)) ',' ...
        num2str(max(patchno)) ',' num2str(median(patchno)) ')\n']);
    fprintf(['(DEBUG:HDM) set BNN to be ',num2str(max(patchno)),'\n']);
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Step 2: Local PCA
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if hdm.debug==1
    fprintf('(DEBUG:HDM) Step 2: find a basis for each tangent plane by PCA\n');
end

lpca.data      = hdm.data;
%%% (CAUTION) theoretically, it should be HDM.epsilon^((D+4)/(D+1))
%%% BUT we don't know the dimension D yet!
%%% throw a random guess instead...
guess = 0.8;
lpca.epsilonpca= hdm.BEpsilon*guess;
lpca.NN        = max(patchno);
lpca.index     = BaseIndex;
lpca.distance  = BaseWeight;
lpca.patchno   = patchno;
%%% lpca.KN is the threshold number gamma chosen by the user in local PCA.
%%% by default choose 0.9
lpca.KN        = 0.9;
lpca.debug     = 1;

[pcaBASIS, D]  = LocalPCA(lpca); %%% D is the estimated dimension
clear lpca patchno

%%% view pcaBASIS and samples on the unit circle
% figure;
% scatter3(hdm.data(1,:),hdm.data(2,:),hdm.data(3,:),1,'b','filled');
% axis equal;
% hold on
% CheckIdx = 1235;
% scatter3(hdm.data(1,CheckIdx),hdm.data(2,CheckIdx),hdm.data(3,CheckIdx),10,'r','filled');
% tangent_basis = pcaBASIS(:,:,CheckIdx);
% basis = pcaBASIS(:,:,CheckIdx)+repmat(hdm.data(:,CheckIdx),1,D);
% scatter3(basis(1,:),basis(2,:),basis(3,:),10,'g','filled');
% line([hdm.data(1,CheckIdx),basis(1,1)],[hdm.data(2,CheckIdx),basis(2,1)],[hdm.data(3,CheckIdx),basis(3,1)],'color','g');
% line([hdm.data(1,CheckIdx),basis(1,2)],[hdm.data(2,CheckIdx),basis(2,2)],[hdm.data(3,CheckIdx),basis(3,2)],'color','g');
% unit_circle = hdm.data(:,BaseIndex(CheckIdx,2:end));
% unit_circle = unit_circle - repmat(hdm.data(:,CheckIdx),1,size(unit_circle,2));
% unit_circle = tangent_basis*tangent_basis'*unit_circle;
% unit_circle = unit_circle./repmat(sqrt(sum(unit_circle.^2)),3,1)+repmat(hdm.data(:,CheckIdx),1,size(unit_circle,2));
% scatter3(unit_circle(1,:),unit_circle(2,:),unit_circle(3,:),10,'k','filled');
% keyboard

%%% exclude the point itself from its nearest neighbors
BaseWeight = BaseWeight(:, 2:end);
BaseIndex = BaseIndex(:, 2:end);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Step 3: Resolving Reflection
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if hdm.debug==1
    fprintf('(DEBUG:HDM) Step 3: deal with the reflection effect\t\t');
end

IdxI = zeros(NB*BNN,1);
IdxJ = zeros(NB*BNN,1);
Vals = zeros(NB*BNN,1);

cback = 0;
for ii=1:NB
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('%4d',ii);
    
    IdxI((ii-1)*BNN+1) = ii;
    IdxJ((ii-1)*BNN+1) = ii;
    Vals((ii-1)*BNN+1) = 1;
    
    for kk = 2:BNN
        jj = BaseIndex(ii,kk);	%% ii and jj are indices of points
        
        IdxI((ii-1)*BNN+kk) = ii;
        IdxJ((ii-1)*BNN+kk) = jj;
        
        Ai = pcaBASIS(:,1:D,ii);
        Aj = pcaBASIS(:,1:D,jj);
        H  = Ai'*Aj; %%% coordinate transform from Aj to Ai
        [U,~,V] = svd(H);
        
        %%% U*V' approximates the parallel transport from Aj to Ai
        if det(U*V')<0
            Vals((ii-1)*BNN+kk) = -1;
        else
            Vals((ii-1)*BNN+kk) = 1;
        end
    end
end
clear U V lambda

REFLECTION = sparse(IdxI,IdxJ,Vals);
[UR, ~] = eigs(sparse(REFLECTION),2,'lm',eigopt);

%%% make all frames coincide with SO(D)
refIdx = find(UR(ii,1)<0);
pcaBASIS(:,1,refIdx) = -pcaBASIS(:,1,refIdx);

clear UR lambdaR IdxI IdxJ Vals REFLECTION

if hdm.debug==1
    fprintf('\n');
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% symmetrize the graph
if hdm.debug==1
    fprintf('(DEBUG:HDM) Step 3'': symmetrize the base heat kernel\t\t');
end

IdxI = repmat((1:NB)',BNN,1);
IdxJ = BaseIndex(:);
tmp = sparse(IdxI,IdxJ,exp(-BaseWeight(:).^2/hdm.BEpsilon));
tmp = min(tmp,tmp');

for ii=1:NB
    BaseWeight(ii,:) = tmp(ii,BaseIndex(ii,:));
end

clear IdxI IdxJ tmp

if hdm.debug==1
    fprintf('\n');
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Step 4: Construct Heat Kernel Matrix
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% This step is still subject to improvement for now.
%%% Should be able to conbime this construction with the symmetrization!!
%%% Idea is to check BaseWeight:
%%% 1. Find the number of non-zero elements in BaseWeight. This number is
%%%    also the number of non-zero blocks in full-sized Hc. Divide this
%%%    number by 2, and only maintain half of the blocks;
%%% 2. Allocate corresponding spaces for IdxI, IdxJ, Vals (all vectors),
%%%    and maintain the alpha-normalization vector D;
%%% 3. Find and symmetrize each block on the fly, store indices in IdxI,
%%%    IdxJ, and values in Vals. Remember to only store blocks with
%%%    BaseWeight(ii,kk)>0 and jj>ii, where jj = BaseIndex(ii,kk), and in
%%%    the meanwhile update two segments in D (one for row ii, the other
%%%    for row jj);
%%% 4. It is then possible to alpha-normalize Vals in place, like before;
%%% 5. Replace D with the row sum of alpna-normalized Vals (slow but in
%%%    place), keep track of two segments at a time. No need to divide H by
%%%    D after this step -- we'll solve a generalized eigenvalue problem so
%%%    as to avoid numerical issues (H*V=sparse(1:NB*NF,1:NB*NF,D)*V*La);
%%% 6. H = sparse(IdxI, IdxJ, Vals, NB*NF, NB*NF, NB*NF*BNN*FNN);
%%%    clear IdxI IdxJ Vals
%%%    H = (H+H')/2;
%%% 7. [~,lambda] = eigs(H, D, hdm.embedmaxdim, 'lm', eigopt);
%%%

if hdm.debug==1
    fprintf('(DEBUG:HDM) Step 4: Construct Heat Kernel\t\t');
end

if strcmpi(hdm.FSampleType,'uniform')
    %%% The tangent space is of dimension D, thus circle_sample should be
    %%% NF uniformly distributed points on S^(D-1).
    %%% For now we choose BNN points on each fibre (so as to be consistent
    %%% with the other sampling strategy on fibres, see below).
    if D==2
        unit_sphere = linspace(0,2*pi,NF+1);
        unit_sphere(end) = [];
        unit_sphere = [cos(unit_sphere);sin(unit_sphere)];
    else
        unit_sphere = randn(D, NF);
        unit_sphere = unit_sphere * diag(1./sqrt(sum(unit_sphere.^2)));
    end
    unit_sphere_kdtree = kdtree_build(unit_sphere');
end

%%% find number of non-zero elements in BaseWeight
nnzBlocks = sum(BaseWeight(:)>0)/2;

%%% pre-allocate space
Vals = zeros(NF*FNN*nnzBlocks,1);
IdxI = zeros(NF*FNN*nnzBlocks,1);
IdxJ = zeros(NF*FNN*nnzBlocks,1);
Da = zeros(NB*NF,1);

count = 0;
cback = 0;
for ii=1:NB
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('%4d',ii);
    
    Aii = pcaBASIS(:,1:D,ii);
    
    for kk=1:BNN
        jj = BaseIndex(ii,kk);
        if ((jj > ii) && (BaseWeight(ii,kk)>0))
            %%% if BaseWeight(ii,kk)==0, we know that ii is not within the
            %%% BNN-nearest-neighborhood of jj
            count = count + 1;
            ss = find(BaseIndex(jj,:) == ii);
            
            blockIdxI = repmat((1:NF)',FNN,1);
            Ajj = pcaBASIS(:,1:D,jj);
            
            %%% extract block (ii,jj)
            H = Ajj'*Aii;
            [U,~,V] = svd(H);
            Pji = V*U'; %%% approximates parallel transport from j to i
            %%% Pij = U*V' approximates the parallel transport from i to j
            %%% but here we really want Pji since we set the unit sphere
            %%% at sample i fixed
            
            %%% fix the unit sphere at sample i and compute how its points
            %%% spread to the unit sphere at sample j
%             [iikkIdxJ, Dist] = kdtree_k_nearest_neighbors_tg(unit_sphere_kdtree, unit_sphere', FNN); %%% for debugging
            [iikkIdxJ, Dist] = kdtree_k_nearest_neighbors_tg(unit_sphere_kdtree, (Pji*unit_sphere)', FNN);
            iikkBlockExpand = sparse(flat(blockIdxI), flat(iikkIdxJ), flat(BaseWeight(ii,kk)*exp(-Dist.^2/hdm.FEpsilon)), NF, NF, NF*FNN);
            
%             %%% testing kdtree_k_nearest_neighbors_tg
%             idx = zeros(NF, hdm.FNN);
%             dist = zeros(NF, hdm.FNN);
%             for th = 1:NF
%                 [idx(th,:), dist(th,:)] = kdtree_k_nearest_neighbors(unit_sphere_kdtree, (Pji*unit_sphere(:,th))', hdm.FNN);
%             end
            
%             %%% plot check nearest neighbors
%             clf; axis equal; hold on
%             scatter(unit_sphere(1,:),unit_sphere(2,:),10,'r','filled');
%             scatter(Pji(1,:)*unit_sphere,Pji(2,:)*unit_sphere,10,'k','filled');
%             rng('shuffle');
%             fixIdx = randi(NF);
%             for fii=1:FNN
%                 line([unit_sphere(1,fixIdx),Pji(1,:)*unit_sphere(:,iikkIdxJ(fixIdx,fii))],[unit_sphere(2,fixIdx),Pji(2,:)*unit_sphere(:,iikkIdxJ(fixIdx,fii))],'color','b');
%             end
            
            %%% extract block (jj,ii)
%             [jjssIdxJ, Dist] = kdtree_k_nearest_neighbors_tg(unit_sphere_kdtree, unit_sphere', FNN);
            [jjssIdxJ, Dist] = kdtree_k_nearest_neighbors_tg(unit_sphere_kdtree, (Pji'*unit_sphere)', FNN);
            jjssBlockExpand = sparse(flat(blockIdxI), flat(jjssIdxJ), flat(BaseWeight(ii,kk)*exp(-Dist.^2/hdm.FEpsilon)), NF, NF, NF*FNN);
            
            %%% symmetrize iikkBlockExpand and jjssBlockExpand
            tmp = min(iikkBlockExpand, jjssBlockExpand');
            tmp = tmp.^2; %%% What is the effect of this?
            tmpc = zeros(NF,FNN);
            for rr = 1:NF
                tmpc(rr,:) = tmp(rr,iikkIdxJ(rr,:));
            end
            
            %%% only store block (ii,jj) into IdxI, IdxJ, Vals
            blockRange = ((count-1)*NF*FNN+1):count*NF*FNN;
%             IdxI(blockRange) = repmat((1:NF)+(ii-1)*NF, FNN, 1);
            IdxJ(blockRange) = flat((jj-1)*NF+iikkIdxJ);
            Vals(blockRange) = flat(tmpc);
            
            %%% update segments ii and jj of Da
            Da(((ii-1)*NF+1):ii*NF) = Da(((ii-1)*NF+1):ii*NF) + sum(tmp,2);
            Da(((jj-1)*NF+1):jj*NF) = Da(((jj-1)*NF+1):jj*NF) + sum(tmp)';
            
            clear blockIdxI iikkIdxJ jjssIdxJ iikkBlockExpand jjssBlockExpand tmp tmpc blockRange
        end
    end
end

if hdm.debug==1
    fprintf('\n');
end

clear pcaBASIS BaseWeight BaseIndex
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% alpha-normalize Hc
if (hdm.alpha > 0)
    if hdm.debug==1
        fprintf('(DEBUG:HDM) Step 4'''': Alpha-Normalizing the Heat Kernel (in place but slow)\t\t');
    end

    Da = 1./Da.^hdm.alpha;
    
    %%% alpha-normalize Hc: slow but in place
    Vals = Da(IdxI).*Vals.*Da(IdxJ);
%     for ii=1:NB*NF
%         validIdx = find(IdxJc(ii,:)>0);
%         Hc(ii,validIdx) = (Dc(ii)*Hc(ii,validIdx)).*Dc(IdxJc(ii,validIdx))';
%     end
    
    if hdm.debug==1
        fprintf('\n');
    end
end
clear Dc
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% construct full-sized Hc and normalization
% sqrtDc = sparse(1:NB*NF, 1:NB*NF, sqrt(1./sum(Hc,2)));
H = sparse(IdxI,IdxJ,Vals,NB*NF, NB*NF, NB*NF*BNN*FNN);
clear IdxI IdxJ Vals Da
H = (H+H')/2;

Da = sum(H);
denom_sqrt_D = sparse(1:NB*NF, 1:NB*NF, 1./sqrt(Da));
clear Da;
H = denom_sqrt_D*H*denom_sqrt_D;
H = (H+H')/2;

%%% symmetrize Hc
%%% not quite necessary due to the constructoin, but helps with numerics
% Hc = (Hc+Hc')/2;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Step 5: Spectral Decomposition
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if hdm.debug==1
    fprintf('(DEBUG:HDM) Step 5: Eigen-Decomposition\n');
end

%%% heavy-lifting eigen-decomposition of the heat kernel
% [~,lambda] = eigs(H, sparse(1:NB*NF,1:NB*NF,Da), hdm.embedmaxdim, 'lm', eigopt);
[~,lambda] = eigs(H, hdm.embedmaxdim, 'lm', eigopt);

% lambdaS		       = diag(lambdaS);
% [lambdaS, sortidx] = sort(lambdaS, 'descend');
% US = US(:, sortidx);	%% these are the eigenvectors for \tilde{S}
% US = sqrtDc*US;	%% these are the eigenvectors for S
% 
% rslt.US = US;
rslt.lambda = diag(lambda);
% clear UC lambdaC sparseC

if hdm.debug==1
    sdlambda = sort(diag(lambda), 'descend');
    sdlambda = 1-sdlambda;
    sdlambda = sdlambda(1:hdm.embedmaxdim);
    bar(sdlambda);
    axis([0,hdm.embedmaxdim,0,max(sdlambda)]);
    title(['NB = ' num2str(NB) ', '...
        'BNN = ' num2str(BNN) ', '...
        'BEpsilon = ' num2str(hdm.BEpsilon) ';    '...
        'NF = ' num2str(NF) ', '...
        'FNN = ' num2str(FNN) ', '...
        'FEpsilon = ' num2str(hdm.FEpsilon) ';   '...
        'alpha = ' num2str(hdm.alpha)]);
end

end
