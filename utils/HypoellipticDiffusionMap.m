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
%   hdm.symmetrize  : symmetrize the graph
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
hdm.symmetrize    = getoptions(hdm, 'symmetrize', 1);
hdm.debug         = getoptions(hdm, 'debug', 0);
hdm.embedmaxdim   = getoptions(hdm, 'embedmaxdim', 100);

eigopt.isreal = 1;
eigopt.issym  = 1;
eigopt.maxit  = 3000;
eigopt.disp   = 0;

%%% number of samples on base and fibre
NB = size(hdm.data,2);
NF = hdm.BNN; %%% subject to change
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
% basis = pcaBASIS(:,:,CheckIdx)+repmat(hdm.data(:,CheckIdx),1,2);
% scatter3(basis(1,:),basis(2,:),basis(3,:),10,'g','filled');
% line([hdm.data(1,CheckIdx),basis(1,1)],[hdm.data(2,CheckIdx),basis(2,1)],[hdm.data(3,CheckIdx),basis(3,1)],'color','g');
% line([hdm.data(1,CheckIdx),basis(1,2)],[hdm.data(2,CheckIdx),basis(2,2)],[hdm.data(3,CheckIdx),basis(3,2)],'color','g');
% unit_circle = hdm.data(:,BaseIndex(CheckIdx,2:end));
% unit_circle = unit_circle - repmat(hdm.data(:,CheckIdx),1,size(unit_circle,2));
% unit_circle = tangent_basis*tangent_basis'*unit_circle;
% unit_circle = unit_circle./repmat(sqrt(sum(unit_circle.^2)),3,1)+repmat(hdm.data(:,CheckIdx),1,size(unit_circle,2));
% scatter3(unit_circle(1,:),unit_circle(2,:),unit_circle(3,:),10,'k','filled');
% pause();

%%% exclude the point itself from its neighbors
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
if hdm.symmetrize==1
    if hdm.debug==1
        fprintf('(DEBUG:HDM) Step 3'': symmetrize the graph\t\t');
    end
    
    IdxI = repmat((1:NB)',BNN,1);
    IdxJ = BaseIndex(:);
    tmp = sparse(IdxI,IdxJ,exp(-BaseWeight(:).^2/hdm.BEpsilon));
    tmp = min(tmp,tmp');
    
    for ii=1:NB
        for jj=1:BNN
            BaseWeight(ii,jj) = tmp(ii,BaseIndex(ii,jj));
        end
    end
    
    clear IdxI IdxJ tmp
else
    BaseWeight = exp(-BaseWeight.^2/hdm.BEpsilon);
end

if hdm.debug==1
    fprintf('\n');
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Step 4: Construct Heat Kernel Matrix
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

%%% Hc is like H but compressed
Hc = zeros(NB*NF,BNN*FNN);
%%% IdxJc is like IdxJ but compressed
IdxJc = zeros(NB*NF,BNN*FNN);
%%% IdxIc = repmat((1:(NB*NF))',1,BNN*FNN);

cback=0;
for ii=1:NB
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('%4d',ii);
    
    Aii = pcaBASIS(:,1:D,ii);
    
    for kk=1:BNN
        if BaseWeight(ii,kk)>0
            jj = BaseIndex(ii,kk);
            Ajj = pcaBASIS(:,1:D,jj);
            
            H = Ajj'*Aii;
            [U,~,V] = svd(H);
            Pji = V*U'; %%% approximates parallel transport from j to i
            %%% Pij = U*V' approximates the parallel transport from i to j
            %%% but here we really want Pji since we set the unit sphere
            %%% at sample i fixed
            
%             %%% testing kdtree_k_nearest_neighbors_tg
%             idx = zeros(NF, hdm.FNN);
%             dist = zeros(NF, hdm.FNN);
%             for th = 1:NF
%                 [idx(th,:), dist(th,:)] = kdtree_k_nearest_neighbors(unit_sphere_kdtree, (Pji*unit_sphere(:,th))', hdm.FNN);
%             end
            
            %%% fix the unit sphere at sample i and compute how its points
            %%% spread to the unit sphere at sample j
            [Idx, Dist] = kdtree_k_nearest_neighbors_tg(unit_sphere_kdtree, (Pji*unit_sphere)', hdm.FNN);
            HcBlockRowRange = ((ii-1)*NF+1):(ii*NF);
            HcBlockColRange = ((kk-1)*FNN+1):(kk*FNN);
            Hc(HcBlockRowRange,HcBlockColRange) = BaseWeight(ii,kk)*exp(-Dist.^2/hdm.BEpsilon);
            IdxJc(HcBlockRowRange,HcBlockColRange) = (jj-1)*BNN*hdm.FNN+Idx;
        end
    end
end
if hdm.debug==1
    fprintf('\n');
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%% symmetrize Hc in place: check this part again and again!
if hdm.debug==1
    fprintf('(DEBUG:HDM) Step 4'': symmetrize the heat kernel in place\t\t');
end

cback=0;
for ii=1:NB
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('%4d',ii);
    
    for kk=1:BNN
        if BaseWeight(ii,kk)>0
            jj = BaseIndex(ii,kk);
            if (jj > ii)
                ss = find(BaseIndex(jj,:) == ii);
                %%% ss is the column index in the jj-th row of BaseWeight that
                %%% contains index ii
                iikkBlockRowRange = ((ii-1)*NF+1):(ii*NF);
                iikkBlockColRange = ((kk-1)*FNN+1):(kk*FNN);
                jjssBlockRowRange = ((jj-1)*NF+1):(jj*NF);
                jjssBlockColRange = ((ss-1)*FNN+1):(ss*FNN);
                
                IdxI = repmat((1:NF)',hdm.FNN,1);
                iikkIdxJ = IdxJc(iikkBlockRowRange,iikkBlockColRange)-(jj-1)*BNN*hdm.FNN;
                iikkBlockExpand = sparse(flat(IdxI), flat(iikkIdxJ), flat(Hc(iikkBlockRowRange,iikkBlockColRange)));
                
                jjssIdxJ = IdxJc(jjssBlockRowRange,jjssBlockColRange)-(ss-1)*BNN*hdm.FNN;
                jjssBlockExpand = sparse(flat(IdxI), flat(jjssIdxJ), flat(Hc(jjssBlockRowRange,jjssBlockColRange)));
                
                tmp = min(iikkBlockExpand, jjssBlockExpand');
                for rr = 1:NF
                    for tt=1:hdm.FNN
                        Hc(iikkBlockRowRange(rr),iikkBlockColRange(tt)) = tmp(rr,iikkIdxJ(rr,tt));
                        Hc(iikkBlockRowRange(rr),iikkBlockColRange(tt)) = tmp(rr,jjssIdxJ(rr,tt));
                    end
                end
            end
        end
    end
end

if hdm.debug==1
    fprintf('\n');
end




cback=0;
for ii=1:NB
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('%4d',ii);
    
    pii = sum(exp(-5*(BaseWeight(ii,:).^2/hdm.BEpsilon)))^hdm.alpha;
    W2i = 0;
    
    for kk = 2:BNN
        jj = BaseIndex(ii,kk);
        
        %%% use Gaussian kernel exp(-5t^2)
        Kij = exp(-5*(BaseWeight(ii,kk).^2/hdm.BEpsilon));
        pjj = sum(exp(-5*(BaseWeight(jj,:).^2/hdm.BEpsilon)))^hdm.alpha;
        
        Kij2 = Kij./(pii*pjj);
        
        W2i = W2i+Kij2;
        
        Ai = pcaBASIS(:,1:D,ii);
        Aj = pcaBASIS(:,1:D,jj);
        H = Ai'*Aj; %%% H is the parallel transport from fibre j to fibre i
        [U,~,V] = svd(H);
        X1 = V*U';
        
        if strcmpi(hdm.FSampleType,'uniform')
            TransDist = zeros(NF,hdm.FNN);
            TransIdxJ = zeros(NF,hdm.FNN);
            for jjj=1:NF
                [TransIdxJ(jjj,:), TransDist(jjj,:)] = kdtree_k_nearest_neighbors(unit_sphere_kdtree,(X1*unit_sphere(:,jjj))',hdm.FNN);
%                 CheckIdx = jjj;
%                 clf;hold on;axis equal;
%                 scatter(unit_sphere(1,:),unit_sphere(2,:),5,'r','filled'); %% this is fibre ii
%                 scatter(X1(1,:)*unit_sphere,X1(2,:)*unit_sphere,5,'k','filled'); %% this is fibre jj
%                 scatter(X1(1,:)*unit_sphere(:,CheckIdx),X1(2,:)*unit_sphere(:,CheckIdx),20,'k','filled');
%                 scatter(unit_sphere(1,TransIdxJ(CheckIdx,:)),unit_sphere(2,TransIdxJ(CheckIdx,:)),20,'b'); %% this is fibre jj
%                 pause();
            end
            TransIdxJ = TransIdxJ(:,end:-1:1);
            TransDist = TransDist(:,end:-1:1);
            %%% view nearest neighbors
%             figure;
%             scatter3(hdm.data(1,:),hdm.data(2,:),hdm.data(3,:),1,'b','filled');
%             axis equal
%             hold on
%             CheckIdx = ii;
%             scatter3(hdm.data(1,CheckIdx),hdm.data(2,CheckIdx),hdm.data(3,CheckIdx),5,'r','filled');
%             tangent_basis = pcaBASIS(:,:,CheckIdx);
%             basis = tangent_basis+repmat(hdm.data(:,CheckIdx),1,2);
%             scatter3(basis(1,:),basis(2,:),basis(3,:),10,'g','filled');
%             line([hdm.data(1,CheckIdx),basis(1,1)],[hdm.data(2,CheckIdx),basis(2,1)],[hdm.data(3,CheckIdx),basis(3,1)],'color','g');
%             line([hdm.data(1,CheckIdx),basis(1,2)],[hdm.data(2,CheckIdx),basis(2,2)],[hdm.data(3,CheckIdx),basis(3,2)],'color','g');
%             unit_sphere_ii = tangent_basis*unit_sphere+repmat(hdm.data(:,CheckIdx),1,FibreSize);
%             scatter3(unit_sphere_ii(1,:),unit_sphere_ii(2,:),unit_sphere_ii(3,:),10,'r','filled');
%             
%             CheckIdx = jj;
%             scatter3(hdm.data(1,CheckIdx),hdm.data(2,CheckIdx),hdm.data(3,CheckIdx),5,'k','filled');
%             tangent_basis = pcaBASIS(:,:,CheckIdx);
%             basis = tangent_basis+repmat(hdm.data(:,CheckIdx),1,2);
%             scatter3(basis(1,:),basis(2,:),basis(3,:),10,'y','filled');
%             line([hdm.data(1,CheckIdx),basis(1,1)],[hdm.data(2,CheckIdx),basis(2,1)],[hdm.data(3,CheckIdx),basis(3,1)],'color','g');
%             line([hdm.data(1,CheckIdx),basis(1,2)],[hdm.data(2,CheckIdx),basis(2,2)],[hdm.data(3,CheckIdx),basis(3,2)],'color','g');
%             unit_sphere_jj = tangent_basis*unit_sphere+repmat(hdm.data(:,CheckIdx),1,FibreSize);
%             scatter3(unit_sphere_jj(1,:),unit_sphere_jj(2,:),unit_sphere_jj(3,:),10,'k','filled');

            TransIdxI = repmat((1:NF)',1,hdm.FNN);
            TransVals = exp(-5*TransDist.^2/hdm.FEpsilon);
            TransBlock = sparse(TransIdxI(:),TransIdxJ(:),TransVals(:),NF,NF);
            TransBlock = max(TransBlock,TransBlock'); %%% need a better way to do symmetrization
            TransBlock = sparse(diag(1./sum(TransBlock,2).^hdm.alpha))*TransBlock*sparse(diag(1./sum(TransBlock).^hdm.alpha));
        else
            %%% estimate the unit sphere from data
        end
        
        Hc((kk-1)*NF+1:kk*NF, (ii-1)*NF+1:ii*NF) = TransBlock*Kij2; 
    end
    Dc((ii-1)*NF+1:ii*NF) = W2i;
end

fprintf('\n');

%%% the following code are used to get the sparse matrix for either the
%%% connection Laplacian or its heat kernel
I = repmat(1:NF*NB, NF*BNN, 1);
J = zeros(NF*BNN,NF*NB);
for ii=1:NB
    H = zeros(NF*BNN,1);
    for jj=1:BNN
        kk = BaseIndex(ii,jj);
        H((jj-1)*NF+1:jj*NF) = (kk-1)*NF+1:kk*NF;
    end
    for jj=1:NF
        J(:,(ii-1)*NF+jj) = H;
    end
end

%%% get A for the heat kernel
sparseA = sparse(I(:),J(:),Hc(:),NF*NB,NF*NB,NF*BNN*NF*NB);
clear Ac Cc I J H REFLECTION

%%% get D^{-1/2}
I = 1:NF*NB;
sparseD = sparse(I(:),I(:),1./Dc,NF*NB,NF*NB,NF*NB);

%%% get \tilde{S}=D^{-1/2}AD^{-1/2}, which is similar to S (for heat kernel)
sparseS = sqrt(sparseD)*sparseA*sqrt(sparseD);

%%% symmetrize sparseS
%%% not quite necessary due to the constructoin, but helps with numerics
sparseS = (sparseS+sparseS')/2;
rslt.sparseS = sparseS;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Step 5: Spectral Decomposition
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if hdm.debug==1
    fprintf('(DEBUG:HDM) step 5: find eigenvalues of the connection Laplacian operator\n');
end

%%% heavy-lifting eigen-decomposition of the heat kernel
[US,lambdaS]	   = eigs(sparseS, hdm.embedmaxdim, 'lm', eigopt);
lambdaS		       = diag(lambdaS);
[lambdaS, sortidx] = sort(lambdaS, 'descend');
US = US(:, sortidx);	%% these are the eigenvectors for \tilde{S}
US = sqrt(sparseD)*US;	%% these are the eigenvectors for S

rslt.US = US;
rslt.lambdaS = lambdaS;
clear UC lambdaC sparseC

if hdm.debug==1
    subplot(2,1,1); bar(lambdaS);
    set(gca,'fontsize',hdm.fontsize);
    title(['(DEBUG) the first ' num2str(hdm.embedmaxdim) ' eigenvalues\newline(note the scale)']);
    axis tight;
    
    subplot(2,1,2); bar((lambdaS./lambdaS(1)).^(2*hdm.T));
    set(gca,'fontsize',hdm.fontsize);
    title('(DEBUG) the diffusion behavior of the eigenvalues\newline(note the scale)');
    axis tight;
end

if hdm.debug==1
    fprintf(['(DEBUG:HDM) The diffusion time T=',num2str(hdm.T),', and the threshold is ',num2str(hdm.delta),'\n']);
end

dimidx = find((lambdaS./lambdaS(1)).^(2*hdm.T) > hdm.delta);
dimno = length(dimidx);
rslt.embeddim = dimno*(dimno+1)./2;

fprintf(['\t\tHypoelliptic diffusion map will embed the dataset into ',num2str(rslt.embeddim),'-dim Euclidean space\n\n']);

end
