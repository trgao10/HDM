function [rslt] = vectordiffusionmap(vdm)
%
%  Vector diffusion map v0.1
%
% INPUT:
%   vdm.data	: pxn matrix which represents n data points in R^p
%   vdm.epsilon	: \epsilon in the vDM
%   vdm.NN	: number of nearest neighbors
%   vdm.T	: diffusion time
%   vdm.delta	: parameter for the truncated VDM
%   vdm.origidx	: the index of the "origin point" you want to measure the affinity
%   vdm.symmetrize : symmetrize the graph
%   vdm.compact : compact support kernel (TODO)
%
% (OPTION)
%   vdm.debug   : debug mode
%   vdm.parallel: parallel computation (TODO)
%   vdm.draw    : draw vector fields
%   vdm.drawNO  : draw the first drawNO vector fields
%
% OUTPUT:
%   rslt.embeddim : the truncated VDM
%   rslt.embed    : the embedded data
%   rslt.vdd	 : the vector diffusion distance of each point to the origidx point
%
% DEPENDENCE:
%   tstoolbox, localpca.m, (TESTING: KN_rankest.m)
%
% by Hau-tieng Wu 2011-06-28
%

if vdm.debug; fprintf('\n(DEBUG:vDM)\t\tStart to work on vector diffusion map\n'); end
clc;
eigopt.isreal = 1;
eigopt.issym  = 1;
eigopt.maxit  = 3000;
eigopt.disp   = 0;


[pp,nn] = size(vdm.data);


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
if vdm.debug; fprintf('(DEBUG:vDM) Step 1: NN search and Data preparation. (For comuptational load issue)\n'); end;

X = (vdm.data).'; atria = nn_prepare(X);
[index,distance] = nn_search(X, atria, [1:nn].', vdm.NN, -1, 0.0);

if vdm.debug
    fprintf(['(DEBUG:vDM) NN=',num2str(vdm.NN),'\n']);
    fprintf(['(DEBUG:vDM) minimal farest distance=',num2str(min(distance(:,end))),'\n']);
    fprintf(['(DEBUG:vDM) maximal farest distance=',num2str(max(distance(:,end))),'\n']);
    fprintf(['(DEBUG:vDM) median farest distance=',num2str(median(distance(:,end))),'\n']);
    fprintf(['(DEBUG:vDM) 1.5*sqrt(min farest distance)=',num2str(1.5*sqrt(min(distance(:,end)))),'.\n']);
    fprintf(['(DEBUG:vDM) 1.5*sqrt(max farest distance)=',num2str(1.5*sqrt(max(distance(:,end)))),'.\n']);
    fprintf(['(DEBUG:vDM) 1.5*sqrt(median farest distance)=',num2str(1.5*sqrt(median(distance(:,end)))),'.\n']);
end

%% patchno is set to convert the NN information to the \sqrt{h} information
patchno = vdm.NN*ones(1,nn);

if vdm.debug==1; fprintf('(DEBUG:vDM) the neighbors with kernel value less than exp(-5*1.5*1.5)=1.3e-5 are trimmed.\n'); end

for ii = 1:nn
    patchno(ii) = length(find(distance(ii,:) <= 1.5*sqrt(vdm.epsilon)));
    distance( ii, find(distance(ii,:)>1.5*sqrt(vdm.epsilon)) ) = inf;
end

distance = distance(:, 1:max(patchno));
index = index(:, 1:max(patchno));
rslt.patchno = patchno; rslt.distance = distance; rslt.index = index;

if vdm.debug==1
    if quantile(patchno,0.9) == vdm.NN
        %% it means the NN is not big enough so that the decay of the kernel will not be enough and the error will be large
        fprintf('(DEBUG:vDM:WARNING) the NN should be chosen larger\n');
    end
    
    fprintf(['(DEBUG:vDM) the number of points with distance less then\n']);
    fprintf(['           1.5*sqrt(epsilon)=',num2str(1.5*sqrt(vdm.epsilon)),' is (min,max,median) = (',num2str(min(patchno)),',',num2str(max(patchno)),',',num2str(median(patchno)),')\n']);
    fprintf(['(DEBUG:vDM) set NN to be ',num2str(max(patchno)),'\n']);
    
end

NN = max(patchno);


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
if vdm.debug; fprintf('(DEBUG:vDM) Step 2: find a basis for each tangent plane by PCA\n'); end


lpca.data      = vdm.data;
%% (CAUTION) theoretically, it should be vdm.epsilon^((D+4)/(D+1)); but we don't know D yet, so I just use 0.8 a random guess. (TODO)
lpca.epsilonpca= vdm.epsilon*0.8;
lpca.NN        = NN;
lpca.index     = rslt.index;
lpca.distance  = rslt.distance;
lpca.patchno   = patchno;
%% lpca.KN is the threshold number gamma chosen by the user in the local PCA algorithm. By default we choose 0.9
lpca.KN        = 0.9;
lpca.debug     = 1;

[pcaBASIS, D]  = LocalPCA(lpca);
clear lpca; clear X;


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
if vdm.debug; fprintf('(DEBUG:vDM) Step 3: deal with the reflection effect\t\t'); end
REFLECTION=diag(ones(1,nn));

cback=0;
for ii=1:nn
    for cc=1:cback; fprintf('\b'); end; cback = fprintf('%4d',ii);
    
    tmpSLICER = REFLECTION(ii,:); 	%% prepare for parallelization
    
    for kk = 2:NN
        jj = index(ii,kk);	%% ii and jj are indices of points
        
        Ai = pcaBASIS(:,1:D,ii);
        Aj = pcaBASIS(:,1:D,jj);
        H  = Ai'*Aj;
        [U, lambda, V] = svd(H);
        X1 = V*U';
        
        if det(X1)<0; tmpSLICER(jj) = -1; else; tmpSLICER(jj) = 1; end
    end
    
    REFLECTION(ii,:)=tmpSLICER;
end
fprintf('\n'); clear tmpSLICER; clear U; clear V; clear lambda; clear X1;

[UR,lambdaR] = eigs(REFLECTION,2,'lm',eigopt);

%% make all frames coincide with SO(d)
for ii=1:nn; if UR(ii,1)<0
        frame = pcaBASIS(:,:,ii); frame(:,1) = -frame(:,1); pcaBASIS(:,:,ii) = frame;
    end; end;

clear UR; clear lambdaR; clear REFLECTION; clear frame;
rslt.pcaBASIS  = pcaBASIS;


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% prepare for the graph
if vdm.symmetrize==1
    if vdm.debug; fprintf('(DEBUG:vDM) Step 4'': symmetrize the graph (very time consuming) \t\t'); end
    count = 0; cback = 0;
    for ii=1:nn;
        for cc = 1:cback; fprintf('\b'); end; cback = fprintf('%4d',ii);
        for kk = 2:NN
            jj = index(ii,kk);
            if ~ismember(ii,index(jj,:)); distance(ii,kk)=inf; count=count+1; end
        end;
    end;
    
    if vdm.debug; fprintf(['(DEBUG:vDM) removed entries=',num2str(100*count./(size(distance,1)*size(distance,2))),'%%\n\n']); end
end
clear count;

%% build up the matrix S=D^{-1}A
Ac = zeros(D*NN, D*nn);
Dc = zeros(D*nn,1);

%% vdm.epsilon=\epsilon_{PCA}. Here is the \epsilon
epsilon = vdm.epsilon;

if vdm.debug; fprintf('(DEBUG:vDM) Step 4: connection Laplacian operator\t\t'); end

%% (TODO) \alpha power
%% work out "entries of A and D". For the sparse matrix purpose
cback=0;
for ii=1:nn
    for cc=1:cback; fprintf('\b'); end; cback = fprintf('%4d',ii);
    
    pii = sum(exp(-5*(distance(ii,:).^2/epsilon)));
    W2i = 0;
    
    for kk = 1:NN
        
        jj = index(ii,kk);
        
        %% here we use kernel exp(-5t^2)
        %% TODO: compact support kernel
        Kij = exp(-5*(distance(ii,kk).^2/epsilon));
        pjj = sum(exp(-5*(distance(jj,:).^2/epsilon)));
        
        Kij2 = Kij./(pii*pjj);
        
        W2i = W2i+Kij2;
        
        Ai = pcaBASIS(:,1:D,ii);
        Aj = pcaBASIS(:,1:D,jj);
        H = Ai'*Aj;
        [U,lambda,V] = svd(H);
        X1 = V*U';
        
        Ac((kk-1)*D+1:kk*D,(ii-1)*D+1:ii*D) = X1*Kij2;
        
    end
    
    Dc((ii-1)*D+1:ii*D) = W2i;
    
end

fprintf('\n');



%% preparation for the construction of the Connection Laplace. The above is the preparation for its heat kernel
if isfield(vdm,'Nabla2'); if vdm.Nabla2
        
        if vdm.debug; fprintf('(DEBUG:vDM) Get Connection Laplacian...\n'); end
        Cc = Ac;
        
        for ii=1:nn
            Cc(1:D,(ii-1)*D+1:ii*D) = (1/vdm.epsilon) * Ac(1:D,(ii-1)*D+1:ii*D)-diag(ones(D,1));
        end
        
    end; end

%% the following code are used to get the sparse matrix for either the connection Laplacian or its heat kernel
I = zeros(D*NN,D*nn);
for ii=1:D*NN; I(ii,:)=1:D*nn; end

J = zeros(D*NN,D*nn);
for ii = 1:nn
    H = zeros(D*NN,1);
    
    for jj = 1:NN; kk = index(ii,jj); H((jj-1)*D+1:jj*D) = (kk-1)*D+1:kk*D; end;
    for jj = 1:D; J(:,(ii-1)*D+jj)=H; end
end


%% get A for the heat kernel
sparseA = sparse(I(:),J(:),Ac(:),D*nn,D*nn,D*NN*D*nn);

%% get A for connection Laplacian
if isfield(vdm,'Nabla2'); if vdm.Nabla2
        sparseC = sparse(I(:),J(:),Cc(:),D*nn,D*nn,D*NN*D*nn);
    end; end

clear Ac; clear Cc; clear I; clear J; clear H; clear REFLECTION;

I=[1:1:D*nn];

%% get D^{-1/2}
sparseD = sparse(I(:),I(:),1./Dc,D*nn,D*nn,D*nn);

%% get \tilde{S}=D^{-1/2}AD^{-1/2}, which is similar to S (for heat kernel)
sparseS = sqrt(sparseD)*sparseA*sqrt(sparseD);

%% since numerically sparseS is almost symmetric by 1e-18, this step won't generate too much error but should help the numerical result
sparseS = (sparseS+sparseS')/2;
rslt.sparseS = sparseS;

%% get S=D^{-1}A, which the true heat kernel
if isfield(vdm,'Nabla2'); if vdm.Nabla2
        sparseC = sparseD*sparseC;
    end; end;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
if vdm.debug; fprintf('(DEBUG:vDM) step 5: find eigenvalues of the connection Laplacian operator\n'); end

[US,lambdaS]	= eigs(sparseS, vdm.lambdaNO, 'lm');
lambdaS		= diag(lambdaS);
[lambdaS,sortidx]= sort(lambdaS, 'descend');
US		= US(:,sortidx);	%% these are the eigenvectors for \tilde{S}
US		= sqrt(sparseD)*US;	%% these are the eigenvectors for S

rslt.US = US;
rslt.lambdaS = lambdaS;

%% for connection laplace
if isfield(vdm,'Nabla2')
    if vdm.Nabla2
        [UC,lambdaC] = eigs(sparseC, vdm.lambdaNO, 'sm');
        lambdaC      = -4*real(diag(lambdaC))-(D-1);
%         lambdaC      = -4*real(diag(lambdaC))-(D-1);
        rslt.lambdaC = lambdaC;
        rslt.UC	 = UC;
    end
end;
clear UC lambdaC sparseC;

if vdm.debug
    subplot(2,1,1); bar(lambdaS);
    set(gca,'fontsize',vdm.fontsize);
    title('(DEBUG) the first 200 eigenvalues\newline(note the scale)');
    axis tight;
    
    subplot(2,1,2); bar((lambdaS./lambdaS(1)).^(2*vdm.T));
    set(gca,'fontsize',vdm.fontsize);
    title('(DEBUG) the diffusion behavior of the eigenvalues\newline(note the scale)');
    axis tight;
    
    figure;
    subplot(2,1,1); bar(rslt.lambdaC);
    set(gca,'fontsize',vdm.fontsize);
    title(['(DEBUG) the first ' num2str(vdm.embedmaxdim) ' eigenvalues\newline(note the scale)']);
    axis tight;
    subplot(2,1,2); bar((rslt.lambdaC./rslt.lambdaC(1)-1).^(2*vdm.T));
    set(gca,'fontsize',vdm.fontsize);
    title(['(DEBUG) the diffusion behavior of the eigenvalues\newline(note the scale) T=' num2str(vdm.T)]);
    axis tight;
end

%
%% TODO: setup some parameters for special manifolds, for example, S^n
%
% lambdaS(1:6)=lambdaS(1); %% if sphere, set the first 6 eigenvalues the same artificially
%

if vdm.debug
    fprintf(['(DEBUG:vDM) The diffusion time T=',num2str(vdm.T),', and the threshold is ',num2str(vdm.delta),'\n']);
end

dimidx = find((lambdaS./lambdaS(1)).^(2*vdm.T) > vdm.delta);
dimno = length(dimidx);
rslt.embeddim = dimno*(dimno+1)./2;

fprintf(['\t\tVector diffusion map will embed the dataset into ',num2str(rslt.embeddim),'-dim Euclidean space\n\n']);

if vdm.debug && rslt.embeddim>500
    fprintf('(DEBUG:vDM:WARNING) the embedding dimension might be too big, really continue?\n'); pause;
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%%
%%% The following code is written to save space. Since the VDM
%%% might embed the data to much higher dimensional space,
%%% be careful if you want to save all the embedded data
%%% To save space, we only save the vector diffusion distance
%%%
if isfield(vdm,'vdd'); if vdm.vdd
        
        if vdm.debug; fprintf('(DEBUG:vDM) step 6: VDM and VDD\t\t'); end
        
        x0pt = zeros(rslt.embeddim, 1);
        
        ss=1;
        for ii = 1:dimno
            for jj = ii:dimno
                if ii ~= jj
                    x0pt(ss) = sqrt(2)*((lambdaS(ii)*lambdaS(jj))^(vdm.T))*...
                        US((vdm.origidx-1)*D+1 : vdm.origidx*D, ii)'*...
                        US((vdm.origidx-1)*D+1 : vdm.origidx*D, jj);
                else
                    x0pt(ss) = ((lambdaS(ii)*lambdaS(jj))^(vdm.T))*...
                        US((vdm.origidx-1)*D+1 : vdm.origidx*D, ii)'*...
                        US((vdm.origidx-1)*D+1 : vdm.origidx*D, jj);
                end
                
                ss=ss+1;
            end
        end
        
        
        rslt.vdd = zeros(1,nn);
        
        xipt = zeros(size(x0pt));
        tic; cback = 0;
        for qq = 1:nn
            for cc=1:cback; fprintf('\b'); end; cback = fprintf('%4d',qq);
            ss = 1;
            for ii = 1:dimno
                for jj = ii:dimno
                    if ii ~= jj
                        xipt(ss) = sqrt(2)*((lambdaS(ii)*lambdaS(jj))^(vdm.T))*...
                            US((qq-1)*D+1:qq*D, ii)'*...
                            US((qq-1)*D+1:qq*D, jj);
                    else
                        xipt(ss) = ((lambdaS(ii)*lambdaS(jj))^(vdm.T))*...
                            US((qq-1)*D+1:qq*D, ii)'*...
                            US((qq-1)*D+1:qq*D, jj);
                    end
                    ss=ss+1;
                end
            end
            rslt.vdd(qq)=norm(xipt-x0pt);
        end
        fprintf('\n'); toc;
        
        if vdm.debug
            if pp==3
                figure; scatter3(vdm.data(1,:),vdm.data(2,:),vdm.data(3,:),10,rslt.vdd(:),'fill');
                hold on; plot3(vdm.data(1,vdm.origidx),vdm.data(2,vdm.origidx),vdm.data(3,vdm.origidx),'. red','markersize',30);
                view([30,90]);	% 90 for torus, 180 for S^2
                colorbar; set(gca,'fontsize',vdm.fontsize); axis tight; axis off
                title('The vector diffusion distance from the red point');
                
            elseif pp==2
                figure; scatter(vdm.data(1,:),vdm.data(2,:),10,rslt.vdd(:),'fill');
                hold on; plot(vdm.data(1,vdm.origidx),vdm.data(2,vdm.origidx),'. red','markersize',30);
                colorbar; set(gca,'fontsize',vdm.fontsize); axis tight; axis off;
                title('The vector diffusion distance from the red point');
                
            elseif pp==1
                figure; plot(vdm.data, rslt.vdd, '.');
                set(gca,'fontsize',vdm.fontsize); hold on;
                plot(vdm.data(vdm.origidx),rslt.vdd(vdm.origidx),'. red','markersize',30);
                title('The vector diffusion distance from the red point');
            else
                fprintf('Sorry, but I don''t know how to visualize the high dimensional manifold\n');
            end
        end
        
    end
    
end
