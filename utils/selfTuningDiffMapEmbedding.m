function [U,lambda] = selfTuningDiffMapEmbedding(D, K)
%SELFTUNINGDIFFMAPEMBEDDING Summary of this function goes here
%   Detailed explanation goes here

D = D - diag(diag(D));
Weights = sort(D,2);
Weights = Weights(:,K+1);
H = exp(-(D.^2)./sqrt(Weights*Weights'))-eye(size(D,1));
sqrtInvD = sparse(1:size(H,1), 1:size(H,1), 1./sqrt(sum(H)));

eigopt.isreal = 1;
eigopt.issym = 1;
eigopt.maxit = 5000;
eigopt.disp = 0;
[U, lambda] = eigs(sqrtInvD*H*sqrtInvD, size(H,1), 'LM', eigopt);
lambda = diag(lambda);

end

