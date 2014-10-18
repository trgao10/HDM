function [mu,Sigma,z] = GMM(data,K)
N = size(data,1);
d = size(data,2);

Sigma = zeros(d,d,K);

tau= zeros(N,K);

Xinit = rand(K,d);
[mu,~] = knnsubsampling(data,Xinit);
z = gmm_knnsearch(data,mu,1);
for i=1:K
    index = z==i;
    num =sum(index);
    
    tau(index,:) =1;% zeros everywhere but in the 'i' position
    
    dd=data(index,:)-repmat(mu(i,:),[num 1]);
    Sigma(:,:,i) = dd'*dd/num;
end 