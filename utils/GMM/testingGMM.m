

%%load image
vector =@(X)reshape(X,size(X,1)*size(X,2),3);

%%load image
name = 'parrot';
X0 = double(imread(['../../data/' name '-1.jpg']))/256;

N0 = size(X0,1);

%% Construct data set X by subsampling X0
%we subsample to compute the optimal mapping
K = 8; %N:number of samples small image
d=3;

X0 = reshape(X0,N0*N0,d);

[muX,SigmaX,PX]=GMM(X0,K);
figure;
subplot(1,2,1);
for k=1:K
    plot_ellipse(SigmaX(1:d,1:d,k)',muX(k,1:d),1.5,'none',rand(3,1));hold on;axis tight;axis equal
end
title('Gaussians')
subplot(1,2,2);plot3(X0(:,1),X0(:,2),X0(:,3),'.b');hold on;
axis tight;axis equal;title('Data')


