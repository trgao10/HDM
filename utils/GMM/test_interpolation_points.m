%% Auxiliar functions
plotpoint=@(P,color)plot3(P(1,:),P(2,:),P(3,:),'o','Color',color,'MarkerSize',5);


%% Initial Data: two Gaussians and a point (in E coordinates)
W = [-1 1 0; 1 1 0 ;0 0 2/sqrt(2)]*sqrt(2)/2;
D2 = [0.1; 0.1 ; 0.1];
Sigma(:,:,1) =W*diag(D2)*W';
m(:,1) = [0 0 0];

V = eye(3,3);
D = [0.1; 0.1 ; 1];
Sigma(:,:,2) = V*diag(D)*V';
m(:,2) = [10 10 10];

N=200;
P = (Sigma(:,:,1)^0.5)*randn(3,N)+repmat(m(:,1),[1 N]);


%% Starting: putting the point in Gaussian 1 coordinates

figure;
plot_ellipse(Sigma(:,:,1),m(:,1)',1,'none',[1 0 0]);hold on;
plot_ellipse(Sigma(:,:,2),m(:,2)',1,'none',[0 0 1]);hold on;

plotpoint(P,[1 0 0]);

options.method='ot'; %type of Gaussian interpolation: 'ot','linear','rao'
options.niter=100;
for t=0:0.25:1
    
   [Si,~,~] = compute_gaussian_barycenter(Sigma,[1-t t],options);
   mi = (1-t)*m(:,1)+t*m(:,2);
      
   Pt=pinv(Sigma(:,:,1))*(Sigma(:,:,1)*Si)^0.5*(P-repmat(m(:,1),[1 N]))+repmat(mi,[1 N]);
%"The Fr�chet distance between multivariate normal distributions". D.C Dowson, B.V Landau
% equation 8
  
   c=rand(3,1);
   plot_ellipse(Si,mi,1,'none',c);hold on;
   plotpoint(Pt,c);
   axis equal 
end 
title(options.method);view([0 90]);axis equal 
