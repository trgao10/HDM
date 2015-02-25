%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Generate N points uniformly on a sphere in Dim dimensions %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load path
clear vars;
path(pathdef);
path(path, genpath('./utils/'));

%% parameters
Dim = 3; %% Dim = 2 or 3
N = 2000;
noise_level = 0;
rng('shuffle');

%%% SPHERE
x = randn(Dim, N);
norm_x = sqrt(sum(x.^2));
x = x * diag(1./norm_x);

%%% visualize
if Dim==3
    plot3(x(1,:), x(2,:), x(3,:), '.');
else
    plot3(x(1,:), x(2,:), zeros(1,N), '.');
end
axis equal;
title('Original Data Set');
set(gca, 'FontSize', 8);
cameratoolbar;
cameratoolbar('SetCoordSys','none');

%% vector diffusion map
% vdm.data	    = x;
% vdm.epsilon	    = 1;
% vdm.NN          = 100;
% vdm.T           = 2;
% vdm.delta       = 0.9;
% vdm.fontsize    = 8;
% vdm.symmetrize  = 0;
% vdm.Nabla2      = 1;
% vdm.debug       = 1;
% vdm.embedmaxdim = 100;
% 
% vdm.lambdaNO = vdm.embedmaxdim; %% compatibility issue with H.-T. Wu's original implementation
% 
% [rslt] = VecDiffMap(vdm);

%% hypoelliptic diffusion map
hdm.data	 = x;
hdm.NF       = 100;
hdm.BEpsilon = 0.2;
hdm.FEpsilon = 1e-5;
hdm.BNN	     = 100;
hdm.FNN	     = 100;
hdm.T	     = 1;
hdm.delta	 = 0.9;
hdm.embedmaxdim = 100;
hdm.symmetrize = 1;
hdm.FSampleType = 'uniform';
hdm.debug = 1;
hdm.fontsize = 8;

[rslt] = HypoellipticDiffusionMap(hdm);

% hdm.NF = 1;
% [rslt2] = HypoellipticDiffusionMap(hdm);
