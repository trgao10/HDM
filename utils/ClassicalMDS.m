function [sv, coord] = ClassicalMDS(D, dim)

% check if D is symmetric
% max(max(abs(D-D')))

N =size(D,1);

DC = -.5*(D.^2 - sum(D.^2)'*ones(1,N)/N -...
   ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2));

%E = ones(N, N);
%DC = -.5*(D.^2-E*D.^2/N-D.^2*E/N+E*D.^2*E/N^2);

% check if DC is symmetric
disp(['symmetric error = ' num2str(max(max(abs(DC-DC'))))]);
DC = (DC+DC')/2;

% [coord, sv] = eigs(DC, dim);
[coord, sv] = eigs(DC,dim-1);