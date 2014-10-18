function [h_in, h_out] = plot_ellipse(T,m,r, color_bound, color_inside)

% plot_ellipse - draw ellipse
%
%   hb = plot_ellipse(T,m,r, color_bound, color_inside);
%
%   T is the tensor of the ellipse
%   m is the center
%   r is the radius
%
%   Copyright (c) 2010 Gabriel Peyre

if nargin<4
    color_bound = [0 0 1];
end
if nargin<5
    color_inside = [1 0 0];
end
if nargin<3
    r = 1;
end


n = 48;
%[V,S] = eig(T);
[V,S,U] = svd(T);


[u,v]=meshgrid(linspace(0,2*pi,n),linspace(0,2*pi,n));
flat=@(x)x(:)';
x=cat(1,flat(cos(u).*cos(v)),flat(cos(u).*sin(v)),flat(sin(u)));
y = r * V*sqrt(S)*x + repmat(m(:), [1 n*n]);
y=real(y);

S = sqrt(S);
quiver3(m(1),m(2),m(3),r*S(1,1)*V(1,1),r*S(1,1)*V(2,1),r*S(1,1)*V(3,1),'Color',color_inside);
quiver3(m(1),m(2),m(3),r*S(2,2)*V(1,2),r*S(2,2)*V(2,2),r*S(2,2)*V(3,2),'Color',color_inside);
quiver3(m(1),m(2),m(3),r*S(3,3)*V(1,3),r*S(3,3)*V(2,3),r*S(3,3)*V(3,3),'Color',color_inside);

%h_out = plot3(y(1,:),y(2,:),y(3,:), 'color', color_bound);
C=cat(3,color_inside(1).*ones(n,n),color_inside(2).*ones(n,n),color_inside(3).*ones(n,n));
surf(reshape(y(1,:),[n n]),reshape(y(2,:),[n n]),reshape(y(3,:),[n n]),C,'EdgeColor',color_bound);
alpha(0.1)

