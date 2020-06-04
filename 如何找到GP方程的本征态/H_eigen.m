function [ H1, E1, U1 ] = H_eigen ( V, N, X, dx, nx )

% Diagonal the potential energy
% V=lambdai/2*(X-x0).^(2*n);% power-law potential ( lambdai, x0, n, N, X, dx, nx )
% V = Ej*cos(2*pi * X);% optical lattice + tilt ( U0, kL, m, a, N, X, dx, nx, hb )
% V = U0/6 * q^2 * (q0-q);

diagV=spdiags(V',0:0,nx,nx);

% Diagonal the kinetic energy
e=ones(nx,1);
H00=spdiags([e -2*e e],-1:1,nx,nx);
H00=-H00/dx^2;

% find sp eigenstates
H=H00/2+diagV;
OPTS1.disp=0;
Nstates = N;
[u1,GF]=eigs(H,Nstates,'SM',OPTS1);
% [u1,GF]=eig(H);
[G,Ind]=sort(diag(GF));

u=zeros(length(X),Nstates);
for ii=1:Nstates
    u(:,ii)=u1(:,Ind(ii))/sqrt(sum(abs(u1(:,Ind(ii))).^2)*dx);
end

H1=H;
E1=G(1:N);
U1=u(:,1:N);

clear u V diagV H00 H GF G Ind u1;

end

