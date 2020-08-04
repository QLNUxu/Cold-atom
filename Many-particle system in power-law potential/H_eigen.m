function [ H1, E1, U1 ] = H_eigen( lambdai, x0, n, N, X, dx, nx )

% global X dx nx 

% Diagonal the potential energy
V = lambdai/2*(X-x0).^(2*n);% power-law potential

diagV = spdiags(V',0:0,nx,nx);

% Diagonal the kinetic energy
e=ones(nx,1);
H00=spdiags([e -2*e e],-1:1,nx,nx);
H00=-H00/dx^2;

%No of single particle states used to diagonalise V 
Nparticles=N;

% find sp eigenstates
OPTS1.disp=0;
Nstates=Nparticles;

H=H00./2+diagV;
[u1,GF]=eigs(H,Nstates,'SM',OPTS1);
[G,Ind]=sort(diag(GF));

u=zeros(length(X),Nparticles);
for ii=1:Nparticles
    u(:,ii)=u1(:,Ind(ii))/sqrt(sum(abs(u1(:,Ind(ii))).^2)*dx);
end

H1=H;
clear u1;
E1=G(1:N);
U1=u(:,1:N);
clear V diagV H00 H GF G Ind u

end

