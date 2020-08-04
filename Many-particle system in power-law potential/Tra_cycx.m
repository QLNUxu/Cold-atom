% close all;%Remove figure
clear;%Clear all define
% clc;%Clear command window
%% Preparation
tic;
% global X dx nx
m=1;
hb=1;
% boundary
a = -10 ;
b = 20 ;

L=b-a; %space lengh
nx = 1 * 10^3;
dx=L/nx;
X=a+L*(1:nx)/nx;%coordinates separation
P=(2*pi*hb/L)*[0:nx/2-1,-nx/2:-1];%%momentum separation 

Nparticles = 10 ;% How many states
q = 2 ;% Power for X

% transport
xii = 0 ;
xif = 5 ;
lambda = 1;
% a = 1/(3^(1/6));% q=2, n=0
a = (19/543)^(1/6);% q=2, n=10
% a = (101/15303)^(1/6);% q=2, n=50
%% Total cycling

%Time cycle
Ti = 1 ;
TT = 10 ;
dT = 0.02 ;
nT = (TT-Ti)/dT;
CylTim = cumsum(Ti:dT:TT,1); % build x-axis for T

% Fermion = zeros(nT+1,1); % fidelity for TG gas
% T_d1 = zeros(nT+1,1);% trace distance for TG gas
% T_d2 = zeros(nT+1,1);% trace distance for Fermion
% Coh_TG = zeros(nT+1,1);% coherence for TG gas
% Coh_Fer = zeros(nT+1,1);% coherence for Fermion



% Particle cycle
% nn = 10 ;% How many states
% CylN = cumsum(1:nn,1) ;% build x-axis for N
% Fermion = zeros(nn,1); % fidelity for TG gas

% parfor Nparticles = 1:nn  

% Diagonal H
[ H_i, E_i, U_i ] = H_eigen( lambda, xii, q, Nparticles, X, dx, nx );
[ H_T, E_T, U_T ] = H_eigen( lambda, xif, q, Nparticles, X, dx, nx );

% dd=5;
% XX=X(1);
% for jj=dd:dd:nx
%     XX(jj/dd+1) = X(jj);
% end

parfor ii=0:nT

% Define time
nt = 10^3 ; % number of time steps
% T = 8 ;
T = Ti + ii * dT;
dt = T/nt ; % step size
dtm = ones(1,nt); 

% define the initial state
psi=U_i;
% psi_Cos=U_i;

%% Cycling

for j=0:nt
% q=1
%     x0 = (120*(j*(T/nt))^3*xif-180*(j*(T/nt))^2*T*xif+60*(j*(T/nt))*T^2*xif+6*(j*(T/nt))^5*lambda*xif-15*(j*(T/nt))^4*T*lambda*xif+10*(j*(T/nt))^3*T^2*lambda*xif-120*(j*(T/nt))^3*xii+180*(j*(T/nt))^2*T*xii-60*(j*(T/nt))*T^2*xii-6*(j*(T/nt))^5*lambda*xii+15*(j*(T/nt))^4*T*lambda*xii-10*(j*(T/nt))^3*T^2*lambda*xii+T^5*lambda*xii)/(T^5*lambda);
% q=2, n=0
%     x0 = (6*(j*(T/nt))^5*(xif-xii))/T^5-(15*(j*(T/nt))^4*(xif-xii))/T^4+(10*(j*(T/nt))^3*(xif-xii))/T^3+xii+(a^2*lambda)/(2^(1/3)*((-120*(j*(T/nt))^3*lambda^2*xif+180*(j*(T/nt))^2*T*lambda^2*xif-60*(j*(T/nt))*T^2*lambda^2*xif+sqrt(2)*T^5*sqrt((lambda^4*(a^6*T^10*lambda^2+1800*(j*(T/nt))^2*(2*(j*(T/nt))^2-3*(j*(T/nt))*T+T^2)^2*(xif-xii)^2))/T^10)+120*(j*(T/nt))^3*lambda^2*xii-180*(j*(T/nt))^2*T*lambda^2*xii+60*(j*(T/nt))*T^2*lambda^2*xii)/T^5)^(1/3))-((-120*(j*(T/nt))^3*lambda^2*xif+180*(j*(T/nt))^2*T*lambda^2*xif-60*(j*(T/nt))*T^2*lambda^2*xif+sqrt(2)*T^5*sqrt((lambda^4*(a^6*T^10*lambda^2+1800*(j*(T/nt))^2*(2*(j*(T/nt))^2-3*(j*(T/nt))*T+T^2)^2*(xif-xii)^2))/T^10)+120*(j*(T/nt))^3*lambda^2*xii-180*(j*(T/nt))^2*T*lambda^2*xii+60*(j*(T/nt))*T^2*lambda^2*xii)/T^5)^(1/3)/(2^(2/3)*lambda);
% q=2, n=10
    x0 = (6*(j*(T/nt))^5*(xif-xii))/T^5-(15*(j*(T/nt))^4*(xif-xii))/T^4+(10*(j*(T/nt))^3*(xif-xii))/T^3+xii+(19*a^2*lambda)/(2^(1/3)*((-120*(j*(T/nt))^3*lambda^2*xif+180*(j*(T/nt))^2*T*lambda^2*xif-60*(j*(T/nt))*T^2*lambda^2*xif+sqrt(2)*T^5*sqrt((lambda^4*(6859*a^6*T^10*lambda^2+1800*(j*(T/nt))^2*(2*(j*(T/nt))^2-3*(j*(T/nt))*T+T^2)^2*(xif-xii)^2))/T^10)+120*(j*(T/nt))^3*lambda^2*xii-180*(j*(T/nt))^2*T*lambda^2*xii+60*(j*(T/nt))*T^2*lambda^2*xii)/T^5)^(1/3))-((-120*(j*(T/nt))^3*lambda^2*xif+180*(j*(T/nt))^2*T*lambda^2*xif-60*(j*(T/nt))*T^2*lambda^2*xif+sqrt(2)*T^5*sqrt((lambda^4*(6859*a^6*T^10*lambda^2+1800*(j*(T/nt))^2*(2*(j*(T/nt))^2-3*(j*(T/nt))*T+T^2)^2*(xif-xii)^2))/T^10)+120*(j*(T/nt))^3*lambda^2*xii-180*(j*(T/nt))^2*T*lambda^2*xii+60*(j*(T/nt))*T^2*lambda^2*xii)/T^5)^(1/3)/(2^(2/3)*lambda);
% q=2, n=50
%     x0 = (6*(j*(T/nt))^5*(xif-xii))/T^5-(15*(j*(T/nt))^4*(xif-xii))/T^4+(10*(j*(T/nt))^3*(xif-xii))/T^3+xii+(101*a^2*lambda)/(2^(1/3)*((-120*(j*(T/nt))^3*lambda^2*xif+180*(j*(T/nt))^2*T*lambda^2*xif-60*(j*(T/nt))*T^2*lambda^2*xif+sqrt(2)*T^5*sqrt((lambda^4*(1030301*a^6*T^10*lambda^2+1800*(j*(T/nt))^2*(2*(j*(T/nt))^2-3*(j*(T/nt))*T+T^2)^2*(xif-xii)^2))/T^10)+120*(j*(T/nt))^3*lambda^2*xii-180*(j*(T/nt))^2*T*lambda^2*xii+60*(j*(T/nt))*T^2*lambda^2*xii)/T^5)^(1/3))-((-120*(j*(T/nt))^3*lambda^2*xif+180*(j*(T/nt))^2*T*lambda^2*xif-60*(j*(T/nt))*T^2*lambda^2*xif+sqrt(2)*T^5*sqrt((lambda^4*(1030301*a^6*T^10*lambda^2+1800*(j*(T/nt))^2*(2*(j*(T/nt))^2-3*(j*(T/nt))*T+T^2)^2*(xif-xii)^2))/T^10)+120*(j*(T/nt))^3*lambda^2*xii-180*(j*(T/nt))^2*T*lambda^2*xii+60*(j*(T/nt))*T^2*lambda^2*xii)/T^5)^(1/3)/(2^(2/3)*lambda);
% reference polynomial 
%     x0 = (6*(j*(T/nt))^5*(xif-xii))/T^5-(15*(j*(T/nt))^4*(xif-xii))/T^4+(10*(j*(T/nt))^3*(xif-xii))/T^3+xii; 
% reference linear
%     x0 = xii-((j*(T/nt))*(-xif+xii))/T;
    
    K=P.^2/(2*m);% kinetic energy
    UK=exp(-1i*K/hb*dt/2);% propagator in momentum space
    V=(X-x0).^(2*q)*lambda/2;% potential energy
%     V_non=X.^(2*nn)*lambda_non/2;
    UV=exp(-1i*V/hb*dt);% propagator in coordinate space
%     UV_non=exp(-1i*V_non/hb*dt);
    
    psi = Dyn_STA( psi, UV, UK, Nparticles, nx );
%     psi_Non = Dyn_STA( psi_Non, UV_non, UK, n);

end

%% Fidelity
% Fermion(Nparticles,1) = abs(det(U_T'*psi*dx))^2;% particle cycle
Fermion(ii+1,1) = abs(det(U_T'*psi*dx))^2;% time cycle
% Fermion(tt,1) = abs(sum(conj(U_T(:,n)).*psi)*dx)^2;
%% reduced density matrix

% phi=psi(1,:);
% UT=U_T(1,:);
% for jj=dd:dd:nx
%     phi(jj/dd+1,:) = psi(jj,:);
%     UT(jj/dd+1,:) = U_T(jj,:);
% end
%% TG gas
% dkx=(2*pi*hb)./L;
% dkx=1;
% [rspdm,nk,Angps]=dmnk(phi',dx*dd,dkx,Nparticles,nx/dd+1);
% [rspdm0,nk0,Angps0]=dmnk(UT',dx*dd,dkx,Nparticles,nx/dd+1);
% rspdm = rspdm/trace(rspdm);
% rspdm0 = rspdm0/trace(rspdm0);
% TG_rsp_fed(ii+1,1) = (trace(sqrtm(sqrtm(rspdm)*rspdm0*sqrtm(rspdm))))^2;
% T_d1(ii+1,1) = 1/2*sum(abs(eig(rspdm - rspdm0)));% trace distance
% Eig_TG = eig(rspdm);
% Coh_TG(ii+1) = max(Eig_TG);% coherence
% MomDen_TG(:,ii+1) = nk';
% MomDen_Fer(:,ii+1) = Angps';
%% Fermion
% rspdmFer=0;
% rspdmFerT=0;
% for jj = 1:Nparticles
%        rspdmFer =rspdmFer + phi(:,jj)*phi(:,jj)';
%        rspdmFerT =rspdmFerT + UT(:,jj)*UT(:,jj)';
% end
% rspdmFer = rspdmFer/trace(rspdmFer);
% rspdmFerT = rspdmFerT/trace(rspdmFerT);
% Fer_rsp_fed(ii+1,1) = (trace(sqrtm(sqrtm(rspdmFer)*rspdmFerT*sqrtm(rspdmFer))))^2;
% T_d2(ii+1,1) = 1/2*sum(abs(eig(rspdmFer - rspdmFerT)));% trace distance
% Eig_Fer = eig(rspdmFer);
% Coh_Fer(ii+1) = max(Eig_Fer);% coherence
end
toc;
% PP=P(1);
% for jj=dd:dd:nx
%     PP(jj/dd+1) = P(jj);
% end
%% Figure
% figure,
% plot(CylN,Fermion)
% xlabel('N');
% ylabel('F');

figure,
plot(CylTim,Fermion);
xlabel('t_f');
ylabel('F');
% axis([Ti TT 0.75 1]);
% set(gca,'LineWidth',2.0,'FontSize',20,'Fontname', 'Times New Roman');

% plot(CylTim,Coh_TG)
% xlabel('t_f');
% ylabel('Cohe');

% figure,
% plot(CylN,T_d)

% figure,
% plot(CylTim,Coh_Fer)
% xlabel('t_f');
% ylabel('Cohe');
% figure,
% plot(CylTim,E_dif,'LineWidth',2.0);
% xlabel('t_f');
% ylabel('\langle V \rangle');
% % axis([Ti TT 0 0.1]);
% set(gca,'LineWidth',2.0,'FontSize',20,'Fontname', 'Times New Roman');

% [t,p]=meshgrid(CylTim,PP);
% figure,
% mesh(p,t,MomDen_TG);
% view(0,90);
% figure,
% mesh(p,t,MomDen_Fer);
% view(0,90);
%  
