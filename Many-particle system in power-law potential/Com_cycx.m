% close all;%Remove figure
clear;%Clear all define
% clc;%Clear command window
%% Preparation
tic;
% global X dx nx
m=1;
hb=1;
% boundary
a = -10;
b = 10;

L=b-a; %space lengh
nx = 1 * 10^3 ;
dx=L/nx;
X=a+L*(1:nx)/nx;%coordinates separation
P=(2*pi*hb/L)*[0:nx/2-1,-nx/2:-1];%%momentum separation 

Nparticles = 10 ;% How many states
q = 3 ;% Power for X

% compression
omegai = 1;
omegaf = 2 * omegai;
lambdai = omegai^(q+1);
lambdaf = omegaf^(q+1);

% nn = 10 ;% How many states
% CylN = cumsum(1:nn,1) ;% build x-axis for N
% Fermion = zeros(nn,1);

% parfor Nparticles = 1:nn

% Initial and target Ham. value state 
x0 = 0;
[ H_i, E_i, U_i ] = H_eigen( lambdai, x0, q, Nparticles, X, dx, nx );
[ H_T, E_T, U_T ] = H_eigen( lambdaf, x0, q, Nparticles, X, dx, nx );
%% Total cycling
Ti = 0.1 ;
TT = 4 ;
dT = 0.02 ;
nT = (TT-Ti)/dT;
CylTim = cumsum(Ti:dT:TT,1); % build x-axis for T

parfor ii=0:nT

% Define time
% T = 1 ;
T = Ti + ii * dT;
nt = 1 * 10^3;% number of time steps
dt = T/nt;
dtm = ones(1,nt);

psi=U_i;% define the initial state

%% Cycling

for j=0:nt
% ground state
    lambda=1/(q*T^5*gamma(1/2+q))*2^(-1-1/(2+2*q))*pi^(1/2-3/(4+4*q))*((4^(1/(1+q))*T^20)/(-6*(j*(T/nt))^5*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))+15*(j*(T/nt))^4*T*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))-10*(j*(T/nt))^3*T^2*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))+T^5*(1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q)))^3+120*pi^(1/(1+q))*(j*(T/nt))^3*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))-180*pi^(1/(1+q))*(j*(T/nt))^2*T*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))+60*pi^(1/(1+q))*(j*(T/nt))*T^2*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q))))*(1/T^5*2^(-(1/(2+2*q)))*pi^(1/(4+4*q))*(-6*(j*(T/nt))^5*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))+15*(j*(T/nt))^4*T*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))-10*(j*(T/nt))^3*T^2*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))+T^5*(1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))))^(1-2*q);
% q=2, n=10
%     lambda = (sqrt(7/221)*T^15*((180*(7/221)^(1/6)*(j*(T/nt))^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6)))/T^4+(120*(7/221)^(1/6)*(j*(T/nt))^3*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6)))/T^5+(60*(7/221)^(1/6)*(j*(T/nt))*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6)))/T^3+(sqrt(221/7)*T^15)/(6*(j*(T/nt))^5*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6))+T^5*(1/lambdai)^(1/6))^3))/(6*(j*(T/nt))^5*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6))+T^5*(1/lambdai)^(1/6))^3;
% q=2, n=20
%     lambda = (T^10*(-1640*3^(1/3)*29^(2/3)*(j*(T/nt))^3*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+2460*3^(1/3)*29^(2/3)*(j*(T/nt))^2*T*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+820*3^(1/3)*29^(2/3)*(j*(T/nt))*T^2*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6))+(841*41^(1/3)*T^20)/(6*(j*(T/nt))^5*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6))+T^5*(1/lambdai)^(1/6))^3))/(841*41^(1/3)*(6*(j*(T/nt))^5*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6))+T^5*(1/lambdai)^(1/6))^3);
% q=2, n=30
%     lambda = (T^10*(-2440*5583^(1/3)*(j*(T/nt))^3*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+3660*5583^(1/3)*(j*(T/nt))^2*T*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))-1220*5583^(1/3)*(j*(T/nt))*T^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+(1861*61^(1/3)*T^20)/(6*(j*(T/nt))^5*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6))+T^5*(1/lambdai)^(1/6))^3))/(1861*61^(1/3)*(6*(j*(T/nt))^5*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6))+T^5*(1/lambdai)^(1/6))^3);    
% q=2, n=40
%     lambda = (T^10*(-1080*3281^(1/3)*(j*(T/nt))^3*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+1620*3281^(1/3)*(j*(T/nt))^2*T*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+540*3281^(1/3)*(j*(T/nt))*T^2*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6))+(3281*T^20)/(6*(j*(T/nt))^5*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6))+T^5*(1/lambdai)^(1/6))^3))/(3281*(6*(j*(T/nt))^5*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6))+T^5*(1/lambdai)^(1/6))^3);
% q=2, n=50
%     lambda = (T^10*(-4040*15303^(1/3)*(j*(T/nt))^3*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+6060*15303^(1/3)*(j*(T/nt))^2*T*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))-2020*15303^(1/3)*(j*(T/nt))*T^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+(5101*101^(1/3)*T^20)/(6*(j*(T/nt))^5*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6))+T^5*(1/lambdai)^(1/6))^3))/(5101*101^(1/3)*(6*(j*(T/nt))^5*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/6)-(1/lambdai)^(1/6))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/6)+(1/lambdai)^(1/6))+T^5*(1/lambdai)^(1/6))^3);
% q=3, n=10
%     lambda = (T^20*(24*sqrt(3345)*(j*(T/nt))^2*T*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+16*sqrt(3345)*(j*(T/nt))^3*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+8*sqrt(3345)*(j*(T/nt))*T^2*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+(223*T^20)/(6*(j*(T/nt))^5*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+T^5*(1/lambdai)^(1/8))^3))/(223*(6*(j*(T/nt))^5*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+T^5*(1/lambdai)^(1/8))^5);
% q=3, n=20
%     lambda = (T^20*(24*sqrt(1405)*(j*(T/nt))^2*T*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+16*sqrt(1405)*(j*(T/nt))^3*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+8*sqrt(1405)*(j*(T/nt))*T^2*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+(281*T^20)/(6*(j*(T/nt))^5*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+T^5*(1/lambdai)^(1/8))^3))/(281*(6*(j*(T/nt))^5*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+T^5*(1/lambdai)^(1/8))^5);    
% q=3, n=30
%     lambda = (T^20*(24*sqrt(115)*(j*(T/nt))^2*T*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+16*sqrt(115)*(j*(T/nt))^3*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+8*sqrt(115)*(j*(T/nt))*T^2*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+(69*sqrt(3)*T^20)/(6*(j*(T/nt))^5*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+T^5*(1/lambdai)^(1/8))^3))/(69*sqrt(3)*(6*(j*(T/nt))^5*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+T^5*(1/lambdai)^(1/8))^5);
% q=3, n=40
%     lambda = (T^20*(24*sqrt(1005)*(j*(T/nt))^2*T*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+16*sqrt(1005)*(j*(T/nt))^3*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+8*sqrt(1005)*(j*(T/nt))*T^2*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+(469*T^20)/(6*(j*(T/nt))^5*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+T^5*(1/lambdai)^(1/8))^3))/(469*(6*(j*(T/nt))^5*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+T^5*(1/lambdai)^(1/8))^5);
% q=3, n=50
%     lambda = (T^20*(24*sqrt(35)*(j*(T/nt))^2*T*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+16*sqrt(35)*(j*(T/nt))^3*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+8*sqrt(35)*(j*(T/nt))*T^2*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+(63*sqrt(3)*T^20)/(6*(j*(T/nt))^5*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+T^5*(1/lambdai)^(1/8))^3))/(63*sqrt(3)*(6*(j*(T/nt))^5*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+T^5*(1/lambdai)^(1/8))^5);
% reference polynomial
%     lambda = (6*(j*(T/nt))^5*(lambdaf-lambdai))/T^5-(15*(j*(T/nt))^4*(lambdaf-lambdai))/T^4+(10*(j*(T/nt))^3*(lambdaf-lambdai))/T^3+lambdai; 
% quasiadiabatic
%     lambda = (omegai/(1-(omegaf-omegai)*(j*(T/nt))/(T*omegaf)))^(q+1);
    
    K=P.^2/(2*m);% kinetic energy
    UK=exp(-1i*K/hb*dt/2);% propagator in momentum space
    V=X.^(2*q)*lambda/2;% potential energy
%     V_non=X.^(2*nn)*lambda_non/2;
    UV=exp(-1i*V/hb*dt);% propagator in coordinate space
%     UV_non=exp(-1i*V_non/hb*dt);
    
    psi = Dyn_STA( psi, UV, UK, Nparticles, nx );
%     psi_Non = Dyn_STA( psi_Non, UV_non, UK, n);

end

%% Fidelity
% Fermion(Nparticles,1) = abs(det(conj(U_T)'*psi*dx))^2;% particle cycle
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
dkx=1;
[rspdm]=dmnk(psi',dx,dkx,Nparticles,nx);
[rspdm0]=dmnk(U_T',dx,dkx,Nparticles,nx);
rspdm = rspdm/trace(rspdm);
rspdm0 = rspdm0/trace(rspdm0);

%% Fermion
rspdmFer=0;
rspdmFerT=0;
for jj = 1:Nparticles
       rspdmFer =rspdmFer + psi(:,jj)*psi(:,jj)';
       rspdmFerT =rspdmFerT + U_T(:,jj)*U_T(:,jj)';
end
rspdmFer = rspdmFer/trace(rspdmFer);
rspdmFerT = rspdmFerT/trace(rspdmFerT);

end
% PP=P(1);
% for jj=dd:dd:nx
%     PP(jj/dd+1) = P(jj);
% end
toc;
%% Figure

figure,
plot(CylN,Fermion)
xlabel('N');
ylabel('Fidelity');

% figure,
% plot(CylTim,Fermion);
% xlabel('t_f');
% ylabel('Fidelity');
% % axis([Ti TT 0.75 1]);
% set(gca,'LineWidth',2.0,'FontSize',20,'Fontname', 'Times New Roman');

