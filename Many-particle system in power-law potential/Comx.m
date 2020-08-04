% close all; % Remove figure
clear; % Clear all define
% clc; % Clear command window
tic
%% Preparation

% global X dx nx
m = 1; % mass
hb = 1;
% Boundary
a = -10 ;
b = 10 ;

L = b - a; % space lengh
nx = 1 * 10^3 ;
dx = L/nx;
X = a+L*(1:nx)/nx; % coordinates separation
P = (2*pi*hb/L)*[0:nx/2-1,-nx/2:-1];%(2*pi*hb/L)*[0:nx/2-1,-nx/2:-1]; % momentum separation 

% Define time
T = 0.5 ;
nt = 1 * 10^3; % number of steps
dt = T/nt;
dtm = ones(1,nt);
t_unitary = [0, cumsum(dtm(1 : end))/length(dtm)*T];

%% Define the physics

% N = 1 ;% particle number
Nparticles = 10 ;
q = 3 ;% Power for X

% from gamma to a
omega0 = 1;
omegaf = 2 * omega0;
lambdai = omega0^(q+1);
lambdaf = omegaf^(q+1);

% Diagonal H
x0=0;
[ H_i, E_i, U_i ] = H_eigen( lambdai, x0, q, Nparticles, X, dx, nx );
[ H_T, E_T, U_T ] = H_eigen( lambdaf, x0, q, Nparticles, X, dx, nx );

% Define the initial state
psi = U_i;
 
for j=0:nt
    
% ground state
%     lambda = 1/(q*T^5*gamma(1/2+q))*2^(-1-1/(2+2*q))*pi^(1/2-3/(4+4*q))*((4^(1/(1+q))*T^20)/(-6*(j*(T/nt))^5*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))+15*(j*(T/nt))^4*T*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))-10*(j*(T/nt))^3*T^2*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))+T^5*(1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q)))^3+120*pi^(1/(1+q))*(j*(T/nt))^3*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))-180*pi^(1/(1+q))*(j*(T/nt))^2*T*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))+60*pi^(1/(1+q))*(j*(T/nt))*T^2*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q))))*(1/T^5*2^(-(1/(2+2*q)))*pi^(1/(4+4*q))*(-6*(j*(T/nt))^5*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))+15*(j*(T/nt))^4*T*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))-10*(j*(T/nt))^3*T^2*((1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))-(1/(q*lambdaf*gamma(1/2+q)))^(1/(2+2*q)))+T^5*(1/(q*lambdai*gamma(1/2+q)))^(1/(2+2*q))))^(1-2*q);
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
    lambda = (T^20*(24*sqrt(3345)*(j*(T/nt))^2*T*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+16*sqrt(3345)*(j*(T/nt))^3*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+8*sqrt(3345)*(j*(T/nt))*T^2*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+(223*T^20)/(6*(j*(T/nt))^5*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+T^5*(1/lambdai)^(1/8))^3))/(223*(6*(j*(T/nt))^5*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+10*(j*(T/nt))^3*T^2*((1/lambdaf)^(1/8)-(1/lambdai)^(1/8))+15*(j*(T/nt))^4*T*(-(1/lambdaf)^(1/8)+(1/lambdai)^(1/8))+T^5*(1/lambdai)^(1/8))^5);
% reference polynomial
%     lambda = (6*(j*(T/nt))^5*(lambdaf-lambdai))/T^5-(15*(j*(T/nt))^4*(lambdaf-lambdai))/T^4+(10*(j*(T/nt))^3*(lambdaf-lambdai))/T^3+lambdai; 
%     lambda = lambdai - (lambdai - lambdaf)*(j*(T/nt))/T;
% quasiadiabatic
%     lambda = (omega0/(1-(omegaf-omega0)*(j*(T/nt))/(T*omegaf)))^(q+1);

    K=P.^2/(2*m);% kinetic energy
    UK=exp(-1i*K/hb*dt/2);% propagator in momentum space
    V=X.^(2*q)*lambda/2;% potential energy
%     V_non=X.^(2*nn)*lambda_non/2;
    UV=exp(-1i*V/hb*dt);% propagator in coordinate space
%     UV_non=exp(-1i*V_non/hb*dt);
    
    psi = Dyn_STA( psi, UV, UK, Nparticles, nx );
end

%% RSPDM
dkx=1;
[rspdmTGT]=dmnk(U_T',dx,dkx,Nparticles,nx);
rspdmTGT = rspdmTGT/trace(rspdmTGT);
[rspdmTG]=dmnk(psi',dx,dkx,Nparticles,nx);
rspdmTG = rspdmTG/trace(rspdmTG);

[x1,x2] = meshgrid(X,X);
% figure,
subplot(2,3,4),
pcolor(x1,x2,rspdmTGT);
shading interp
axis([-1.7 1.7 -1.7 1.7]);
% title('\lambda^{n=9}(t)');
xlabel(' x ');
ylabel(' x^\prime ');
set(gca,'LineWidth',2.0,'FontSize',22,'Fontname', 'Times New Roman');
set(gca,'XTick',-1.5: 1.5 : 1.5);
set(gca,'XTicklabel',{'-1.5','0','1.5'});
set(gca,'YTick',-1.5: 1.5 : 1.5);
set(gca,'YTicklabel',{'-1.5','0','1.5'});
% set(gca,'XTick',[]);
% set(gca,'XTicklabel',{});
% set(gca,'YTick',[]);
% set(gca,'YTicklabel',{});
hold on
contour(x1,x2,real(rspdmTGT),3,'w','LineWidth',2);
contour(x1,x2,real(rspdmTG),3,'r--','LineWidth',2);
%% Figure

% wavefunction or fidelity
% figure,
% plot(X,conj(psi_Gau(:,1)').*ifft(K.*fft(psi_Gau(:,1)')))
% axis([a b 0 1.2]);
% axis([a b 0 1]);
% xlabel(' x ');
% ylabel(' |\psi|^2 ');
% set(gca,'LineWidth',2.0,'FontSize',20,'Fontname', 'Times New Roman','xticklabel',[]);

% reduced single particle density matrix

% [x1,x2] = meshgrid(XX,XX);
% figure,
% mesh(x1,x2,real(rspdmTG));
% view(90,-90);
% hold on
% contour(x1,x2,real(rspdmTG),2,'k','LineWidth',1.5);
% hold on
% contour(x1,x2,real(rspdmTGT),2,'w','LineWidth',1.5);
% axis([-3.5 3.5 -3.5 3.5]);
% set(gca,'LineWidth',2.0,'FontSize',22,'Fontname', 'Times New Roman');
% % h=text(2.5, -2.5, '(a)','Color','white','FontSize',20);

% 
% [x1,x2] = meshgrid(XX,XX);
% figure,
% mesh(x1,x2,real(rspdmFer));
% view(90,-90);
% hold on
% contour(x1,x2,real(rspdmFer),3,'k','LineWidth',1.5);
% hold on
% contour(x1,x2,real(rspdmFerT),3,'w','LineWidth',1.5);
% axis([-3.5 3.5 -3.5 3.5]);
% set(gca,'LineWidth',2.0,'FontSize',22,'Fontname', 'Times New Roman');
% h=text(2.5, -2.5, '(a)','Color','white','FontSize',20);

toc