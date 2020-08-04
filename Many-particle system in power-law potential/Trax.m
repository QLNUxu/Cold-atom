% close all; % Remove figure
clear; % Clear all define
% clc; % Clear command window
tic
%% Preparation
% global X dx nx
m = 1; % mass
hb = 1;% Plank
% Boundary
a = -10 ;% left boundary
b = 20 ;% right boundary
L = b-a; % space lengh
nx = 1 * 10^3 ;
dx = L/nx;
X = a+L*(1:nx)/nx; % coordinates separation
P = (2*pi*hb/L)*[0:nx/2-1,-nx/2:-1]; % momentum separation 

% Define time
T = 10 ;
nt = 1 * 10^3 ; % number of steps
dt = T/nt;
dtm = ones(1,nt);
t_unitary = [0, cumsum(dtm(1 : end))/length(dtm)*T];

%% Define the physics

Nparticles = 10 ;
q = 2 ;% Powe

% transport
xii = 0 ; % initial position
xif = 5 ; % final position
lambda = 1;
% a = 1/(3^(1/6));% q=2, n=0
a = (19/543)^(1/6);% q=2, n=10
% a = (101/15303)^(1/6);% q=2, n=50

% Diagonal H
[ H_i, E_i, U_i ] = H_eigen( lambda, xii, q, Nparticles, X, dx, nx ); % initial state
[ H_T, E_T, U_T ] = H_eigen( lambda, xif, q, Nparticles, X, dx, nx ); % final state
%% Define the matrix
% Define the initial state
psi = U_i;

for j=0:nt
% q=1
    x0 = (120*(j*(T/nt))^3*xif-180*(j*(T/nt))^2*T*xif+60*(j*(T/nt))*T^2*xif+6*(j*(T/nt))^5*lambda*xif-15*(j*(T/nt))^4*T*lambda*xif+10*(j*(T/nt))^3*T^2*lambda*xif-120*(j*(T/nt))^3*xii+180*(j*(T/nt))^2*T*xii-60*(j*(T/nt))*T^2*xii-6*(j*(T/nt))^5*lambda*xii+15*(j*(T/nt))^4*T*lambda*xii-10*(j*(T/nt))^3*T^2*lambda*xii+T^5*lambda*xii)/(T^5*lambda);
% q=2, n=0
%     x0 = (6*(j*(T/nt))^5*(xif-xii))/T^5-(15*(j*(T/nt))^4*(xif-xii))/T^4+(10*(j*(T/nt))^3*(xif-xii))/T^3+xii+(a^2*lambda)/(2^(1/3)*((-120*(j*(T/nt))^3*lambda^2*xif+180*(j*(T/nt))^2*T*lambda^2*xif-60*(j*(T/nt))*T^2*lambda^2*xif+sqrt(2)*T^5*sqrt((lambda^4*(a^6*T^10*lambda^2+1800*(j*(T/nt))^2*(2*(j*(T/nt))^2-3*(j*(T/nt))*T+T^2)^2*(xif-xii)^2))/T^10)+120*(j*(T/nt))^3*lambda^2*xii-180*(j*(T/nt))^2*T*lambda^2*xii+60*(j*(T/nt))*T^2*lambda^2*xii)/T^5)^(1/3))-((-120*(j*(T/nt))^3*lambda^2*xif+180*(j*(T/nt))^2*T*lambda^2*xif-60*(j*(T/nt))*T^2*lambda^2*xif+sqrt(2)*T^5*sqrt((lambda^4*(a^6*T^10*lambda^2+1800*(j*(T/nt))^2*(2*(j*(T/nt))^2-3*(j*(T/nt))*T+T^2)^2*(xif-xii)^2))/T^10)+120*(j*(T/nt))^3*lambda^2*xii-180*(j*(T/nt))^2*T*lambda^2*xii+60*(j*(T/nt))*T^2*lambda^2*xii)/T^5)^(1/3)/(2^(2/3)*lambda);
% q=2, n=10
%     x0 = (6*(j*(T/nt))^5*(xif-xii))/T^5-(15*(j*(T/nt))^4*(xif-xii))/T^4+(10*(j*(T/nt))^3*(xif-xii))/T^3+xii+(19*a^2*lambda)/(2^(1/3)*((-120*(j*(T/nt))^3*lambda^2*xif+180*(j*(T/nt))^2*T*lambda^2*xif-60*(j*(T/nt))*T^2*lambda^2*xif+sqrt(2)*T^5*sqrt((lambda^4*(6859*a^6*T^10*lambda^2+1800*(j*(T/nt))^2*(2*(j*(T/nt))^2-3*(j*(T/nt))*T+T^2)^2*(xif-xii)^2))/T^10)+120*(j*(T/nt))^3*lambda^2*xii-180*(j*(T/nt))^2*T*lambda^2*xii+60*(j*(T/nt))*T^2*lambda^2*xii)/T^5)^(1/3))-((-120*(j*(T/nt))^3*lambda^2*xif+180*(j*(T/nt))^2*T*lambda^2*xif-60*(j*(T/nt))*T^2*lambda^2*xif+sqrt(2)*T^5*sqrt((lambda^4*(6859*a^6*T^10*lambda^2+1800*(j*(T/nt))^2*(2*(j*(T/nt))^2-3*(j*(T/nt))*T+T^2)^2*(xif-xii)^2))/T^10)+120*(j*(T/nt))^3*lambda^2*xii-180*(j*(T/nt))^2*T*lambda^2*xii+60*(j*(T/nt))*T^2*lambda^2*xii)/T^5)^(1/3)/(2^(2/3)*lambda);
% q=2, n=50
%     x0 = (6*(j*(T/nt))^5*(xif-xii))/T^5-(15*(j*(T/nt))^4*(xif-xii))/T^4+(10*(j*(T/nt))^3*(xif-xii))/T^3+xii+(101*a^2*lambda)/(2^(1/3)*((-120*(j*(T/nt))^3*lambda^2*xif+180*(j*(T/nt))^2*T*lambda^2*xif-60*(j*(T/nt))*T^2*lambda^2*xif+sqrt(2)*T^5*sqrt((lambda^4*(1030301*a^6*T^10*lambda^2+1800*(j*(T/nt))^2*(2*(j*(T/nt))^2-3*(j*(T/nt))*T+T^2)^2*(xif-xii)^2))/T^10)+120*(j*(T/nt))^3*lambda^2*xii-180*(j*(T/nt))^2*T*lambda^2*xii+60*(j*(T/nt))*T^2*lambda^2*xii)/T^5)^(1/3))-((-120*(j*(T/nt))^3*lambda^2*xif+180*(j*(T/nt))^2*T*lambda^2*xif-60*(j*(T/nt))*T^2*lambda^2*xif+sqrt(2)*T^5*sqrt((lambda^4*(1030301*a^6*T^10*lambda^2+1800*(j*(T/nt))^2*(2*(j*(T/nt))^2-3*(j*(T/nt))*T+T^2)^2*(xif-xii)^2))/T^10)+120*(j*(T/nt))^3*lambda^2*xii-180*(j*(T/nt))^2*T*lambda^2*xii+60*(j*(T/nt))*T^2*lambda^2*xii)/T^5)^(1/3)/(2^(2/3)*lambda);
% reference polynomial
%     x0 = (6*(j*(T/nt))^5*(xif-xii))/T^5-(15*(j*(T/nt))^4*(xif-xii))/T^4+(10*(j*(T/nt))^3*(xif-xii))/T^3+xii; 
% reference linear
%     x0 = xii-((j*(T/nt))*(-xif+xii))/T;
% couping
%     lambda = -((2*(j*(T/nt))^2*(2*epsilon-lambdaf-lambdai))/T^2)+lambdai-((j*(T/nt))*(-4*epsilon+lambdaf+3*lambdai))/T;
%     x0 = Expression1(j+1);
    
    K=P.^2/(2*m); % kinetic energy
    UK=exp(-1i*K/hb*dt/2); % propagator in momentum space
    V=(X-x0).^(2*q)*lambda/2; % potential energy
%     V_non=X.^(2*nn)*lambda_non/2;
    UV=exp(-1i*V/hb*dt); % propagator in coordinate space
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
subplot(2,3,3),
pcolor(x1,x2,rspdmTGT);
shading interp
axis([2 8 2 8]);
% title('\lambda^{n=9}(t)');
% xlabel(' x ');
% ylabel(' x^\prime ');
set(gca,'LineWidth',4.0,'FontSize',55,'Fontname', 'Times New Roman');
% set(gca,'XTick',-1.5: 1.5 : 1.5);
% set(gca,'XTicklabel',{'-1.5','0','1.5'});
% set(gca,'YTick',-1.5: 1.5 : 1.5);
% set(gca,'YTicklabel',{'-1.5','0','1.5'});
set(gca,'XTick',[]);
set(gca,'XTicklabel',{});
set(gca,'YTick',[]);
set(gca,'YTicklabel',{});
hold on
contour(x1,x2,real(rspdmTGT),3,'w','LineWidth',4);
contour(x1,x2,real(rspdmTG),3,'r--','LineWidth',4);
%% Figure
% figure,
% plot(t_unitary,W_STA);
% axis([a b 0 1]);
% xlabel(' x ');
% ylabel(' |\psi|^2 ');
% set(gca,'LineWidth',2.0,'FontSize',20,'Fontname', 'Times New Roman','xticklabel',[]);

% [x1,x2] = meshgrid(XX,XX);
% figure,
% mesh(x1,x2,real(rspdm));
% view(90,-90);
% hold on
% contour(x1,x2,real(rspdm),[2,3,4],'k','LineWidth',1);
% hold on
% contour(x1,x2,real(rspdm0),[2,3,4],'w','LineWidth',1);
% axis([0 10 0 10]);
% set(gca,'LineWidth',2.0,'FontSize',22,'Fontname', 'Times New Roman');
% h=text(2.5, -2.5, -2,'(a)','Color','white','FontSize',20);

%  t1=0:dt*dStep:T;
%  [t,x]=meshgrid(t1,X);
%  figure,
%  mesh(x,t,abs(PP1).^2);
%  view(0,90);
% 
%  t1=0:dt*dStep:T;
%  [t,p]=meshgrid(t1,PP);
%  figure,
%  mesh(p,t,QQ1);
%  view(0,90);
toc