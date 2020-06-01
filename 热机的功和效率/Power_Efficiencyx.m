close all;%Remove figure
clear;%Clear all define
% clc;%Clear command window
tic
%%
% global Control Timeslots;

m  = 1;
hb = 1;
a  = -40;
b  = 40;
L  = b-a;             %space lengh
dx = 10^-1/4;
nx = L / dx;
X  = a+L*(1:nx)/nx;                   %coordinates separation 
P  = (2*pi*hb/L)*[0:nx/2-1,-nx/2:-1]; %momentum separation 

Control.m  = m;
Control.hb = hb;
Control.a  = a;
Control.b  = b;
Control.L  = L;
Control.nx = nx;
Control.dx = dx;
Control.X  = X;
Control.P  = P;

jj = 1;
% F_c   = zeros();
% F_e   = zeros();
% Power = zeros();
% eta   = zeros();
Tf    = zeros();
for gN = -2:0.1:2
syms a
% gN    = 1;
w_0   = 20*2*pi;
w_i   = 20*2*pi/w_0;
w_f   = 5*2*pi/w_0;                      % intial and final frequency 
C     = 0.407052;%0.168088;           % For the 10th or 100th eigenstate as ansatz
a_i   = double(vpasolve(a*w_i^2-1/(a^3)-C*gN/(a^2*sqrt(pi*2))==0,a,[0,inf]));%intial width
a_f   = double(vpasolve(a*w_f^2-1/(a^3)-C*gN/(a^2*sqrt(pi*2))==0,a,[0,inf]));% final width
w_ix  = 1/a_i^2;
w_fx  = 1/a_f^2;
beta  = a_i; 
gamma = a_f;
d     = w_i^2; % delta
C1    = -d*beta^2 + 1/(beta^2)+2*C*gN/(sqrt(2*pi)*beta);
C2    = d*gamma^2 + 1/(gamma^2)+2*C*gN/(sqrt(2*pi)*gamma);%integrated cofficienr
x1B   = sqrt((C2-C1)/(2*d)); % xB
fun1  = @(s)sqrt(1./(C1+d.*s.^2-1./(s.^2)-2*C*gN./(sqrt(2*pi).*s)));
t1    = integral(fun1,beta,x1B); % time for first segment
fun2  = @(s)sqrt(1./(C2-d.*s.^2-1./(s.^2)-2*C*gN./(sqrt(2*pi).*s)));
t2    = integral(fun2,x1B,gamma); % time for second segment
T     = t1+t2;                       % total time
% T     = 40;
betacold  = 0.5;                         % cold inverse temperature
betahot   = 0.05;                        % hot

phys.gN    = gN;
phys.C     = C;
phys.w_i   = w_i;
phys.w_f   = w_f;
phys.a_i   = a_i;
phys.a_f   = a_f;
phys.w_ix  = w_ix;
phys.w_fx  = w_fx;
phys.state = 100;                % How many states
phys.delta = d;
phys.betahot  = betahot;
phys.betacold = betacold;

Timeslots.t1  = t1;
Timeslots.t2  = t2;
Timeslots.T   = T;
Timeslots.nt  = 3 * 10^3 ;
% Timeslots.dt  = 10^-3;
Timeslots.dt  = T / Timeslots.nt;

%%
%Compression
[ Output_cmp ]=Compression(Control, Timeslots, phys);
%Expansion
[ Output_exp ]=Expansion(Control, Timeslots, phys);

F.cmp(jj,1)      = Output_cmp.Fide;
F.exp(jj,1)      = Output_exp.Fide;
F.c(:,jj)        = Output_cmp.F;
F.e(:,jj)        = Output_exp.F;
E.cmp_Ei(jj,1)   = Output_cmp.Ei;
E.cmp_Ef(jj,1)   = Output_cmp.Ef;
E.exp_Ei(jj,1)   = Output_exp.Ei;
E.exp_Ef(jj,1)   = Output_exp.Ef;
W.cmp(jj,1)      = Output_cmp.W_cmp;
W.exp(jj,1)      = Output_exp.W_exp;
W.cmp_ad(jj,1)   = Output_cmp.W_ad;
W.exp_ad(jj,1)   = Output_exp.W_ad;
W.cmp_irr(jj,1)  = Output_cmp.W_irr;
W.exp_irr(jj,1)  = Output_exp.W_irr;
W.Power(jj,1) = -(Output_cmp.W_cmp + Output_exp.W_exp) ./ (2 * T);
W.eta(jj,1)   = -(Output_cmp.W_cmp + Output_exp.W_exp) ./ (Output_exp.Ei - Output_cmp.Ef);
Tf(jj,1)    = T;

jj = jj+1
end
%%
g = linspace(-2,2,jj-1)';
% figure,
% yyaxis left
% plot(g,W.Power);
% xlabel('gN');
% ylabel('P');
% % % axis([0 TT 0 4000]);
% % set(gca,'LineWidth',2.0,'FontSize',20,'Fontname', 'Times New Roman');
% yyaxis right
% plot(g,W.eta,'--');
% xlabel('gN');
% ylabel('\eta');
% % set(gca,'LineWidth',2.0,'FontSize',20,'Fontname', 'Times New Roman');
% 
% 
% vqF = interp1(g , F.cmp, linspace(-2,2,201),'pchip');
% vpW = interp1(g , W.cmp_irr, linspace(-2,2,201),'pchip');
% figure,
% yyaxis left
% % semilogy(g,1-F.cmp)
% semilogy(linspace(-2,2,201),1-vqF)
% xlabel('gN');
% ylabel('1-F');
% yyaxis right
% semilogy(linspace(-2,2,201),vpW)
% % semilogy(g,W.cmp_irr)
% xlabel('gN');
% ylabel('W_{irr}');
%%
runtime=toc/3600