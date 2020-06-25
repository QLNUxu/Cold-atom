% close all;%Remove figure
clear;%Clear all define
% clc;%Clear command window
%%
tic;
m  = 1;
hb = 1;
a  = -20;
b  = 20;
L  = b-a;             %space lengh
dx = 10^-1/4 ;
nx = L / dx;
X  = a+L*(1:nx)/nx;                   %coordinates separation 
P  = (2*pi*hb/L)*[0:nx/2-1,-nx/2:-1]; %momentum separation 

%% bang-bang
syms a
gN    = 0.2;
C     = 0.407052;
w_0   = 20*2*pi;
w_i   = 20*2*pi/w_0;
w_f   = 5*2*pi/w_0; %intial and final frequency 
% a_i   = double(vpasolve(a*w_i^2-1/(a^3)-C*gN/(a^2*sqrt(pi*2))==0,a,[0,inf]));%intial width
% a_f   = double(vpasolve(a*w_f^2-1/(a^3)-C*gN/(a^2*sqrt(pi*2))==0,a,[0,inf]));% final width
% w_ix  = 1/a_i^2;
% w_fx  = 1/a_f^2;
% beta  = a_i; 
% gamma = a_f;
% d     = w_i^2; % delta
% C1    = -d*beta^2 + 1/(beta^2)+2*C*gN/(sqrt(2*pi)*beta);
% C2    = d*gamma^2 + 1/(gamma^2)+2*C*gN/(sqrt(2*pi)*gamma);%inytegrated cofficienr
% x1B   = sqrt((C2-C1)/(2*d)); % xB
% fun1  = @(s)sqrt(1./(C1+d.*s.^2-1./(s.^2)-2*C*gN./(sqrt(2*pi).*s)));
% t1    = integral(fun1,beta,x1B); % time for first segment
% fun2  = @(s)sqrt(1./(C2-d.*s.^2-1./(s.^2)-2*C*gN./(sqrt(2*pi).*s)));
% t2    = integral(fun2,x1B,gamma); % time for second segment
% T     = t1+t2;                       % total time
T     = 10;
nt    = 1 * 10^3 ;
dt    = T / nt;
% nt    = T / dt;
% delta = w_i^2;
%% inverse engineering
% syms a(t) w2(t) t
% a(t) = a_f + (6*(a_i-a_f)*t^5)/T^5-(15*(a_i-a_f)*t^4)/T^4+(10*(a_i-a_f)*t^3)/T^3;
% w2(t) = (1/a(t)^3 + gN*C/(sqrt(2*pi)*a(t)^2) - diff(a(t),2))/a(t);
%% Eigenvalue and Eigenstate
% operator_a    = diag(sqrt(1:nx-1),1);                   % creation operator
% operator_adag = operator_a';                            % annihilate operator
% operator_x    = 1/sqrt(2)*(operator_adag + operator_a);
% operator_p    = 1i/sqrt(2)*(operator_adag - operator_a);
% operator_H    = operator_p^2/(2*m) + m*w_fx^2*operator_x^2/2;
eigval_States = 10 ;                            % How many state
% Eig_Val       = eigs(operator_H,eigval_States,'smallestabs');
% V_i   = m*w_fx^2*X.^2/2;                 % harmonic oscillator
% V_f   = m*w_ix^2*X.^2/2;
% [ EigStai.H, EigStai.E, EigStai.U ] = H_eigen( V_i, eigval_States, X, dx, nx );
% [ EigStaf.H, EigStaf.E, EigStaf.U ] = H_eigen( V_f, eigval_States, X, dx, nx );

an      = 3;
K       = P.^2/(2*m);                % kinetic energy
V       = m*w_f^2*X.^2/2;
EigStat = zeros(nx,eigval_States);
E_eig   = zeros(eigval_States,1);
for n = 1:eigval_States
  % n     = 3;
  psi_0   = sqrt(1/(sqrt(pi)*2^(n-1)*factorial(n-1)*an))*exp(-X.^2/(2*an^2)).*hermiteH(n-1,X/an);
%    psi   = sqrt(2/an)*cos(n*pi*X/an);
%   [psix,mu] = ImgTimEvolution(psi,dt,dx,gN,X,P,m,V_i);
  [psix,psi_old,mu,mu_old,mu_error,j] = ImgTimEvolutionx(psi_0,dt,dx,gN,X,P,m,V,n,EigStat);
  EigStat(:,n) = psix';
  E_eig(n,1) = sum(conj(EigStat(:,n)).*ifft(K'.*fft(EigStat(:,n)))).*dx+sum(conj(EigStat(:,n)).*(m*w_f^2*X.^2/2)'.*EigStat(:,n)).*dx+sum(gN*abs(EigStat(:,n)).^4).*dx;
end

% plot(X,abs(EigStat(:,1)).^2,X,abs(EigStat(:,2)).^2)
% plot(X,abs(EigStai.U(:,1)).^2,X,abs(EigStat(:,1)).^2)
% sum(abs(EigStai.U(:,1).^2))*dx


% a1=abs(sum(conj(EigStat(:,1)).*psi')*dx)
% a2=abs(sum(conj(EigStat(:,2)).*psi')*dx)
% psi1=psi-a1*EigStat(:,1)'-a2*EigStat(:,2)';
% psi1=psi1/sqrt((dx*norm(psi1).^2));
% ImgTimEvolution(psi1,dt,dx,gN,X,P,m,V_i);
% plot(X,abs(psix).^2)
%% 
% K   = P.^2/(2*m);                % kinetic energy
UK  = exp(-1i*K/hb*dt/2);        % propagator in momentum space
psi = EigStat(:,eigval_States);
Nframe = 25;
spacetime=[];
% betacold = 0.5;                  % inverse temperature
% Proba = Probability( betacold, Eig_Val(1:eigval_States), eigval_States );

% Etx   = zeros(nt+1,1);
% w     = zeros(nt+1,1);
% x     = linspace(w_i,w_f,nt+1);
for j = 1 : nt

%     if (j*(T/nt))==0
%         w(j+1,1) = w_f^2;
%     elseif (j*(T/nt)) < t2
%         w(j+1,1) = delta;
%     elseif (j*(T/nt) > t2) && (j*(T/nt) < T)
%         w(j+1,1) = -delta;
%     elseif (j*(T/nt)) == T
%         w(j+1,1) = delta;
%     end    
%     
%     w(j+1,1) = w_fx^2 - (w_fx^2 - w_ix^2)*(j/nt);   % linear
%     w(j+1,1) = double(w2(j/nt*T));                  % inverse engineering
%     w(j+1,1) = (w_fx/(1-(w_ix-w_fx)*(j*(T/nt))/(T*w_ix)))^2;    % quasi-adiabatic
%     w(j+1,1) = 

    V        = m*w_f^2*X.^2/2 + gN*abs(psi(:,1)').^2;
    UV       = exp(-1i*V/hb*dt);
    psi(:,1) = FFT( psi( :, 1 )', UV, UK )';
    if mod(j,nt/Nframe) == 0 %Save wavefunction every Nframe steps
      spacetime=vertcat(spacetime,abs(psi'.^2));
    end
%     Et = zeros(eigval_States,1);
%     for k = 1:eigval_States
%         V        = m*w_f*X.^2/2 + gN*abs(psi(:,k)').^2;
%         UV       = exp(-1i*V/hb*dt);
%         psi(:,k) = FFT( psi( :, k )', UV, UK )';
%         Et(k,1)  = Proba(k,1)*(sum(conj(psi(:,k)).*ifft(K'.*fft(psi(:,k)))).*dx+sum(conj(psi(:,k)).*V'.*psi(:,k)).*dx);%+sum(gN*abs(psi(:,k)).^4).*dx);
%     end
%     Etx(j+1,1) = sum(Et);
%     at(j+1,1)  = sqrt(sum(conj(psi(:,1)).*X'.^2.*psi(:,1))*dx);
end
%%
subplot(1,3,1); %Plot potential
plot(X,V,'k'); xlabel('x (m)'); ylabel('V (J/hbar)');
subplot(1,3,2); %Plot initial and final density
plot(X,abs(psi_0).^2,'k',X,abs(psi).^2,'b');
legend('\psi(x,0)','\psi(x,T)');xlabel('x (m)');ylabel('|\psi|^2 (m^{-1})');
subplot(1,3,3); % Plot spacetime evolution as pcolor plot
dt_large=dt*double(nt/Nframe);
pcolor(X,dt_large*(1:1:Nframe),spacetime); shading interp;
xlabel('(m)'); ylabel('t (s)');
% Ei      = 0;         % initial energy
% Ef      = 0;         % final energy
% Etau    = 0;         % target energy
% rho_tau = 0;         % target state
% rho_tf  = 0;         % final state
% for i = 1:eigval_States
%     Ei   = Ei + Proba(i,1)*(sum(conj(EigStai.U(:,i)).*ifft(K'.*fft(EigStai.U(:,i)))).*dx+sum(conj(EigStai.U(:,i)).*(m*w_f^2*X.^2/2)'.*EigStai.U(:,i)).*dx+sum(gN*abs(EigStai.U(:,i)).^4).*dx);
%     Etau = Etau + Proba(i,1)*(sum(conj(EigStaf.U(:,i)).*ifft(K'.*fft(EigStaf.U(:,i)))).*dx+sum(conj(EigStaf.U(:,i)).*(m*w_i^2*X.^2/2)'.*EigStaf.U(:,i)).*dx+sum(gN*abs(EigStaf.U(:,i)).^4).*dx);
%     Ef   = Ef + Proba(i,1)*(sum(conj(psi(:,i)).*ifft(K'.*fft(psi(:,i)))).*dx+sum(conj(psi(:,i)).*V'.*psi(:,i)).*dx+sum(gN*abs(psi(:,i)).^4).*dx);
%     rho_tau = rho_tau + Proba(i)*EigStaf.U(:,i)*EigStaf.U(:,i)';
%     rho_tf  = rho_tf  + Proba(i)*psi(:,i)*psi(:,i)';
% end

% Fide        = trace(sqrtm(sqrtm(rho_tau)*rho_tf*sqrtm(rho_tau))) * dx ;
% F           = Fidelity( psi, EigStaf.U, eigval_States, dx  );
% W_exp       = Etx(end) - Etx(1);
% W_ad        = Etau - Etx(1);
% W_irr       = W_exp - W_ad;
%%
% time = linspace(0,T,nt+1);
% figure,
% plot(time,w)
% xlabel('t');
% ylabel('\omega^2(t)');
% figure,
% plot(time,Etx-Etx(1))
% xlabel('t');
% ylabel('\langle W(t) \rangle');
% text(1.25,23,'(b)')
toc