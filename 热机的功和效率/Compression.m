function [ Output ] = Compression( Control, Timeslots, phys )
m   = Control.m;
hb  = Control.hb;
dx  = Control.dx;
nx  = Control.nx;
X   = Control.X;
P   = Control.P;

nt  = Timeslots.nt;%number of steps
dt  = Timeslots.dt;
% t1  = Timeslots.t1;
% t2  = Timeslots.t2;
T   = Timeslots.T;

gN      = phys.gN;
w_i     = phys.w_f;
w_f     = phys.w_i;
a_i     = phys.a_f;
a_f     = phys.a_i;
w_ix    = phys.w_fx;
w_fx    = phys.w_ix;
delta   = phys.delta;
beta    = phys.betacold;
C       = phys.C;

%% Eigenvalue and Eigenstate
operator_a    = diag(sqrt(1:nx-1),1);                   % creation operator
operator_adag = operator_a';                            % annihilate operator
operator_x    = 1/sqrt(2)*(operator_adag + operator_a);
operator_p    = 1i/sqrt(2)*(operator_adag - operator_a);
operator_H    = operator_p^2/(2*m) + m*w_ix^2*operator_x^2/2;
eigval_States = phys.state ;                            % How many state
Eig_Val       = eigs(operator_H,eigval_States,'smallestabs');
V_i           = m*w_ix^2*X.^2/2;              % harmonic oscillator
V_f           = m*w_fx^2*X.^2/2;
[ EigStai.H, EigStai.E, EigStai.U ] = H_eigen( V_i, eigval_States, X, dx, nx );
[ EigStaf.H, EigStaf.E, EigStaf.U ] = H_eigen( V_f, eigval_States, X, dx, nx );
%% inverse engineering
% syms a(t) w2(t) t
% a(t) = a_i + (6*(a_f-a_i)*t^5)/T^5-(15*(a_f-a_i)*t^4)/T^4+(10*(a_f-a_i)*t^3)/T^3;
% w2(t) = (1/a(t)^3 + gN*C/(sqrt(2*pi)*a(t)^2) - diff(a(t),2))/a(t);
%% 
K     = P.^2/(2*m);                %kinetic energy
UK    = exp(-1i*K/hb*dt/2);        %propagator in momentum space
psi   = EigStai.U;
Proba = Probability( beta, Eig_Val, eigval_States );
% x     = linspace(w_i,w_f,nt+1);
w     = zeros(nt+1,1);
for j = 0 : nt
    % bang-bang
%     if (j*(T/nt))==0
%         w(j+1,1) = w_i^2;
%     elseif (j*(T/nt)) < t2
%         w(j+1,1) = delta;
%     elseif (j*(T/nt) > t2) && (j*(T/nt) < T)
%         w(j+1,1) = -delta;
%     elseif (j*(T/nt)) == T
%         w(j+1,1) = delta;
%     end    
    
%     w(j+1,1) = x(j+1);    % linear
%     w(j+1,1) = double(w2(j/nt*T));                  % inverse engineering 
%     w(j+1,1) = (w_i/(1-(w_f-w_i)*(j*(T/nt))/(T*w_f)))^2;
    for k = 1:eigval_States
        V        = m*w(j+1,1)*X.^2/2 + gN*abs(psi(:,k)').^2;
        UV       = exp(-1i*V/hb*dt);
        psi(:,k) = FFT( psi( :, k )', UV, UK )';
    end

end

Ei      = 0;         % initial energy
Ef      = 0;         % final energy
Etau    = 0;         % target energy
rho_tau = 0;         % target state
rho_tf  = 0;         % final state
for i = 1:eigval_States
    Ei   = Ei   + Proba(i,1)*(sum(conj(EigStai.U(:,i)).*ifft(K'.*fft(EigStai.U(:,i)))).*dx + sum(conj(EigStai.U(:,i)).*(m*w_i^2*X.^2/2)'.*EigStai.U(:,i)).*dx + sum(gN*abs(EigStai.U(:,i)).^4).*dx);
    Etau = Etau + Proba(i,1)*(sum(conj(EigStaf.U(:,i)).*ifft(K'.*fft(EigStaf.U(:,i)))).*dx + sum(conj(EigStaf.U(:,i)).*(m*w_f^2*X.^2/2)'.*EigStaf.U(:,i)).*dx + sum(gN*abs(EigStaf.U(:,i)).^4).*dx);
    Ef   = Ef   + Proba(i,1)*(sum(conj(psi(:,i)).*ifft(K'.*fft(psi(:,i)))).*dx             + sum(conj(psi(:,i)).*(m*w_f^2*X.^2/2)'.*psi(:,i)).*dx             + sum(gN*abs(psi(:,i)).^4).*dx);%+sum(gN*abs(psi(:,i)).^4).*dx);
    rho_tau = rho_tau + Proba(i)*EigStaf.U(:,i)*EigStaf.U(:,i)';
    rho_tf  = rho_tf  + Proba(i)*psi(:,i)*psi(:,i)';
end

Fide        = trace(sqrtm(sqrtm(rho_tau)*rho_tf*sqrtm(rho_tau))) * dx ;
F           = Fidelity( psi, EigStaf.U, eigval_States, dx  );
W_cmp       = Ef - Ei;
W_ad        = Etau - Ei;
W_irr       = W_cmp - W_ad;
%% output
Output.Ei       = Ei;
Output.Ef       = Ef;
Output.Etau     = Etau;
Output.Fide     = Fide;
Output.F        = F;
Output.W_cmp    = W_cmp;
Output.W_ad     = W_ad;
Output.W_irr    = W_irr;
Output.w        = w;
end
