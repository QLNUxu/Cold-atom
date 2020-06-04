% SOLVES THE 1D GPE IN IMAGINARY TIME USING THE SPLIT-STEP METHOD
function [psi,psi_old,mu,mu_old,mu_error,j] = ImgTimEvolutionx(psi,dt,dx,g1d,x,k,m,V,n,EigStat)
hbar=1; 
% dx=x(2)-x(1); 
% dk=2*pi/(x(end)-x(1));
N=dx*norm(psi).^2; 
Nx=length(x); 
psi_old=psi;
mu_old=1;  mu_error=1; mu = 1;
j=1;
  while mu_error > 1e-10
  psi=psi.*exp(-0.5*dt*(V+(g1d/hbar)*abs(psi).^2 ));
  psi_k=fft(psi);
  psi_k=psi_k.*exp(-0.5*dt*((hbar/m)*k.^2));
  psi=ifft(psi_k);
  psi=psi.*exp(-0.5*dt*(V+(g1d/hbar)*abs(psi).^2));
  psi=psi*sqrt(N)/sqrt((dx*norm(psi).^2));

  psi_redu = 0;
  if n > 1
    for i = 1:n-1
      a(i) = sum(conj(EigStat(:,i)).*psi')*dx;
      psi_redu = psi_redu + a(i)*EigStat(:,i)';
    end
    psi = psi-psi_redu;
    psi = psi/sqrt((dx*norm(psi).^2));
  else
    psi;
  end
  
  mu=log(abs(sum(conj(psi_old).*psi)*dx)^2)/dt; mu_error=abs(abs(mu-mu_old)/mu);
    if mod(j,5000) == 0
      mu_error
    end
    if j > 1e8
      'no solution found'
    break
    end
    psi_old=psi; mu_old=mu; j=j+1;
  end
end