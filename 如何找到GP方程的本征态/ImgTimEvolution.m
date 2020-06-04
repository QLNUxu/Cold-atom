function [psi,mu] = ImgTimEvolution(psi,dt,dx,g1d,x,k,m,V)
hbar=1; 
% dx=x(2)-x(1); 
% dk=2*pi/(x(end)-x(1));
N=dx*norm(psi).^2; 
Nx=length(x); 
psi_mid_old=psi((Nx)/2);
mu_old=1; j=1; mu_error=1;
  while mu_error > 1e-8
  psi=psi.*exp(-0.5*dt*(V+(g1d/hbar)*abs(psi).^2));
  psi_k=fftshift(fft(psi))/Nx;
  psi_k=psi_k.*exp(-0.5*dt*(hbar/m)*k.^2);
  psi=ifft(ifftshift(psi_k))*Nx;
  psi=psi.*exp(-0.5*dt*(V+(g1d/hbar)*abs(psi).^2));
  psi_mid=psi((Nx)/2);
  mu=log(psi_mid_old/psi_mid)/dt; mu_error=abs(mu-mu_old)/mu;
  psi=psi*sqrt(N)/sqrt((dx*norm(psi).^2));
    if mod(j,5000) == 0
      mu_error
    end
    if j > 1e8
      'no solution found'
    break
    end
psi_mid_old=psi((Nx)/2); mu_old=mu; j=j+1;
  end
end