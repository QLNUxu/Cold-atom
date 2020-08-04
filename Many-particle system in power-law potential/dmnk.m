
%%%part of crpations program for the rspdm of the tonks gas

function [rspdm]=dmnk(psi,dx,dkx,Nparticles,Npts)%,nk,Angps

rspdm=zeros(Npts);
for i=1:Npts
    [rspdm(i,i)]=rspdmij(psi,Nparticles,dx,i,i);
    for j=i+1:Npts
        [rspdm(i,j)]=rspdmij(psi,Nparticles,dx,i,j);
        rspdm(j,i)=conj(rspdm(i,j));
    end
end
% nk=Npts*fftshift(diag(fft((ifft(rspdm)).'))*dx^2/(2*pi));
% nk=nk.';
% 
% uF=zeros(size(psi));
% for j=1:Nparticles
%     uF(j,:)=fftshift(abs(fft(psi(j,:))).^2)/Npts;
% end
% Angps=ones(1,Nparticles)*uF*dx/dkx;
