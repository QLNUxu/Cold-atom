
%%part of the croation code

function [rspdmij]=rspdmij(psi,Nparticles,dx,ii,jj)

rspdmij=0;

if (jj<ii)
    return
end
    
if (ii==jj)
    P=eye(Nparticles);
else
    psi_pom0=psi(1:Nparticles,ii:jj);
    psi_pom1=psi_pom0;
%    psi_pom(1:Nparticles,ii:jj)=-psi_pom(1:Nparticles,ii:jj);
    psi_pom1(1:Nparticles,1)=0.5*psi_pom1(1:Nparticles,1);
    psi_pom1(1:Nparticles,jj-ii)=0.5*psi_pom1(1:Nparticles,jj-ii);
    P=eye(Nparticles)-2*conj(psi_pom0)*(psi_pom1).'*dx;
end
    
A=(inv(P).')*det(P);

rspdmij=(psi(:,ii)')*A*psi(:,jj);
%rspdmij=(psi(:,jj).')*A*conj(psi(:,ii));
