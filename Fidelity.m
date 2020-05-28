function [ F ] = Fidelity( psi_T, psi_t, n, dx  )

% global dx

F=zeros(n,1);

for j=1:n
    
    F(j,1)=abs(sum(conj(psi_T(:,j)).*psi_t(:,j))*dx)^2;

end

