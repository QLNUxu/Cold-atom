function [ Prb ] = Probability( beta, eigenvalue, n )

Prb=zeros(n,1);
for j=1:n
    
    Prb(j,1)=exp(-beta*eigenvalue(j,1));
    
end

NN=sum(Prb);
Prb=Prb/NN;



