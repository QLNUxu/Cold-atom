function [ psi_out ] = Dyn_STA( psi, UV, UK, n, nx )

% global nx

psi_out=zeros(nx,n);

for j=1:n
    
    psi_out( : , j )=FFT( psi( :, j )', UV, UK )';
    
end

% psi_out_STA=FFT( psi_STA', UV_STA, UK )';
