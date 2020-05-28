function [ psi_out ] = FFT( psi_in, UV, UK )

    psi_1=fft(psi_in);
    psi_2=UK.*psi_1;
    psi_3=ifft(psi_2);
    psi_4=UV.*psi_3;
    psi_5=fft(psi_4);
    psi_6=UK.*psi_5;
    psi_7=ifft(psi_6);
    psi_out=psi_7;

end

