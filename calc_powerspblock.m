function [f,P,T] = calc_powerspblock(X,sampling_f)
% function [f,P,T] = calc_powersp(X,sampling_f)
% Calculates the powerspectrum P as function of frequency f
% of imput data X sampled at frequency sampling_f

nBlock = 40; %number of blocks
lengthBlock = floor(length(X)/nBlock);

fNyq    =   sampling_f / 2;
delta_t =   1 / sampling_f;
time    =   [0 : delta_t : (lengthBlock-1)*delta_t]';
T       =   max(time);
f       =   ([1 : lengthBlock] / T)';
ind     =   find(f <= fNyq); % only to the Nyquist f
f       =   f(ind);
PSD     =   zeros(3500,1);


for i = 1:40;
    FT      =   delta_t*fft(X((lengthBlock*(i-1)+1):lengthBlock*i));
    P       =   FT .* conj(FT) / T;
    P       =   P(ind);
    PSD     =   PSD + P;
    
end
P      =   PSD/nBlock;

end
