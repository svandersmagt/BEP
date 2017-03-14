
X = bead(1).x;
sampling_f = 2800;

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
loglog(f,P,'r-');