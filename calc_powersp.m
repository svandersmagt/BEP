function [f,P,T] = calc_powersp(X,sampling_f)
% function [f,P,T] = calc_powersp(X,sampling_f)
% Calculates the powerspectrum P as function of frequency f
% of imput data X sampled at frequency sampling_f



fNyq    =   sampling_f / 2;
delta_t =   1 / sampling_f;

time    =   [0 : delta_t : (length(X)-1)*delta_t]';
T       =   max(time);
f       =   ([1 : length(X)] / T)';

FT      =   delta_t*fft(X);
P       =   FT .* conj(FT) / T;

ind     =   find(f <= fNyq); % only to the Nyquist f
f       =   f(ind);
P       =   P(ind);


end
