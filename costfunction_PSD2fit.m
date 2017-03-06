function costfunc = costfunction_PSD2fit(Fmag,fs,f,L,R,PSD)
% Costfunction for PSD Max Likelihood fit according to
% the model by Lansdorp and Saleh (RSI 2012)

eta = 1; % As far as I know, this is constant anyway, so why calc?
PSD2model = analytical_PSD2_overdamped_bead(Fmag,fs,f,L,R);

costfunc = eta .* (PSD./PSD2model + log(PSD2model));
costfunc = sum(costfunc);

end