function costfunc = costfunction_PSD2fit(Fmag,fs,f,L,R,PSD,kT,eta)
% Costfunction for PSD Max Likelihood fit according to
% the model by Lansdorp and Saleh (RSI 2012)

eta1 = 1; % As far as I know, this is constant anyway, so why calc?
PSD2model = analytical_PSD2_overdamped_bead(Fmag,fs,f,L,R,kT,eta);

costfunc = eta1 .* (PSD./PSD2model + log(PSD2model));
costfunc = sum(costfunc);

end