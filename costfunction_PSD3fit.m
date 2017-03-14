function costfunc = costfunction_PSD3fit(Fmag,fs,f,L,R,PSD,kT,eta)
% Costfunction for PSD Max Likelihood fit according to
% the model by Lansdorp and Saleh (RSI 2012)

eta1 = 1; % As far as I know, this is constant anyway, so why calc?
PSD3model = analytical_PSD3_overdamped_bead(Fmag,fs,f,L,R,kT,eta);

costfunc = eta1 .* (PSD./PSD3model + log(PSD3model));
costfunc = sum(costfunc);

end