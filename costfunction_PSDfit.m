function costfunc = costfunction_PSDfit(f,alpha,kappa,fs,PSD,kT)
% Costfunction for PSD Max Likelihood fit according to
% the model by Lansdorp and Saleh (RSI 2012)

eta = 1; % As far as I know, this is constant anyway, so why calc?
PSDmodel = analytical_PSD_overdamped_bead(alpha,kappa,fs,f,kT);

costfunc = eta .* (PSD./PSDmodel + log(PSDmodel));
costfunc = sum(costfunc);

end