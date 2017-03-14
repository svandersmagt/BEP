function costfunc = costfunction_PSD2fit2(Fmag,fs,f,L,R,PSD,kT,eta)
% Costfunction for PSD Max Likelihood fit according to
% the model by Lansdorp and Saleh (RSI 2012)

PSD2model = analytical_PSD2_overdamped_bead(Fmag,fs,f,L,R,kT,eta);

costfunc = (PSD2model - PSD).^2.;
costfunc = sum(costfunc);

end