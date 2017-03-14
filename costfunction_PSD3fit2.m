function costfunc = costfunction_PSD3fit2(Fmag,fs,f,L,R,PSD,kT,eta)
% Costfunction for PSD Max Likelihood fit according to
% the model by Lansdorp and Saleh (RSI 2012)

PSD3model = analytical_PSD3_overdamped_bead(Fmag,fs,f,L,R,kT,eta);

costfunc = (PSD3model - PSD).^2;
costfunc = sum(costfunc);

end