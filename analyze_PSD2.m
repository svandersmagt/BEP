function [MLfit, Fmag, Rfit] = analyze_PSD2(fs,Ext,x,y,R,kT,eta)
%%% Get an estimate of the force and the corner frequency
F_est = kT*Ext/std(y-mean(y))^2; %pN
f_c = calc_fcorner(F_est,Ext,R,eta); %Hz
fitgood = f_c < fs/2;
L = Ext; %nm
R_est = R; %nm

%%% Add a line to x to make the length 2^integer, change to nm
%%% Subtract the mean (offset)
x(end+1) = x(1);

%%% Calc PSD, find data points below 1/20 of f_c (the corner frequency)
%%% Also throw away first point, which is merely mean(x)
%%% NOTE: 1/20 is hard-coded, seems to work fine for all data
[f, PSD, ~] = calc_powerspblock(x,fs);
f(1) = []; PSD(1) = [];
goodinds = f > f_c(1)/20;

%%% Max likelihood fit to PSD (using goodinds)
%%% According to the model by Lansdorp and Saleh (RSI 2012)
[par, ~, ~] = fminsearch(@(par) costfunction_PSD2fit2(par(1),fs,f(goodinds),L,par(2),PSD(goodinds),kT,eta), [F_est R_est]); %(Fmag,fs,f,L,R,PSD)
Fmag = par(1); Rfit = par(2); MLfit = [Fmag Rfit fitgood];
PSDmodel = analytical_PSD2_overdamped_bead(Fmag,fs,f,L,Rfit,kT,eta);
PSD2modelTest = analytical_PSD2_overdamped_bead(F_est,fs,f,L,R_est,kT,eta);

% [par, ~, ~] = fmincon(@(par) costfunction_PSD2fit(par(1),fs,f(goodinds),L,par(2),PSD(goodinds)), [F_est R_est],[] ,[] ,[] ,[] , [F_est-20 R_est-500],[F_est+20 R_est+500]);
% Fmag = par(1); Rfit = par(2); MLfit = [Fmag Rfit fitgood];
% PSDmodel = analytical_PSD2_overdamped_bead(Fmag,fs,f,L,Rfit);

%%% Fit a Lorentzian to find the corner frequency (using the 'goodinds')
% fcornerMin_0 =1; fcornerPlus_0 =1; Amp1_0 = 10^(-3); Amp2_0 = 10^(-3);
% [par, res, ef] = fminsearch(@(par) norm(PSD(goodinds) - (par(1)*((1+(f(goodinds)/par(3)).^(2)).^(-1)))...
%     + par(2)*((1+(f(goodinds)/par(4)).^(2)).^(-1))), [Amp1_0 Amp2_0 fcornerMin_0 fcornerPlus_0]);
% fcornerMin = par(3); Amp1  = par(1);
% fcornerPlus = par(4); Amp2 = par(2);
% fcornerMin = abs(fcornerMin);
% fcornerPlus = abs(fcornerPlus);
% Lorentzianfit =  (Amp1*(1+(f/fcornerMin).^(2)).^(-1) + Amp2*(1+(f/fcornerPlus).^(2)).^(-1));

figure;
loglog(f(goodinds),PSD(goodinds),'r-');
hold on
loglog(f(goodinds),PSDmodel(goodinds),'b-');
loglog(f(goodinds),PSD2modelTest(goodinds),'y-');
% plot(freq(goodinds),Lorentzianfit(goodinds),'g-');

title('Fitting in x direction using Daldrop method');
xlabel('freq (Hz)');
ylabel('Power Spectrum (nm^2/Hz)')
legend('Power spectrum', 'PSD2fit','PSD2model with estimated values');
hold off
%pause;


end