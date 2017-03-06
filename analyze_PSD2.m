function [MLfit, Fmag, Rfit] = analyze_PSD2(fs,mean_z,y,R)
%%% Get an estimate of the force and the corner frequency
kT = 4/1000; %pN um
F_est = kT*mean_z/std(y-mean(y))^2; %pN
f_c = calc_fcorner(F_est,mean_z); %Hz
fitgood = f_c < fs/2;
L = mean_z * 1000; %nm
R_est = R*1000; %nm

%%% Add a line to x to make the length 2^integer, change to nm
%%% Subtract the mean (offset)
y(end+1) = y(1);
y = y*1000; %nm
%x = x - mean(x);
N = length(y);
logN2 = log(N)/log(2);


%%% Calc PSD, find data points below 1/20 of f_c (the corner frequency)
%%% Also throw away first point, which is merely mean(x)
%%% NOTE: 1/20 is hard-coded, seems to work fine for all data
[f, PSD, ~] = calc_powersp(y,fs);
f(1) = []; PSD(1) = [];
goodinds = f > f_c(1)/20;

%%% Max likelihood fit to PSD (using goodinds)
%%% According to the model by Lansdorp and Saleh (RSI 2012)
% [par, ~, ~] = fminsearch(@(par) costfunction_PSD2fit(par(1),fs,f(goodinds),L,par(2),PSD(goodinds)), [F_est R_est]); %(Fmag,fs,f,L,R,PSD)
% Fmag = par(1); Rfit = par(2); MLfit = [Fmag Rfit fitgood];
% PSDmodel = analytical_PSD2_overdamped_bead(Fmag,fs,f,L,Rfit);

[par, ~, ~] = fmincon(@(par) costfunction_PSD2fit(par(1),fs,f(goodinds),L,par(2),PSD(goodinds)), [F_est R_est],[] ,[] ,[] ,[] , [F_est-20 R_est-500],[F_est+20 R_est+500]);
Fmag = par(1); Rfit = par(2); MLfit = [Fmag Rfit fitgood];
PSDmodel = analytical_PSD2_overdamped_bead(Fmag,fs,f,L,Rfit);

%%% Fit a Lorentzian to find the corner frequency (using the 'goodinds')
% fcornerMin_0 =1; fcornerPlus_0 =1; Amp1_0 = 10^(-3); Amp2_0 = 10^(-3);
% [par, res, ef] = fminsearch(@(par) norm(PSD(goodinds) - (par(1)*((1+(f(goodinds)/par(3)).^(2)).^(-1)))...
%     + par(2)*((1+(f(goodinds)/par(4)).^(2)).^(-1))), [Amp1_0 Amp2_0 fcornerMin_0 fcornerPlus_0]);
% fcornerMin = par(3); Amp1  = par(1);
% fcornerPlus = par(4); Amp2 = par(2);
% fcornerMin = abs(fcornerMin);
% fcornerPlus = abs(fcornerPlus);
% Lorentzianfit =  (Amp1*(1+(f/fcornerMin).^(2)).^(-1) + Amp2*(1+(f/fcornerPlus).^(2)).^(-1));

freq = 1:length(PSD);
figure;
hold on
loglog(freq(goodinds),PSD(goodinds),'r-');
loglog(freq(goodinds),PSDmodel(goodinds),'b-');
% plot(freq(goodinds),Lorentzianfit(goodinds),'g-');
hold off
%pause;


end