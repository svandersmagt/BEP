function [MLfit, MLforce, fcorner] = analyze_PSD(fs,mean_z,x)
% Analyzes the Power Spectral Density of a single trace
% and runs fits
%
% Input:
% - fs: measurement frequency in Hz
% - mean_z: mean extension in um
% - x: data trace in um
% - plotflag: to plot set 1 - to not plot set 0
%
% Output:
% - MLfit: (1) alpha in pN s/nm (2) kappa in pN/nm (3) good/bad (1/0)
% - MLforce: force according to ML PSD fit in pN
% - fcorner: corner frequency (according to Lorentzian fit) in Hz

%%% Get an estimate of the force and the corner frequency
kT = 4/1000; %pN um
F_est = kT*mean_z/std(x-mean(x))^2; %pN
f_c = calc_fcorner(F_est,mean_z); %Hz
fitgood = f_c < fs/2;

%%% Add a line to x to make the length 2^integer, change to nm
%%% Subtract the mean (offset)
x(end+1) = x(1);
x = x*1000; %nm
%x = x - mean(x);
N = length(x);
logN2 = log(N)/log(2);

%%% Calc PSD, find data points below 1/20 of f_c (the corner frequency)
%%% Also throw away first point, which is merely mean(x)
%%% NOTE: 1/20 is hard-coded, seems to work fine for all data
[f, PSD, ~] = calc_powersp(x,fs);
f(1) = []; PSD(1) = [];
goodinds = f > f_c(1)/20;

%%% Max likelihood fit to PSD (using goodinds)
%%% According to the model by Lansdorp and Saleh (RSI 2012)
alpha_0 = 1E-5; kappa_0 = 4E-4;
[par, res, ef] = fminsearch(@(par) costfunction_PSDfit(f(goodinds), par(1), par(2),fs, PSD(goodinds)), [alpha_0 kappa_0]);
alpha = par(1); kappa = abs(par(2)); MLfit = [alpha kappa fitgood];
PSDmodel = analytical_PSD_overdamped_bead(alpha,kappa,fs,f);
MLforce = kappa*1000*mean_z;

%%% Fit a Lorentzian to find the corner frequency (using the 'goodinds')
fcorner_0 =1; Amp_0 = 10^(-2);
[par, res, ef] = fminsearch(@(par) norm(PSD(goodinds) - (par(1)*((1+(f(goodinds)/par(2)).^(2)).^(-1)))), [Amp_0 fcorner_0]);
fcorner = par(2); Amp  = par(1);
fcorner = abs(fcorner);
Lorentzianfit =  (Amp*(1+(f/fcorner).^(2)).^(-1));
[i ii] = min(abs(f-fcorner));

%PSD2modelTest = analytical_PSD2_overdamped_bead(MLforce,fs,f,mean_z,1400);

figure;
loglog(f(goodinds),PSD(goodinds),'r-');
hold on
loglog(f(goodinds),PSDmodel(goodinds),'b-');
loglog(f(goodinds),Lorentzianfit(goodinds),'g-');
%loglog(freq(goodinds),PSD2modelTest(goodinds),'y-');

title('Fitting in y direction using Zhongbo Yus method');
xlabel('freq (Hz)');
ylabel('Power Spectrum (nm^2/Hz)')
legend('Power spectrum', 'PSDfit', 'Lorentzianfit');
hold off
%pause;

%%% Plot if required
if 0;
    figure(1); clf; hold on; box on;
    plot(f,PSD,'b','linewidth',1)
    plot(f, PSDmodel, 'k--', 'linewidth', 2)
    %plot(f, Lorentzianfit, 'g--', 'linewidth', 2)
    plot(fcorner,Lorentzianfit(ii),'x','markersize',14,'color',[0 0.8 0],'linewidth',5)
    %line([f_c f_c],[min(PSD) max(PSD)],'color','r','linewidth',2,'linestyle','--')
    plot(f(~goodinds),PSD(~goodinds),'r','linewidth',1.5)
    
    legend('PSD','ML fit','Theoretical F_{corner}','Non-fitted data','location','southwest')
    legend(gca,'boxoff')
    
    set(gca, 'fontsize', 14, 'linewidth', 1, 'fontweight', 'bold', 'TickLength',[0.02 0.02]);
    xlabel('Frequency (Hz)')
    ylabel('PSD (nm^2/Hz)')
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    set(gca,'xlim',[f(1)/2 f(end)*2])
    title(['F_{var} = ' num2str(F_est,2) ' pN,   F_{PSD} = ' num2str(MLforce,2) ' pN'])
    
    if 0
        %%% Calculate the 'blocked' PSD, divide the data into blocks
        %%% of m points. Try to find a decent automated value of m
        %%% NOTE: mfac is thus hard-coded
        mfac = 3/4;
        m = 2^(round(mfac*logN2));
        b = N/m-1; %number of bins
        for i = 1:b
            xbin = x((i-1)*m/2+1:m+(i-1)*m/2);
            %%% 'Window' the data using a Hann window
            if 1
                Hannwindow = hann(m);
                xbin = xbin.*Hannwindow*sqrt(8/3);
            end
            [fbin PSDbin(:,i) c] = calc_powersp(xbin,fs);
        end
        %%% Average to find PSD_blocked, throw away first point = mean(xbin)
        PSD_blocked = mean(PSDbin');
        fbin(1) = []; PSD_blocked(1) = [];
        plot(fbin, PSD_blocked, '--', 'linewidth', 2,'color',[0.5 0 0.5])
        plot(f, PSDmodel, 'k--', 'linewidth', 2)
        plot(fcorner,Lorentzianfit(ii),'x','markersize',14,'color',[0 0.8 0],'linewidth',5)
    end
    %print(gcf,'-r300','-dpng','testPSD.png');
end



end
