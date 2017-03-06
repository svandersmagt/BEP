function [Ext, Fx_real, Fy_real, PSDfit, PSDforce, fcorner, MLfit, MLforce, Rfit] = analyze_one_trace2(time, x, y, z, fs, R)
%
% Function to analyze magnetic tweezers time traces
%
% Input: (time, x, y, z, filenr, F_IND, fs, plotflag)
% - time in s (vector)
% - trace x in mum (vector, same length)
% - trace y in mum (vector, same length)
% - trace z in mum (vector, same length)
% - beadnumber (integer; display in the title of the (x,y,z) plot)
% - F_IND - set 1 to analyze x fluctuations, set 2 to analyze y fluctuations
% - fs in Hz (number, sampling frequency)
% if plotflag == 1: Show the x,y,z time trace and histograms, as well as
% the fits for the PSD, AV and SA methods

foo = 1;

%%% Global, "hard-coded" variables
binWidth_XY = 0.01; % in mum
binWidth_Z  = 0.01; % in mum
kT = 4/1000;

%%% Analyze the real space fluctuations
mean_x = mean(x);
mean_y = mean(y);
mean_z = mean(z);
std_x  = std(x);
std_y  = std(y);
%std_z  = std(z);

Fx_real = kT*mean_z./std_x^2;
Fy_real = kT*mean_z./std_y^2;
Ext = mean_z;

%%% Analyze the fluctuations in freq domain using Power Spectral Density (PSD) max likelihood fit
%%% Also analyze the real space fluctuations using the Allan variance
%%% These implement the two methods of Landsdorp and Saleh (RSI 2012)
%%% The third method is that by te Velthuis et al. (Bioph J 2010)

    [PSDfit, PSDforce, fcorner] = analyze_PSD(fs,mean_z,y);
    [MLfit, MLforce, Rfit] = analyze_PSD2(fs,mean_z,x,R);


if 0;
    %%% If plotflag = 1, show the x,y,z-traces with histograms
    %%% and the x-y, x-z and y-z plots
    %%% Code from show_xyz
    warning off
    figure(5);clf;
    
    subplot(3,3,1); hold on; box on;
    title([ 'F_x = ' num2str(Fx_real, 2) ' pN'], 'fontsize', 14);
    plot(time, x)
    ylabel('X (\mum)', 'fontsize', 14);
    xlabel('Time  (s)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    subplot(3,3,2); hold on; box on;
    title([ 'F_y = ' num2str(Fy_real, 2) ' pN'], 'fontsize', 14);
    plot(time, y, 'r')
    ylabel('Y (\mum)', 'fontsize', 14);
    xlabel('Time  (s)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    
    subplot(3,3,3); hold on; box on;
    title([ '<z> = ' num2str(mean_z, 3) ' \mum'], 'fontsize', 14);
    plot(time, z, 'Color', [0 0.5 0])
    ylabel('Z (\mum)', 'fontsize', 14);
    xlabel('Time  (s)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    
    %%% Plot X vs. Y
    subplot(3,3,7); hold on; box on;
    title([ 'Bead number: ' num2str(beadnumber)], 'fontsize', 14);
    
    plot(x, y, 'k')
    
    maxval = max([max(x) abs(min(x)) max(y) abs(min(y))]);
    
    if maxval > 2.1
        axis([-3 3 -3 3])
    elseif maxval > 1.1
        axis([-2 2 -2 2])
    else
        axis([-1 1 -1 1])
    end
    axis square
    xlabel('X (\mum)', 'fontsize', 14);
    ylabel('Y (\mum)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    %%% Plot X vs. Z
    subplot(3,3,8); hold on; box on;
    plot(x, z, 'k')
    
    xlabel('X (\mum)', 'fontsize', 14);
    ylabel('Z (\mum)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    %%% Plot Y vs. Z
    subplot(3,3,9); hold on; box on;
    plot(y, z, 'k')
    
    xlabel('Y (\mum)', 'fontsize', 14);
    ylabel('Z (\mum)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    
    
    %%% Subtract out the means
    x = x - mean_x;
    y = y - mean_y;
    z = z - mean_z;
    binWidth = 0.001;
    
    
    %%% Analyze X %%%
    %%%--- Fit an exponential distribution ---%%%
    binCenters = min(x):binWidth:max(x);
    [MUHAT,SIGMAHAT,MUCI,SIGMACI] = normfit(x);
    xfit = min(x):(binWidth*0.1):max(x);
    yfit = length(x)* binWidth * normpdf(xfit, MUHAT, SIGMAHAT);
    
    
    subplot(3,3,4); hold on; box on;        %%% X histogramm
    hist(x, binCenters)
    h1 = findobj(gca,'Type','patch');
    set(h1,'FaceColor', [0 0 0],'EdgeColor','k')
    
    plot(xfit, yfit, 'color', [0 0 1], 'linewidth', 2)
    ylabel('Counts', 'fontsize', 14);
    xlabel('X (\mum)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    title(['\sigma  = ' num2str(SIGMAHAT*1000,3) ' nm ' ])
    
    %%% Analyze Y %%%
    %%%--- Fit an exponential distribution ---%%%
    binCenters = min(y):binWidth:max(y);
    [MUHAT,SIGMAHAT,MUCI,SIGMACI] = normfit(y);
    xfit = min(y):(binWidth*0.1):max(y);
    yfit = length(y)* binWidth * normpdf(xfit, MUHAT, SIGMAHAT);
    
    
    subplot(3,3,5); hold on; box on;        %%% X histogramm
    hist(y, binCenters)
    h1 = findobj(gca,'Type','patch');
    set(h1,'FaceColor', [0 0 0],'EdgeColor','k')
    
    plot(xfit, yfit, 'color', [1 0 0], 'linewidth', 2)
    ylabel('Counts', 'fontsize', 14);
    xlabel('Y (\mum)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    title(['\sigma  = ' num2str(SIGMAHAT*1000,3) ' nm ' ])
    
    
    %%% Analyze Z %%%
    %%%--- Fit an exponential distribution ---%%%
    binCenters = min(z):binWidth:max(z);
    [MUHAT,SIGMAHAT,MUCI,SIGMACI] = normfit(z);
    xfit = min(z):(binWidth*0.1):max(z);
    yfit = length(z)* binWidth * normpdf(xfit, MUHAT, SIGMAHAT);
    
    
    subplot(3,3,6); hold on; box on;        %%% X histogramm
    hist(z, binCenters)
    h1 = findobj(gca,'Type','patch');
    set(h1,'FaceColor', [0 0 0],'EdgeColor','k')
    
    plot(xfit, yfit, 'color', [0 0.5 0], 'linewidth', 2)
    ylabel('Counts', 'fontsize', 14);
    xlabel('Z (\mum)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    title(['\sigma  = ' num2str(SIGMAHAT*1000,3) ' nm ' ])
    %print(gcf,'-r300','-dpng','testtrace.png');
    warning on
end

end