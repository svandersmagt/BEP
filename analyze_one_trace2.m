function [Ext, Fx_real, Fy_real, PSDfit, PSDforce, fcorner, MLfitx, MLforcex, Rfitx, MLfity, MLforcey, Rfity] = analyze_one_trace2(time, x, y, z, fs, R, kT, eta)
%
% Function to analyze magnetic tweezers time traces
%
% Input: (time, x, y, z, filenr, F_IND, fs, plotflag)
% - time in s (vector)
% - trace x in nm (vector, same length)
% - trace y in nm (vector, same length)
% - trace z in nm (vector, same length)
% - fs in Hz (number, sampling frequency)
% - bead radius R in nm
% - kT in pN nm
% - eta in pN s/nm^2 (viscosity
% the fits for the PSD, AV and SA methods

%%% Analyze the real space fluctuations
mean_x = mean(x);
mean_y = mean(y);
mean_z = mean(z);
std_x  = std(x);
std_y  = std(y);

Fx_real = kT*mean_z./std_x^2;
Fy_real = kT*mean_z./std_y^2;
%Ext = mean(sqrt((x-mean_x).^2 + (y-mean_y).^2 + z.^2)); haha
Ext = mean_z;

%%% Analyze the fluctuations in freq domain using Power Spectral Density (PSD) max likelihood fit
%%% Also analyze the real space fluctuations using the Allan variance
%%% These implement the two methods of Landsdorp and Saleh (RSI 2012)
%%% The third method is that by te Velthuis et al. (Bioph J 2010)

    [PSDfit, PSDforce, fcorner] = analyze_PSD(fs,Ext,y-mean_y,kT,R,eta);
    [MLfitx, MLforcex, Rfitx] = analyze_PSD2(fs,Ext,x-mean_x,y-mean_y,R,kT,eta);
    [MLfity, MLforcey, Rfity] = analyze_PSD3(fs,Ext,y,R,kT,eta);

end