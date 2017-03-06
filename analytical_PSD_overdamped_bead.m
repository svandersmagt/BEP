function PSDmodel = analytical_PSD_overdamped_bead(alpha,kappa,fs,f)
% Calculates the theoretical Power Spectral Density for Brownian motion
% in a harmonic well (Lansdorp and Saleh, RSI 2012)
% Input: alpha in pN s/nm (number)
%        kappa in pN/nm (number)
%        fs in s (number)
%        f in s (vector)
% Output: PSDmodel in nm^2/Hz (vector)

kT = 4; %pN nm
prefac = 2*kT*alpha/kappa^2;

num = 2*alpha/kappa*fs*(sin(pi*f/fs).^2)*sinh(kappa/(alpha*fs));
denom = cos(2*pi*f/fs)-cosh(kappa/(alpha*fs));

PSDmodel = prefac*(1+num./denom);

end