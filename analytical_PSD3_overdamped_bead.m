function PSD3model = analytical_PSD3_overdamped_bead(Fmag,fs,f,L,R)
% Calculates the theoretical Power Spectral Density for Brownian motion
% in a harmonic well (Daldrop 2015)
% Input: alphaX in pN s/nm (number)
%        alphaPhi in pN s nm (number)
%        Fmag in pN (number)
%        fs in Hz (number)
%        f in Hz (vector)
%        L in nm (number)
%        R in nm (number)
% Output: PSD2model in nm^2/Hz (vector)

kT = 4.1; %pN nm
eta = 9.2E-10; %viscosity in pN s/nm^2
Cpar = (1-9/16*(1+L/R)^(-1)+1/8*(1+L/R)^(-3)-45/256*(1+L/R)^(-4)-1/16*(1+L/R)^(-5))^(-1); %Daldrop eq(S10)
Crot = 1 + 5/16*(1+L/R)^(-3); %Daldrop eq(S12)
alphaY = 6*pi*eta*R*Cpar; %Daldrop eq(11)
kappa = Fmag/L;

fc = kappa/(2*pi*alphaY);

PSD1 = 4*kT*alphaY/(kappa)^2;
PSD2 = zeros(length(f),1);
PSD3 = zeros(length(f),1);
for n=[-1 0];
    PSD2 = PSD2 + 1./(1+(abs(f+n*fs)./fc).^2);
    PSD3 = PSD3 + ((sin(pi/fs*abs(f+n*fs))).^2)./(pi/fs*abs(f+n*fs)).^2;
end



PSD3model = PSD1*(PSD2.*PSD3); %Daldrop eq (S13)

end