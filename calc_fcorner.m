
function f_c = calc_fcorner(forces,mean_z)
% Runs a simple calculation to estimate the corner frequency in MT
% Hard-coded are the viscosity and bead-radius
% Input are the forces (vector in pN) and mean_z (vector in um)
% Output is the theoretical corner frequency f_c (vector in Hz)

eta      = 0.001;      %%% Viscosity of water, in Pa * s = N * s / m^2
%R_bead   = 0.5*10^(-6);    %%% Bead radius in m - MYONE
R_bead   = 1.4*10^(-6);    %%% Bead radius in m - M270

Length = mean_z*10^-6; %%% Tether extension in m
forces = forces*10^-12; %%% forces in N



f_c = 1/(2*pi).* forces ./ Length ./(6*pi*eta*R_bead);

end