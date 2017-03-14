
function f_c = calc_fcorner(forces,mean_z,R,eta)
% Runs a simple calculation to estimate the corner frequency in MT
% input
% forces in pN
% mean_z, extension in nm
% R, bead radius in nm
% eta, water viscosity in pN s/nm^2

f_c = 1/(2*pi).* forces ./ mean_z ./(6*pi*eta*R);

end