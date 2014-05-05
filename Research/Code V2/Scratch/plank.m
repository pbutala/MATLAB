% File: http://web.mit.edu/8.13/matlab/Examples/planck.m
% Date: 2008-June-13
% Author: Scott Sewell


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Planck Radiation Law
% See Hecht 'Optics: 4th Edition, Page 584-586' for Details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 6.626e-34;   	% Planck's Constant = 4.135 x 10^-15 eV s
c = 3e8;  		% speed of light (MKS)
T= 6000; 		% kelvin
k= 1.38066e-23; 	% Boltzmann constant in J/K
l=0:20e-9:2000e-9;
p=2*3.14*h*c*c./(l.^5);	% 
b6000=p./(exp(h*c./(l*k*T)-1));
b5000=p./(exp(h*c./(l*k*5000)-1));
b4000=p./(exp(h*c./(l*k*4000)-1));
plot(l,b6000,'.');
title('Planck Radiation Law');
xlabel('Wavelength [m]')
ylabel('Spectral Irradiance [W m^{-2} sr^{-1} nmm^{-1}]');
hold on;
plot(l,b5000,'--');
plot(l,b4000,'+');
legend('6000K','5000K','4000K');