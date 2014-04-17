%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation Script for LOS Illumination and Received Power
%Modified From Ghasemlooy Book
%Program 3.1
%
% Note: This assumes the source is pointed directly downward and the
% receiver is pointed directly upward.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clearvars;
clc;
dlinelw = get(0,'DefaultLineLineWidth');
set(0,'DefaultLineLineWidth',2);
daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
daxesfontsize = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter Declarations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Room Layout
%--------------------------------------------------------------------------
room_res=5;                     % Resolution of rx plane (divisions/meter)
lx=7; ly=5;                      % Room dimension (m X m)
XT=0.5; YT=0.5;                      % position of LED.(0,0) is center
h=2;                             % Dist from source to rx plane (meter)
Lux_desired = 400;               % Desired Illumination Level

% Transmitter Properties
%--------------------------------------------------------------------------
N_LEDs      = 1;               % Number of LEDs (assume single point source)
theta       = 60;                % Semi-angle at half power
m =-log10(2)/log10(cosd(theta)); % Lambertian order of emission

% Use either total lumens OR total power PER LED. The other will be calculated.
Lum_total   = 1000;              % LED Luminous Flux (lm)
%P_total     = 20;                % LED transmitted optical power (W)

% Power Spectral Density (PSD)
LAMBDAMIN=200; LAMBDAMAX=1100; DLAMBDA=1;
lambda=LAMBDAMIN:DLAMBDA:LAMBDAMAX;

% This is a simple base PSD... Update for LEDs to be used
s1=18; m1=450; a1=1; s2=60; m2=555; a2=2.15*a1; s3=25; m3=483; a3=-0.2*a1;
Sprime = a1/(sqrt(2*pi)*s1)*exp(-(lambda-m1).^2/(2*s1^2)) + ...
         a2/(sqrt(2*pi)*s2)*exp(-(lambda-m2).^2/(2*s2^2)) + ...
         a3/(sqrt(2*pi)*s3)*exp(-(lambda-m3).^2/(2*s3^2));
Sprime=Sprime/sum(Sprime);  %Normalized PSD

% Receiver Properties
%--------------------------------------------------------------------------
Adet=pi*((1e-3)^2)/4;                      % Detector physical area of a PD (m^w)
Ts=0.5;                           % Optical filter gain
index=1.5;                      % Refractive index of a lens at a PD
FOV=90*pi/180;                  % Receiver FOV
G_Con=(index^2)/sin(FOV);       % Optical concentrator gain


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Conversions between radiant power (W) and luminous power (lm)
%--------------------------------------------------------------------------
V = func_eye_sensitivity(LAMBDAMIN, LAMBDAMAX, DLAMBDA, 1978);
if (exist('Lum_total', 'var'))
    P_total=Lum_total/(683*sum(Sprime.*V)*DLAMBDA);
elseif (exist('P_total', 'var'))
    Lum_total = P_total*683*sum(Sprime.*V)*DLAMBDA;
end
%Scale by the number of LEDs
P_total = P_total*N_LEDs;
Lum_total = Lum_total*N_LEDs;

% Calculate receiver plane divisions
%--------------------------------------------------------------------------
Nx=lx*room_res; Ny=ly*room_res; 
x=-lx/2:lx/Nx:lx/2;
y=-ly/2:ly/Ny:ly/2;
[XR,YR]=meshgrid(x,y);  % receiver plane grid

% Calculate received power
%--------------------------------------------------------------------------
D1=sqrt((XR-XT(1,1)).^2+(YR-YT(1,1)).^2+h^2);     % distance vector
cosphi_A1=h./D1;                                  % angle vector
H_0=(m+1)*Adet.*cosphi_A1.^(m+1)./(2*pi.*D1.^2);  % DC channel gain
P_rec=P_total.*H_0.*Ts.*G_Con; % received power from source to rx plane (W)
P_rec_dBm=10*log10(P_rec);     % received power in decibels

% Calculate Illumination
%--------------------------------------------------------------------------
P_optical = P_rec / (Adet*Ts.*G_Con); % Received optical power
Lux_rec   = (P_optical*683*sum(Sprime.*V)*DLAMBDA); % Surface Illum (lux)
coverage  = 100 * sum(sum(Lux_rec > Lux_desired)) / (Nx * Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Emission Pattern
%--------------------------------------------------------------------------
phi = 0:0.01:(pi/2);
cos_phi = cos(phi);
cos_phi = cos_phi.^m;
cos_phi = 100 * cos_phi / norm(cos_phi);
phi = 180*phi/pi;

figure();
plot(phi, 100*cos_phi/cos_phi(1));
hold on;
plot([0 60 60],[50 50 0],'r');
axis([0 90 0 100]);
grid on;
% title(['LED Distribution for: Semi Angle ' num2str(theta) char(176) ', Lambertian Order ' num2str(m)]);
xlabel('\phi: Angle (degrees)');
ylabel('Relative Power (%)');

% Power Spectral Distribution
%--------------------------------------------------------------------------
figure();
plot(lambda, 100*Sprime/max(Sprime));
axis([300 750 0 100]);
% title('LED SPD');
xlabel('(\lambda): Wavelength (nm)');
ylabel('Relative Power (%)');
grid;

% Eye Sensitivity Function
%--------------------------------------------------------------------------
figure();
semilogy(lambda, V);
axis([350 750 1e-4 1]);
title('Eye Sensitivity Function V(\lambda)');
xlabel('wavelength (nm)');
grid;

% Received Power
%--------------------------------------------------------------------------
figure();
meshc(x,y,P_rec_dBm);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Received power (dBm)');
tStr = sprintf('Received Power- Max:%0.0fdBm, Min:%0.0fdBm', max(P_rec_dBm(:)), min(P_rec_dBm(:)));  % PB
title(tStr);
axis([-lx/2 lx/2 -ly/2 ly/2 min(min(P_rec_dBm)) max(max(P_rec_dBm))]);
rotate3d on

% Received Illumination
%--------------------------------------------------------------------------
figure();
meshc(x,y,Lux_rec);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Surface Illumination (Lux)');
tStr = sprintf('Ilumination- Max:%0.0flx, Min:%0.0flx', max(Lux_rec(:)), min(Lux_rec(:)));  % PB
title(tStr);                                                                                % PB
axis([-lx/2 lx/2 -ly/2 ly/2 min(min(Lux_rec)) max(max(Lux_rec))]);
rotate3d on

% Desired Illumination
%--------------------------------------------------------------------------
if (coverage > 0)
    figure();
    image(x,y,100*(Lux_rec > Lux_desired));
    title(['Locations where illumination > ' num2str(Lux_desired) ', Coverage: ' num2str(coverage) '%']);
    xlabel('X (m)');
    ylabel('Y (m)');
end

set(0,'DefaultLineLineWidth',dlinelw);
set(0,'DefaultAxesFontName',daxesfontname);
set(0,'DefaultAxesFontSize',daxesfontsize);