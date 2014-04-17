% script to calculate specs for receiver
% 

close all;
clearvars;
clc;

%% INIT
SNRdb = 25;                 % SNR target (dB)
Brx = 20e6;                 % Receiver bandwidth  (Hz)
R = 0.1;                    % Responsivity (A/W)
Il = 400;                   % Target illumination (lx = lm/m2)
Ntx = 9;                    % Number of transmitters
Tb = 0.2;                   % Relative power in blue times transmission of blue filter
K = 4e-3;                   % Optical to electrical power conversion ratio (W/lm)
In = (5e-12)*sqrt(Brx);     % Noise current (A)

%% CALC
% Calculate collection area required to achieve target SNR
Ao = ((In*Ntx)/(R*Il*K*Tb))*power(10,(SNRdb/40));
% Calculate aperture diameter
Do = sqrt(4*Ao/pi);

%% DISP
str = sprintf('Ao = %0.2e',Ao);
disp(str);
str = sprintf('Do = %0.2e',Do);
disp(str);