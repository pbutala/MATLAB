% scrMLTest
close all;
clearvars;
clc;

h = modem.pammod('M', 4);              % Modulator object
g = modem.pamdemod('M', 4);            % Demodulator object

% msg = randi([0 3],4,1);               % Modulating message
msg = [0 1 2 3]'; 

modSignal = modulate(h,msg);           % Modulate signal
w = 2*rand(4,1);

txSignal = modSignal + w;      % add noise
demodSignal = demodulate(g,txSignal); % Demodulate signal

errSym = numel(find(msg ~= demodSignal));
