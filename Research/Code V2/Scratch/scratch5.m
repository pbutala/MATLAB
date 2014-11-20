% Raised cosine filtering
close all;
clearvars;
clc;
rng('Default');

N = 128;
CLKIN = 0.5e9;
CLKOUT = 1e9;
USF = CLKOUT/CLKIN;
ns = 1:N;
uN = 1:1/USF:(N+1-1/USF);
sig = randi([0 1], N, 1);
% sig(sig==0) = -1;

uSig = updnClock(sig, CLKIN, CLKOUT, 'raisedcosine', true);



figure;
stem(ns,sig);
hold on;
plot(uN,uSig,'r');


