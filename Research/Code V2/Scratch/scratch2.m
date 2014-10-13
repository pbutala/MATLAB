close all;
clearvars;
clc;
rng('default');

Xs = [1 0; 0 1]; 
Kx = cov(Xs',1);
H = [1 0; 0 1];
Ys = H*Xs;
Ky = cov(Ys',1);
Ky2 = H*Kx*H';
