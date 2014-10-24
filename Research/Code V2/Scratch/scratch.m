close all;
clearvars;
clc;
rng('default');

snrdb = 20;
N = 32;
NC = 2;
n = 1:N*NC;

A = 1;

sig = A*sin(2*pi*n/N);

