close all;
clearvars;
clc;

load('\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Matfiles\LEDPSD\CIE1931_JV_1978_2deg\Gaussian\R_703_0_1_G_564_0_1_B_429_0_1\res_0.10000.mat');
RGB.obs.showGamut();
xV = [RGB.obs.gmtX];
yV = [RGB.obs.gmtY];
lV = RGB.obs.gmtL;
lV = (lV-lV(1))/(lV(end)-lV(1));

p = patch(xV,yV,lV);