close all;
clearvars;
clc;

ofdmType = 'dcoofdm';
N = power(2,6);
M = 8;
ofstSDsclDco = 3.2;
ofstSDsclAco = 0; 
RES = 1e-3;
frmRes = '%0.4f';
PERS = power(2,10);
% PERS = ceil(1/RES);
MAXMC = 1e9;

[iMN, arMN, MN] = getOFDMMean(ofdmType, N, M, ofstSDsclDco, ofstSDsclAco, RES, PERS, MAXMC,gca);
fprintf(['MN = ' frmRes '\n'],MN);
fprintf('iMN = %d\n',iMN);
fprintf(['max(arMN) = ' frmRes '\n'],max(arMN));
fprintf(['min(arMN) = ' frmRes '\n'],min(arMN));