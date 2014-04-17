% scrSNRvsBER
close all;
clearvars;
clc;

% DEFAULT COSMETIC SETTINGS
dlinelw = get(0,'DefaultLineLineWidth');
set(0,'DefaultLineLineWidth',2);
daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
daxesfontsize = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',16);
dfigvis = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','On');

%% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fCLOSEFIGS = false;

%% VARIABLE PARAMS
rngSNRdb = 1:15;
rngSNR = power(10,rngSNRdb/10);

%% SET CONSTELLATION
k = 3; % bits per symbols
M = 2^k; % # of symbols
syms = getPAMsyms(M,1,1); %PAM modulation
if ~isempty(find(syms<0,1))
    error('Modulation symbol cannot have a negative value');
end
syms = M*syms/sum(syms(:));
Pavg = mean(syms(:));
% display(sum(syms(:)));
% display(mean(syms(:)));

%%  Set variables
% ITERBITS = 1e3;
% ITERSYMS = ceil(ITERBITS/k);
% ITERBITS = ITERSYMS*k;
ITERSYMS = 1e3;
%% SETUP
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Results\0. Other\1. BERvsSNR\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrBERvsSNR.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Results/0. Other/1. BERvsSNR/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrBERvsSNR.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Results\\0. Other\\1. BERvsSNR\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrBERvsSNR.m';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes 'scrBERvsSNR.m'];
ctFileVars = [ctDirRes 'dataBERvsSNR.mat'];

sym_err = zeros(length(rngSNRdb),1);
for idSNRdb = 1:length(rngSNRdb)
    vSNRdb = rngSNRdb(idSNRdb);
    vSNR = power(10,vSNRdb/10);
    vSNRrt = sqrt(vSNR);
    err = 0;
    for ii = 1:1:ITERSYMS
        bit_stream = randi([0,1], [1,k]);
        sym_index = bin2dec(int2str(bit_stream));
        x = syms(sym_index+1);
        w = (Pavg/vSNRrt)*randn(1,1);
        y = x + w;
        
        dysym = abs(y-syms);
        xh = syms(find(dysym == min(dysym),1));
%         xhb = dec2bin(xh-1,k);         % generate at least log2(M) 
%         for xx = 1:k
%             sym_bits(xx)= bin2dec(xhb(xx)); 
%         end
%         err = err + biterr2(bit_stream,sym_bits);
        err = err + (x ~= xh);
    end
%     bit_err(idSNRdb) = err/ITERBITS;
    sym_err(idSNRdb) = err/ITERSYMS;
end
semilogy(rngSNRdb,sym_err);                % plot bit error vs snr
grid on;
hold all;

sym_err_aly = qfunc(sqrt(rngSNR)/(M-1))*(2*(M-1)/M);
semilogy(rngSNRdb,sym_err_aly,'-.');                % plot bit error vs snr
% restore defaults
set(0,'DefaultLineLineWidth',dlinelw);
set(0,'DefaultAxesFontName',daxesfontname);
set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);






















