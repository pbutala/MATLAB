close all;
clearvars;
clc;

%% 2x2
% H = [0.4 0; 0.1 0.4]
% K = [1 0; 0 1]
% PrxM = H*K*H'
% Prx = sum(diag(PrxM))

%% 1x2
% H = [0.4; 0.1]
% K = 2
% PrxM = H*K*H'
% Prx = sum(diag(PrxM))
a = sqrt(2);
CSTL = (1/sqrt(2))*[1+1i;-1+1i;1-1i;-1-1i];
M = numel(CSTL); bpsym = log2(M);
SIGPTX = (CSTL'*CSTL)/M;

SNRdb = 0:5:25;                 SNRCNT = numel(SNRdb);
SNR = power(10,SNRdb/10);

% h= 0.1;

% theoretical BER
Qf = qfunc(sqrt(SNR));
BERth = Qf;

% MC BER
SYMCNT = ceil((1e5)/bpsym);
TOTALBITS = bpsym*SYMCNT;

BERmc = zeros(1,SNRCNT);

for iSNR = 1:SNRCNT
    snr = SNR(iSNR);
    snrrt = sqrt(snr);
    err = 0;
    for iSYM = 1:SYMCNT
        bits = randi([0 1],[1,bpsym]);
        symId = bin2decMat(bits) + 1;
        
        x = CSTL(symId);
        w = (a/(2*snrrt))*randn(1) + 1i*(a/(2*snrrt))*randn(1);
        y = x + w;
        
        % nearest neighbor
        vyx = repmat(y,size(CSTL)) - CSTL;
        dyx = vyx.*vyx.'';
        symIdh = find(dyx == min(dyx),1,'first');
        bitsh = dec2binMat(symIdh-1,bpsym);
        
        err = err + biterr2(bits,bitsh);
        
    end
    BERmc(iSNR) = err/TOTALBITS;
end

% plot
figure(1);
% set(gca,'YScale','log');
semilogy(SNRdb,BERth,'b-');
hold on;
semilogy(SNRdb,BERmc,'r:');
grid on;
axis([SNRdb(1) SNRdb(end) 2*power(10,-3) 0.5]);
