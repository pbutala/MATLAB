close all;
clearvars;
clc;

% CSTL = (1/sqrt(2))*[1+1i;-1+1i;1-1i;-1-1i];
CSTL = getQAMsyms(64,true);
M = numel(CSTL); bpsym = log2(M);
Mrt = sqrt(M);
SIGPTX = (CSTL'*CSTL)/M;
h=1;
SIGPRX = h*SIGPTX;
SIGPRXrt = sqrt(SIGPRX);

SNRdb = 0:5:25;                 SNRCNT = numel(SNRdb);
SNR = power(10,SNRdb/10);

% theoretical BER
BERth = ((Mrt-1)/(Mrt*log2(Mrt)))*erfc(sqrt((3*SNR)/(2*(M-1))));

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
        wsd = SIGPRXrt/(sqrt(2)*snrrt);
        w = wsd*randn(1) + 1i*wsd*randn(1);
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
xlabel('SNR_{e} (dB)');
ylabel('BER');
legend(gca,'Theoretical','Simulation');