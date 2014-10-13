close all;
clearvars;
clc;
rng('default');

Nt = 2;
Nr = 1;

SIGOF = -1;
SIGON = 1;
SIGCORR = (SIGON^2+SIGOF^2)/2;
Kx = diag(SIGCORR*ones(Nt,1));

CNST = [SIGOF, SIGON, SIGOF, SIGON;...
        SIGOF, SIGOF, SIGON, SIGON];
CNSTBITASGN = [0 0;0 1;1 1;1 0];
BITSHAMD = getHammDistMat(CNSTBITASGN);

SYMNUM = size(CNST,2);
bpsym = log2(SYMNUM);

H = [1 0.5];
HxCNST = H*CNST;
HxCNSTDIST = getConstDistMat(HxCNST);

SIGPRX = H*Kx*H';

SNRdb = 0:5:25;                 
SNRCNT = numel(SNRdb);
SNR = power(10,SNRdb/10);

% theoretical BER
BERup = zeros(1,SNRCNT);

% MC BER
SYMCNT = ceil(2e5/bpsym);
BITCNT = bpsym*SYMCNT;

BERmc = zeros(1,SNRCNT);
for iSNR = 1:SNRCNT
    snr = SNR(iSNR);
    Wsd = sqrt(SIGPRX/snr);
    
    for s1 = 1:SYMNUM
        for s2 = 1:SYMNUM
            pe = qfunc(HxCNSTDIST(s1,s2)/(2*Wsd));
            BERup(iSNR) = BERup(iSNR) + BITSHAMD(s1,s2)*pe/bpsym;
        end
    end
    
    err = 0;
    for iSYM = 1:SYMCNT
        bits = randi([0 1],[1,bpsym]);
        symId = bin2decMat(bits) + 1;
        
        X = CNST(:,symId);
        W = Wsd*randn(1);
        Y = H*X+W;
%         Y = H*X;
        
        vXhX = Y-HxCNST;
        dXhX = sum(vXhX.^2,1);
        symIdh = find(dXhX == min(dXhX),1,'first');
        bitsh = dec2binMat(symIdh-1,bpsym);
        
        err = err + biterr2(bits,bitsh);
    end
    BERmc(iSNR) = err/BITCNT;
end

% plot
figure(1);
% set(gca,'YScale','log');
semilogy(SNRdb,BERup,'b-');
hold on;
semilogy(SNRdb,BERmc,'r:');
grid on;
axis([SNRdb(1) SNRdb(end) 2*power(10,-5) 0.5]);
xlabel('SNR_{e} (dB)');
ylabel('BER');
legend(gca,'Theoretical (Bound_{up})','Simulation');