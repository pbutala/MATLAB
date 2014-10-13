close all;
clearvars;
clc;
rng('default');

Nt = 1;
Nr = 2;

% SIGPAM = [-3;-1;1;3];
SIGPAM = [-1;1];
% SIGCORR = sum(SIGPAM.^2,1)/numel(SIGPAM);
% Kx = diag(SIGCORR*ones(Nt,1));

% CNST = [SIGPAM' zeros(size(SIGPAM')) zeros(size(SIGPAM')) zeros(size(SIGPAM'));...
%         zeros(size(SIGPAM')) SIGPAM' zeros(size(SIGPAM')) zeros(size(SIGPAM'));...
%         zeros(size(SIGPAM')) zeros(size(SIGPAM')) SIGPAM' zeros(size(SIGPAM'));...
%         zeros(size(SIGPAM')) zeros(size(SIGPAM')) zeros(size(SIGPAM')) SIGPAM'];
% CNST = [SIGPAM' zeros(size(SIGPAM'));...
%         zeros(size(SIGPAM')) SIGPAM'];
CNST = SIGPAM';

Kx = cov(CNST',1);
SYMNUM = size(CNST,2);
bpsym = log2(SYMNUM);

CNSTBITASGN = getGrayBits(SYMNUM);
BITSHAMD = getHammDistMat(CNSTBITASGN);

% H = eye(Nr);
H = [0.1;0.1];

HxCNST = H*CNST;
HxCNSTDIST = getConstDistMat(HxCNST);
HxCNSTDISTSRT = sort(HxCNSTDIST);
HxCNSTDMIN = HxCNSTDISTSRT(2,:);
SIGPRX = sum(diag(H*Kx*H'),1);

SNRdb = 0:5:25;                 
SNRCNT = numel(SNRdb);
SNR = power(10,SNRdb/10);

% theoretical BER
BERup = zeros(1,SNRCNT);

% MC BER
SYMCNT = ceil(1e4/bpsym);
BITCNT = bpsym*SYMCNT;

BERmc = zeros(1,SNRCNT);
for iSNR = 1:SNRCNT
    snr = SNR(iSNR);
    Wsd = sqrt(SIGPRX/snr);
    
    for iS = 1:SYMNUM
        pMsg = 1/SYMNUM;
        pMsErr = 1 - ((1 - qfunc(HxCNSTDMIN(iS)/(2*Wsd)))^Nr);
        BERup(iSNR) = BERup(iSNR) + pMsg*pMsErr;                    % assume all bits in error if msg err. i.e. errbits/bpsym = 1.
    end
    
    err = 0;
    for iSYM = 1:SYMCNT
        bits = randi([0 1],[1,bpsym]);
        symId = bin2decMat(bits) + 1;
        
        X = CNST(:,symId);
        W = (Wsd/sqrt(Nr))*randn(Nr,1);
        Y = H*X+W;
%         Y = H*X;
        
        vXhX = repmat(Y,1,SYMNUM)-HxCNST;
        dXhX = sum(vXhX.^2,1);
        symIdh = find(dXhX == min(dXhX),1,'first');
        bitsh = dec2binMat(symIdh-1,bpsym);
        
        if biterr2(bits,bitsh) > 0
            err = err + 1;
        end
%         err = err + biterr2(bits,bitsh);
    end
%     BERmc(iSNR) = err/BITCNT;
    BERmc(iSNR) = err/SYMCNT;
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
% ylabel('BER');
ylabel('SER');
legend(gca,'Theoretical_{up}','Simulation');