close all;
clearvars;
clc;

a = sqrt(2);
CSTL = (1/sqrt(2))*[1+1i;-1+1i;1-1i;-1-1i];
M = numel(CSTL); bpsc = log2(M);
SIGPTX = (CSTL'*CSTL)/M;

SNRdb = 0:5:25;                 SNRCNT = numel(SNRdb);
SNR = power(10,SNRdb/10);

% h= 0.1;

% theoretical BER
% Qf = qfunc(sqrt(SNR));
Qf = qfunc(sqrt(2*SNR)); % for DCO, data repeated over 2 subcarriers, thus 2x power per stream.
% Qf = qfunc(sqrt(SNR/2));
BERth = Qf;

% MC BER
N = 64;
% bpsym = bpsc*(N/2-1); % dco
bpsym = bpsc*(N/4); % aco

SYMCNT = ceil((1e4)/bpsym);
TOTALBITS = bpsym*SYMCNT;

BERmc = zeros(1,SNRCNT);

for iSNR = 1:SNRCNT
    snr = SNR(iSNR);
    snrrt = sqrt(snr);
    err = 0;
    for iSYM = 1:SYMCNT
%         bits = randi([0 1],[N/2-1,bpsc]); % dco
        bits = randi([0 1],[N/4,bpsc]); % aco
        symId = bin2decMat(bits) + 1;
        
        X = genOFDMsignal(...
            'data',symId,...
            'OFDMtype','acoofdm',...
            'N',N,...
            'Symbols',CSTL,...
            'OffsetDcoStddev', 3.2,...
            'OffsetAcoStddev', 0);                    
%         fSig = [0;CSTL(symId)];
%         X = ifft(fSig,N,1)*sqrt(N);
        
%         Y = X;
%         W = (a/(2*snrrt))*randn(N,1) + 1i*(a/(2*snrrt))*randn(N,1);
        W = (a/(2*snrrt))*randn(N,1);
        Y = X + W;
        
        symIdh = decodeOFDMsignal(Y,...
            'OFDMtype','acoofdm',...
            'N',N,...
            'Symbols',CSTL);
        bitsh = dec2binMat(symIdh-1,bpsc);
        err = err + biterr2(bits,bitsh);
%         fSigh = fft(Y,N,1)/sqrt(N);
%         for iSC = 2:N
%             % nearest neighbor
%             y = fSigh(iSC);
%             vyx = repmat(y,size(CSTL)) - CSTL;
%             dyx = vyx.*vyx.'';
%             symIdh = find(dyx == min(dyx),1,'first');
%             bitsh = dec2binMat(symIdh-1,bpsc);
%             err = err + biterr2(bits(iSC-1,:),bitsh);
%         end
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
