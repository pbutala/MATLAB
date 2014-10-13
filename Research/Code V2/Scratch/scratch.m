close all;
clearvars;
clc;

%% 2x2
H = [0.4 0; 0.1 0.4]
K = [1 0; 0 1]
PrxM = H*K*H'
Prx = sum(diag(PrxM))

%% 1x2
% H = [0.4; 0.1]
% K = 2
% PrxM = H*K*H'
% Prx = sum(diag(PrxM))

% TOTALBITS = 1e5;
% SIGOFF = 0;
% SIGON = sqrt(2);
% SIG = [SIGOFF SIGON];
% SIGPTX = 0.5*(SIGOFF^2) + 0.5*(SIGON^2);
% 
% SNRdb = 0:2.5:25;                 SNRCNT = numel(SNRdb);
% SNR = power(10,SNRdb/10);
% 
% h= 0.1;
% 
% % theoretical BER
% BERth = qfunc(sqrt(SNR/2));
% 
% % MC BER
% BERmc = zeros(1,SNRCNT);
% for iSNR = 1:SNRCNT
%     snr = SNR(iSNR);
%     err = 0;
%     for iBITS = 1:TOTALBITS
%         bit = randi([0 1],1);
%         x = SIG(bit+1);
%         w = sqrt(((SIGON^2)*(h^2))/(2*snr))*randn(1);
%         y = h*x + w;
%         
%         if(y>h*SIGON/2)
%             bith = 1;
%         else 
%             bith = 0;
%         end
%         
%         if(bith ~= bit)
%             err = err + 1;
%         end
%     end
%     BERmc(iSNR) = err/TOTALBITS;
% end
% 
% % plot
% figure(1);
% % set(gca,'YScale','log');
% semilogy(SNRdb,BERth,'b-');
% hold on;
% semilogy(SNRdb,BERmc,'r:');
% grid on;
% axis([SNRdb(1) SNRdb(end) 2*power(10,-3) 0.5]);
