function DEL = getDeltaSNR(BERTH,BERLAST,RATIOS,DSNR)
% persistent pDSNR;
% if isempty(pDSNR)
%     pDSNR = DMAX;
% end
% ALPHA = -10;
% BETA = -0.5;
% if ~isequal(DBER,-inf)
%     DSNR = DMAX/(1+ALPHA*DBER);
% else
%     DSNR = BETA*abs(pDSNR);
% end
% if DSNR >= 0
%     if exist('DFMIN','var')
%         DSNR = max(abs(DFMIN),DSNR);
%     end
% else
%     if exist('DRMIN','var')
%         DSNR = min(-1*abs(DRMIN),DSNR);
%     end
% end
% pDSNR = DSNR;
RATIOS = sort(RATIOS(:));
DSNR = sort(DSNR(:));
BETA = BERLAST/BERTH;
LENR = numel(RATIOS);
DEL = DSNR(end);
for i=1:LENR
    if BETA<RATIOS(i)
        DEL = DSNR(i);
        break;
    end
end
end