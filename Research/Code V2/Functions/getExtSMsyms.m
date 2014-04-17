function SYMS = getExtSMsyms(Nt,M,Iavg)
nSYMS = M^Nt;
MM = dec2binMat(0:nSYMS-1)';

fSYMS = zeros(Nt,nSYMS);
m = log2(M);

for iNt = 1:Nt
    nbeg = (iNt-1)*m + 1;
    nend = iNt*m;
    fSYMS(iNt,:) = bin2decMat(MM(nbeg:nend,:)');
end
% Xs = (1:M)*2*Iavg/(M+1);
% MM = zeros(Nt,M^Nt);
% for R = 1:Nt
%     MM(R,:) = reshape(repmat(0:M-1,M^(R-1),M^(Nt-R)),1,M^Nt);
% end
% tsym = repmat(MM,1,2^Nt);
% ftx = reshape(repmat(dec2binMat(0:(power(2,Nt)-1))',M^Nt,1),Nt,(M^Nt)*(2^Nt));
% tsym = Xs(tsym.*ftx + 1);
% SYMS = removeDuplicates(tsym);

Xs = (1:M)*2*Iavg/(M+1);
SYMS = Xs(1+fSYMS);
end

% function M = removeDuplicates(V)
% [~,nC] = size(V);
% fOK = ones(1,nC);
% for iC0 = 1:nC
%     C0 = V(:,iC0);
%     for iC1 = iC0+1:nC
%         if fOK(iC1)
%             C1 = V(:,iC1);
%             fOK(iC1) = ~isequal(C0,C1);
%         end
%     end
% end
% M = V(:,logical(fOK));
% end