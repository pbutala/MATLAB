% scrtestcodeopt
close all;
clearvars;
clc;

Nt= 4;
M = 4;
Xs = getPAMsyms(M,1,1);   % get PAM symbols
SYMBSNT = [zeros(Nt,1) reshape(repmat(dec2binMat(1:(power(2,Nt)-1)).',M-1,1),Nt,(M-1)*(power(2,Nt)-1))];
SYMBSMOD = [zeros(Nt,1) repmat(repmat(Xs(2:end).',Nt,1),1,power(2,Nt)-1)];
tSYMBS = SYMBSNT.*SYMBSMOD;

nBPSYM = floor(log2(size(tSYMBS,2)));
nSYM = power(2,nBPSYM);
optSYM = getCodeOpt(tSYMBS);
SYMBS = optSYM(:,1:nSYM);

% calculate progressively dmin
Dmin = zeros(nSYM-1,1);

for iD = 1:nSYM-1
    nCODE = iD+1;
    C = SYMBS(:,1:nCODE);
    for iC0 = 1:nCODE               % iter all codes
        C0 = C(:,iC0);              % select first
        for iC1 = iC0 : nCODE       % iter all 'next' codes
            C1 = C(:,iC1);          % select next
            V = C1-C0;              % calc vector from C0 to C1
            D = sqrt(sum(V.*V));    % calc distance from C0 to C1
            mD(iC0,iC1) = D;        % store distance in matrix
            mD(iC1,iC0) = D;        % store distance in matrix (across diag)
        end
    end
    Dmin(iD) = min(mD(mD ~= 0));
end

plot(1:nSYM-1,Dmin);