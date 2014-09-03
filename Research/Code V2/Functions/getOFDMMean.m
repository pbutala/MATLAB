function [iMN, arMN, MN] = getOFDMMean(ofdmType, N, M, ofstSDsclDco, ofstSDsclAco, RES, PERS, MAXMC, Hax)
if ~exist('RES','var')
    RES = 1e-3;
end
if ~exist('PERS','var')
    PERS = power(2,10);
end
if ~exist('MAXMC','var')
    MAXMC = 1e9;
end
fPLOT = exist('Hax','var');

switch lower(ofdmType)
    case 'acoofdm'
        d = N/4;
    case 'dcoofdm'
        d = N/2-1;
end
SYMS = getQAMsyms(M);
% BPS = d*log2(M);

if fPLOT
    axes(Hax);
    hold on;
    xlabel('Iteration'); ylabel('Mean'); title(['Mean computation for ' ofdmType]);
    grid on;
end

arMN = -1*ones(PERS,1);
iLOOP = 0; iMN = PERS;
LOOPDONE = false;

while ~LOOPDONE
    iLOOP = iLOOP + 1;
    oMN = iMN;
    iMN = rem(iMN,PERS) + 1;
    
    bits = randi([0 1],[d,log2(M)]);
    symIdx = bin2decMat(bits)+1;
    
    tSig = genOFDMsignal(...
        'data',symIdx,...
        'OFDMtype',ofdmType,...
        'N',N,...
        'Symbols',SYMS,...
        'OffsetDcoStddev', ofstSDsclDco,...
        'OffsetAcoStddev', ofstSDsclAco);
    
    tMN = mean(tSig);
    
    arMN(iMN) = (arMN(oMN)*(iLOOP-1) + tMN)/iLOOP;
    dMN = max(arMN) - min(arMN);
    if dMN <= RES
        LOOPDONE = true;
    end
    
    if iLOOP >= MAXMC
        warning('Escape from monte-carlo loop. Could not satisfy loop terminate criterion');
        LOOPDONE = true;
    end
    
    if fPLOT
        scatter(Hax,iLOOP,arMN(iMN),'.b');
    end
end
MN = arMN(iMN);
if fPLOT
    axis([1 iLOOP 0 ceil(2*MN)]);
end
end