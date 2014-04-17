function [oNt oM obpS] = bpS2NtM(typOSM,bpS,vNt,vM,fEXACT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('fEXACT','var')
    fEXACT = true;
end
NTYPES = 8;
% order - SSK, nSM, gSM, eSM
% 1:4 - Ideal
% 5:8 - PAM
idIDSSK = 1;idIDnSM = 2;idIDgSM = 3;idIDeSM = 4;
idPAMSSK = 5;idPAMnSM = 6;idPAMgSM = 7;idPAMeSM = 8;

fPLOT = (nargout == 0);

nNt = length(vNt);
nM = length(vM);
nTyp = length(typOSM);
mbpS = zeros(nM,nNt,NTYPES);
oNt = cell(1,nTyp);
oM = cell(1,nTyp);
obpS = cell(1,nTyp);
% clLgd = cell(lenM*4+1,1);

for idM = 1:nM
    M = vM(idM);
    %%%% SSK %%%%
    mbpS(idM,:,idIDSSK) = floor(log2(vNt));
    mbpS(idM,:,idPAMSSK) = floor(log2(vNt));
    %%%% nSM %%%%
    mbpS(idM,:,idIDnSM) = floor(log2(vNt*M));
    mbpS(idM,:,idPAMnSM) = floor(log2(1+vNt*(M-1)));
    %%%% gSM %%%%
    mbpS(idM,:,idIDgSM) = floor(vNt + log2(M));
    mbpS(idM,:,idPAMgSM) = floor(log2(1+(power(2,vNt)-1)*(M-1)));
    %%%% eSM %%%%
    seeSYM = ones(length(vNt),1);
    seeSYMp = ones(length(vNt),1);
    for iNt = 1:length(vNt)
        nt = vNt(iNt);
        for k=1:nt
            seeSYM(iNt) = seeSYM(iNt) + nchoosek(nt,k)*power(M,k);
            seeSYMp(iNt) = seeSYMp(iNt) + nchoosek(nt,k)*power(M-1,k);
        end
    end
    mbpS(idM,:,idIDeSM) = floor(log2(seeSYM));
    mbpS(idM,:,idPAMeSM) = floor(log2(seeSYMp));
end
if fPLOT
    figure;
    view(3);
%     etc vars
    clLC = {'g','k','b','r',};
    lenLC = length(clLC);
    clLS = {'--','-',':','-.'};
    lenLS = length(clLS);
    clMK = {'h','o','s','d','v','^','<','>','p'};
    lenMK = length(clMK);
    
    clLgd = {'SSK' 'nSM'       'gSM'        'eSM'...
             'SSK' 'nSM M-PAM' 'gSM M-PAM'  'eSM M-PAM'};
    iLgd = 1;
    hold all;
end
for iTyp = 1:nTyp
    TYP = typOSM(iTyp);
    if (fEXACT)
        [R,C] = find(mbpS(:,:,TYP) == bpS);
    else
        [R,C] = find(mbpS(:,:,TYP) >= bpS);
    end
    Nts = vNt(C);
    Ms = vM(R);
%     mbpSTYP = mbpS(:,:,TYP);
%     bpSs = mbpSTYP(logical(V));
    bpSs = mbpS(sub2ind(size(mbpS),R,C,TYP*ones(size(R))));
    oNt{iTyp} = Nts;
    oM{iTyp} = Ms;
    obpS{iTyp} = bpSs;
    
    if fPLOT
        iLC = rem(iLgd,lenLC)+1;
        iLS = rem(iLgd,lenLS)+1;
        iMK = rem(iLgd,lenMK)+1;
        stem3(Ms,Nts,bpSs,[clLC{iLC} clLS{iLS} clMK{iMK}]);
        iLgd = iLgd + 1;
    end
end

if fPLOT
    xlabel('M');
    ylabel('N_{t}');
    zlabel('#-bits / symbol');
    grid on;
    axis([min(vM(:))-1, max(vM(:))+1, ...
        min(vNt(:))-1, max(vNt(:))+1]);
    legend(clLgd{typOSM},'Location','NorthEast');
end

end % end function