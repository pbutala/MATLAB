% scratch3
figure(figOFST);
clLgdOfst = {};
hold all;
iMK = 0;
for iM = 1:lenM
    iLC = rem(1,lenLC)+1;
    iLS = rem(1,lenLS)+1;
    iMK = rem(iMK+1,lenMK)+1;
    plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
    plot(rngOfstSD,SNRD(:,iM),plStyle);
    STRM = sprintf('M:%d',rngMdco(iM));
    clLgdOfst{end+1} = ['DCO ' STRM];
    
    iLC = rem(2,lenLC)+1;
    iMK = rem(iMK+1,lenMK)+1;
    plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
    plot(rngOfstSD,SNRA(:,iM),plStyle);
    STRM = sprintf('M:%d',rngMaco(iM));
    clLgdOfst{end+1} = ['ACO ' STRM];
end
xlabel('Offset (Std Dev)');
yStr = sprintf('SNR (BER=10^{%d})',log10(BERTH));
ylabel(yStr);
tStr = sprintf('SNR vs Offset, Target BER = 10^{%d}',log10(BERTH));
title(tStr);
grid on;
axis([rngOfstSD(1) rngOfstSD(end) rngSNRdb(1) rngSNRdb(end)]);
legend(gca,clLgdOfst,'Location','SouthWest');

fname = [ctDirRes STRPREFIX 'SNR vs Offset' CHARIDXARCHIVE];
saveas(figOFST,[fname '.png'],'png');
saveas(figOFST,[fname '.fig'],'fig');