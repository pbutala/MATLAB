% scrExtOSM_All
close all;
clearvars;
clc;
% 1. comment out close calls from all files. Keep clear vars uncommented.
% 2. remove figure call from SM and GSM files.
% 3. start iLgd at different indexes
% 4. compute clLgd array here and do not clear it
% 5. comment out clLgd creation in all files
% 6. comment out assining clLgd in all files. Do it here.

clLgd = cell(3,1);
scrExtOSM_SSK;
scrExtOSM_SM;
scrExtOSM_GSM;
legend(clLgd,'Location','SouthWest');

% FLAGS
fSTATION = 3;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fCLOSEFIGS = false;
fARCHIVE = true;

CHAROVERWRITE = '~';
if(fARCHIVE)
    CHARIDXARCHIVE = '1';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end

% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Results\5. Ext OSM All\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrExtOSM_All.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Results/5. Ext OSM All/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrExtOSM_All.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Results\\5. Ext OSM All\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrExtOSM_All.m';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes 'scrExtOSM_All' CHARIDXARCHIVE '.m'];
ctFileVars = [ctDirRes 'ExtOSM_All' CHARIDXARCHIVE '.mat'];      % file to store workspace

% delete old files
if fSAVEALL
    delete([ctDirRes '*' CHAROVERWRITE '.*']);
end

% save data
if fSAVEALL
    copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
    save(ctFileVars);                       % save workspace
end

if fSAVEALL
    fname = [ctDirRes ' OSM All BER vs SNR' CHARIDXARCHIVE];
    %     fname = [ctDirRes ' OSM BER vs SNR']; % Use this to archive
    saveas(gcf,[fname '.png'],'png');
    saveas(gcf,[fname '.fig'],'fig');
end

if fCLOSEFIGS
    close;
end




































