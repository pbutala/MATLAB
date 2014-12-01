% scrPilotTest
close all;
clearvars;
clc;



%% TEST PILOTBARKER
% clkPlt = 1;
% clkDAC = 1;
% clkADC = 1;
% % create pilot object
% plt = cPilotBarker('BARKER13', clkPlt, clkDAC);
% % generate frame
% frm = plt.PILOT;
% 
% for SHIFT = 0:numel(frm)-1
% % for SHIFT = 1
%     % Shift frame
%     sigSh = [frm(end-SHIFT+1:end); frm(1:end-SHIFT)];
%     
%     % Sample frame at receiver (clkADC) (to simulate sampling)
%     sigShSmp = updnClock(sigSh, clkDAC, clkADC);
%     
%     % Up-sample received signal to transmit clock
%     sigShSmpUP = updnClock(sigShSmp, clkADC, clkDAC, plt.FILTER);
%     
%     % Get alignment index
%     idx = plt.alignPilot(sigShSmpUP - min(sigShSmpUP), clkDAC);
%     
%     % Print alignment information
%     alnErr = idx-1-SHIFT;
%     NFRM = numel(frm);
%     alnErr(alnErr<=-NFRM/2) = alnErr(alnErr<=-NFRM/2) + NFRM;
%     if(abs(alnErr) >= ceil(clkDAC/(2*clkADC)))
%         fprintf('Shift = %d, idx = %d, err = %d\n',SHIFT,idx,alnErr);
%     end
% end

%% TEST PILOTTONE
nCycl = 1;
clkPlt = 25e6;
clkDAC = 1e9;
clkADC = 125e6;
% create pilot object
plt = cPilotTone(nCycl,clkPlt,clkDAC);

% generate frame
frm = plt.PILOT;

for SHIFT = 0:numel(frm)-1
% for SHIFT = 1
    % Shift frame
    sigSh = [frm(end-SHIFT+1:end); frm(1:end-SHIFT)];
    
    % Sample frame at receiver (clkADC) (to simulate sampling)
    sigShSmp = updnClock(sigSh, clkDAC, clkADC);
    
    % Get alignment index
    idx = plt.alignPilot(sigShSmp - min(sigShSmp), clkADC);
    
    % Print alignment information
    SHIFTDS = round(SHIFT*clkADC/clkDAC);
    alnErr = idx-1-SHIFTDS;
    NFRM = numel(frm);
    fprintf('Shift = %d, idx = %d, err = %d\n',SHIFTDS,idx,alnErr);
    
end
























































