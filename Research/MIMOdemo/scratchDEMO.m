close all;
clearvars;
clc;
rng('Default');

% CONSTANT
MOD_OFDM = 1;
MOD_OOK = 2;
MODULATION = MOD_OOK; % SELECT CORRECT MODULATION SCHEME
% MODULATION = MOD_OFDM; % SELECT CORRECT MODULATION SCHEME

FIGTITLE = 'Off';
% SYSTEM
dFs = 25e6;
dFp = 25e6;

switch(MODULATION)
    case MOD_OFDM
        fprintf('--OFDM--\n');
        demo = cDemoOFDM(dFs, dFp);
        spFrm = 1;
        BPFrm = spFrm*demo.mod.BPSYM;
        txBits = randi([0 1],demo.mod.BPSYM,1);
        
    case MOD_OOK
        %*************************** OOK *********************************%
        fprintf('--OOK--\n');
        spFrm = 4;
        BPFrm = spFrm;
        demo = cDemoOOK(dFs, dFp, spFrm);
        txBits = randi([0 1],BPFrm,1);
%         txBits = ones(BPFrm,1);
    otherwise
        error('Modulation not defined');
end

% TRANSMIT
demo.mod.write(txBits);
txSig = demo.mod.read(demo.mod.COUNTOUT);
txSig(txSig>demo.DAC.dSIGMAX) = demo.DAC.dSIGMAX;
txSig(txSig<demo.DAC.dSIGMIN) = demo.DAC.dSIGMIN;

txFrm = [demo.DAC.setRail2Rail(demo.plt.PILOT); txSig];
% txFrm = demo.DAC.setRail2Rail(plt.PILOT);
% Plot transmit frame
FIGTX = figure('Name', sprintf('Simulate Transmit'), 'NumberTitle', FIGTITLE);
plot(1:demo.frmTxNSmp16, txFrm);
axis([1 demo.frmTxNSmp16  min(txFrm) max(txFrm)]);
xlabel('Normalized time');
ylabel('Signal value');
title('Transmit frame');

% CHANNEL
h = 1;
% for SHIFT = 0:NFRM-1
for SHIFT = 0:9
    %% TRANSMIT
    % Transmit frame with DAC
    txFrmCh = demo.DAC.getOutput(txFrm);
    NFRM = demo.frmTxNSmp16;
    
    % Shift frame start (sampling offset at receiver)
    sftIs = NFRM-SHIFT+1:2*NFRM-SHIFT;
    while ~isempty(sftIs(sftIs<=0))
        sftIs(sftIs<=0) = sftIs(sftIs<=0) + NFRM;
    end
    while ~isempty(sftIs(sftIs>NFRM))
        sftIs(sftIs>NFRM) = sftIs(sftIs>NFRM) - NFRM;
    end
    txFrmCh = txFrmCh(sftIs);
    
    % Scale frame (channel gain)
    frmCh = h*txFrmCh;
    
    %% RECEIVE
    % Create figure to plot received signal and processing chain
    FIGRX = figure('Name', 'Simulate Receive', 'NumberTitle', FIGTITLE);
    FIGNR = 2;
    FIGNC = 2;
    
    % Sample transmitted frame at receiver with ADC
    rxFrmCh = updnClock(frmCh,demo.DAC.dCLKs,demo.ADC.dCLKs);
    rxFrm = demo.ADC.getOutput(rxFrmCh);
    
    % Plot receive Frame (@ADC clock)
    figure(FIGRX);
    subplot(FIGNR,FIGNC,1);
    plot(1:demo.frmRxNSmp16, rxFrm);
    axis([1 demo.frmRxNSmp16 min(rxFrm) max(rxFrm)]);
    xlabel('Normalized time');
    ylabel('Signal value');
    title('Receive frame (@ADC clock)');
        
    % Find frame starting index by aligning pilot
    idx = demo.plt.alignPilot(rxFrm-min(rxFrm), demo.ADC.dCLKs);
    
    % Find pilot in received frame
    NFRM = demo.frmRxNSmp16;
    pltIs = idx:idx+demo.plt.getLength(demo.ADC.dCLKs)-1;
    while ~isempty(pltIs(pltIs<=0))
        pltIs(pltIs<=0) = pltIs(pltIs<=0) + NFRM;
    end
    while ~isempty(pltIs(pltIs>NFRM))
        pltIs(pltIs>NFRM) = pltIs(pltIs>NFRM) - NFRM;
    end
    rxPltUS = rxFrm(pltIs);
    
    % Show pilot start and stop cursors on figure
    figure(FIGRX);
    subplot(FIGNR,FIGNC,1);
    drCursor([pltIs(1), pltIs(end)],'Vertical',['g-';'r-']);
        
    % Estimate channel gain from pilot
    hhat = demo.plt.getScale(rxPltUS, demo.ADC.dCLKs);
    scl = 1/(hhat*demo.mod.SCALE*demo.DAC.dGAIN);
    
    % Find signal in up-sampled received frame
    if pltIs(end) < pltIs(1) % if index roll over
        sigIs = pltIs(end)+1:pltIs(1)-1;
    else
        sigIs = [pltIs(end)+1:NFRM, 1:pltIs(1)-1];
    end
    rxSig = rxFrm(sigIs);
    
    % Show signal start and stop cursors on figure
    figure(FIGRX);
    subplot(FIGNR,FIGNC,1);
    drCursor([sigIs(1), sigIs(end)],'Vertical',['g:';'r:']);
    
    % Plot receive signal (@fSig*2 clock)
    figure(FIGRX);
    subplot(FIGNR,FIGNC,3);
    plot(1:demo.SPFRM*demo.demod.NPSYM, rxSig);
    axis([1 demo.SPFRM*demo.demod.NPSYM min(rxSig) max(rxSig)]);
    xlabel('Normalized time');
    ylabel('Signal value');
    title('Recovered signal');
        
    % Scale received signal
    rxSig = rxSig*scl;
    
    % Demodulate received signal
    rxSig = rxSig - min(rxSig);
    demo.demod.write(rxSig);
    rxBits = demo.demod.read(demo.demod.COUNTOUT);
    
    % Find bit errors
    berr = biterr2(txBits,rxBits);
    
    % Find align error (for known shifts) and print
    % TODO: Block enabled only for simulated hardware
    alnErr = idx-1-SHIFT;
    alnErr(alnErr<=-NFRM/2) = alnErr(alnErr<=-NFRM/2) + NFRM;
    fprintf('BERRs = %3d\n',berr);
end
























































