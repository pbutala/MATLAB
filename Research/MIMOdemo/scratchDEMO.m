close all;
clearvars;
clc;
rng('Default');

% CONSTANT
MOD_OFDM = 1;
MOD_OOK = 2;
% MODULATION = MOD_OOK; % SELECT CORRECT MODULATION SCHEME
MODULATION = MOD_OFDM; % SELECT CORRECT MODULATION SCHEME


% SYSTEM
dFs = 25e6;
dFp = 25e6;

switch(MODULATION)
    case MOD_OFDM
        fprintf('--OFDM--\n');
        demo = cDemoOFDM(dFs, dFp);
        BPFrm = demo.mod.BPSYM;
        txBits = randi([0 1],demo.mod.BPSYM,1);
        
    case MOD_OOK
        %*************************** OOK *********************************%
        fprintf('--OOK--\n');
        BPFrm = 128;
        demo = cDemoOOK(dFs, dFp, BPFrm);
        txBits = randi([0 1],BPFrm,1);
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
% CHANNEL
h = 1;
% for SHIFT = 0:NFRM-1
for SHIFT = 0:9
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
    
    % Sample transmitted frame at receiver with ADC
    rxFrmCh = updnClock(frmCh,demo.DAC.dCLKs,demo.ADC.dCLKs);
    rxFrm = demo.ADC.getOutput(rxFrmCh);
    
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
























































