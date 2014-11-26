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
    
    % Shift frame start (sampling offset at receiver)
    sftIs = NFRM-SHIFT+1:2*NFRM-SHIFT;
    while ~isempty(sftIs(sftIs<=0))
        sftIs(sftIs<=0) = sftIs(sftIs<=0) + NFRM;
    end
    while ~isempty(sftIs(sftIs>NFRM))
        sftIs(sftIs>NFRM) = sftIs(sftIs>NFRM) - NFRM;
    end
    txFrmCh = txFrmCh(sftIs);
    
%     figure;
%     plot(txFrmCh);
    
    % Scale frame (channel gain)
    frmCh = h*txFrmCh;
    
    % Sample transmitted frame at receiver with ADC
    rxFrmCh = updnClock(frmCh,demo.DAC.dCLKs,demo.ADC.dCLKs);
    rxFrm = demo.ADC.getOutput(rxFrmCh);
    
    % Upsample received frame to transmit (DAC) clock for better alignment
    rxFrmUP = updnClock(rxFrm,demo.ADC.dCLKs,demo.DAC.dCLKs);
    
    % Find frame starting index by aligning pilot
    idx = demo.plt.alignPilot(rxFrmUP-min(rxFrmUP), demo.DAC.dCLKs);
    
    % Find pilot in up-sampled received frame
    pltIs = idx:idx+demo.plt.LENGTH-1;
    while ~isempty(pltIs(pltIs<=0))
        pltIs(pltIs<=0) = pltIs(pltIs<=0) + NFRM;
    end
    while ~isempty(pltIs(pltIs>NFRM))
        pltIs(pltIs>NFRM) = pltIs(pltIs>NFRM) - NFRM;
    end
    rxPltUS = rxFrmUP(pltIs);
    
    % Estimate channel gain from pilot
    hhat = demo.plt.getScale(rxPltUS, demo.DAC.dCLKs);
    scl = 1/(hhat*demo.mod.SCALE*demo.DAC.dGAIN);
    
    % Find signal in up-sampled received frame
    if pltIs(end) < pltIs(1) % if index roll over
        sigIs = pltIs(end)+1:pltIs(1)-1;
    else
        sigIs = [pltIs(end)+1:NFRM, 1:pltIs(1)-1];
    end
    
    rxSigUS = rxFrmUP(sigIs);
    rxSig = updnClock(rxSigUS,demo.DAC.dCLKs,demo.ADC.dCLKs);
    
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
    fprintf('Shift = %4d, BERRs = %3d, ALNERR = %4d\n',SHIFT,berr,alnErr);
%     if(abs(alnErr) >= 1/(2*demo.demod.SMPLF))
%         fprintf('Shift = %d, idx = %d, err = %d\n',SHIFT,idx,alnErr);
%     end
end
























































