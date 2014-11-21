close all;
clearvars;
clc;
rng('Default');

% CONSTANT
MOD_OFDM = 1;
MOD_OOK = 2;
MODULATION = MOD_OOK; % SELECT CORRECT MODULATION SCHEME
% MODULATION = MOD_OFDM; % SELECT CORRECT MODULATION SCHEME


% SYSTEM
dFs = 25e6;
dFp = 25e6;
dac = cFMC204_ETH();
adc = cFMC116_ETH();
plt = cPilotBarker('BARKER13', dFp, dac.dCLKs);

switch(MODULATION)
    case MOD_OFDM
        %************************** OFDM *********************************%
        fprintf('--OFDM--\n');
        ofdmTyp = 'DCOOFDM';            % type of OFDM
        ofst = 3.2;                     % offset applied to time domain signal
        nsc = 8;                       % # total subcarriers
        msc = 4;                       % M-QAM
        cliph = 3.2;                    % clip high
        clipl = 0;                      % clip low
        syms = getQAMsyms(msc,true);    % QAM symbols
        
%         %******************* TO SET SIGNAL SCALE *******************%
%         [cdf,bins,lo,hi] = getOFDMdist(ofdmTyp, nsc, syms, ofst, 0);
%         I = find(cdf>0.99,1,'first');
%         SIGMAX = bins(I);
        %***********************************************************%
        % SIGMAX | OFDM | NSC  | MSC  | cdf>TH |
        % 5.5102 |  DCO |  64  |  16  |   0.99 |
        %        |  DCO |   8  |   4  |   0.99 |
        SIGMAX = 10; 
        SIGMIN = 0;
        SIGSCALE = (dac.dSIGMAX - dac.dSIGMIN)/(SIGMAX-SIGMIN);
        mod = cModOFDM('DCOOFDM', nsc, syms, ofst, SIGMIN, SIGMAX,...
            SIGSCALE,...
            dac.dSIGMIN - SIGMIN,...
            dFs, dac.dCLKs);
        
        demod = cDemodOFDM('DCOOFDM', nsc, syms, ofst, adc.dCLKs, dFs);
        UPDNFLTTYP = 'IDEALRECT';
        txBits = randi([0 1],mod.BPSYM,1);
        
    case MOD_OOK
        %*************************** OOK *********************************%
        fprintf('--OOK--\n');
        BPFrm = 128;
        SIGSCALE = (dac.dSIGMAX - dac.dSIGMIN);
        mod = cModOOK(dac.dSIGMAX, dac.dSIGMIN, dFs, dac.dCLKs,...
            SIGSCALE, dac.dSIGMIN, BPFrm);
        demod = cDemodOOK(1, 0, adc.dCLKs, dFs, BPFrm*ceil(adc.dCLKs/dFs));
        UPDNFLTTYP = 'RAISEDCOSINE';
        txBits = randi([0 1],BPFrm,1);
    otherwise
        error('Modulation not defined');
end
mod.FILTER = UPDNFLTTYP;

% TRANSMIT
mod.write(txBits);
txSig = mod.read(mod.COUNTOUT);
txSig(txSig>dac.dSIGMAX) = dac.dSIGMAX;
txSig(txSig<dac.dSIGMIN) = dac.dSIGMIN;

txFrm = [dac.setRail2Rail(plt.PILOT); txSig];
% txFrm = dac.setRail2Rail(plt.PILOT);
NFRM = numel(txFrm);
% CHANNEL
h = 1;

% for SHIFT = 0:NFRM-1
for SHIFT = 0:9
    % Transmit frame with DAC
    txFrmCh = dac.getOutput(txFrm);
    
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
    rxFrmCh = updnClock(frmCh,dac.dCLKs,adc.dCLKs);
    rxFrm = adc.getOutput(rxFrmCh);
    
    % Upsample received frame to transmit (DAC) clock for better alignment
    rxFrmUP = updnClock(rxFrm,adc.dCLKs,dac.dCLKs);
    
    % Find frame starting index by aligning pilot
    idx = plt.alignPilot(rxFrmUP-min(rxFrmUP), dac.dCLKs);
    
    % Find pilot in up-sampled received frame
    pltIs = idx:idx+plt.LENGTH-1;
    while ~isempty(pltIs(pltIs<=0))
        pltIs(pltIs<=0) = pltIs(pltIs<=0) + NFRM;
    end
    while ~isempty(pltIs(pltIs>NFRM))
        pltIs(pltIs>NFRM) = pltIs(pltIs>NFRM) - NFRM;
    end
    rxPltUS = rxFrmUP(pltIs);
    
    % Estimate channel gain from pilot
    hhat = plt.getScale(rxPltUS, dac.dCLKs);
    scl = 1/(hhat*mod.SCALE*dac.dGAIN);
    
    % Find signal in up-sampled received frame
    if pltIs(end) < pltIs(1) % if index roll over
        sigIs = pltIs(end)+1:pltIs(1)-1;
    else
        sigIs = [pltIs(end)+1:NFRM, 1:pltIs(1)-1];
    end
    
    rxSigUS = rxFrmUP(sigIs);
    rxSig = updnClock(rxSigUS,dac.dCLKs,adc.dCLKs);
    
    % Scale received signal
    rxSig = rxSig*scl;
    
    % Demodulate received signal
    rxSig = rxSig - min(rxSig);
    demod.write(rxSig);
    rxBits = demod.read(demod.COUNTOUT);
    
    % Find bit errors
    berr = biterr2(txBits,rxBits);
    
    % Find align error (for known shifts) and print
    % TODO: Block enabled only for simulated hardware
    alnErr = idx-1-SHIFT;
    alnErr(alnErr<=-NFRM/2) = alnErr(alnErr<=-NFRM/2) + NFRM;
    fprintf('Shift = %4d, BERRs = %3d, ALNERR = %4d\n',SHIFT,berr,alnErr);
%     if(abs(alnErr) >= 1/(2*demod.SMPLF))
%         fprintf('Shift = %d, idx = %d, err = %d\n',SHIFT,idx,alnErr);
%     end
end
























































