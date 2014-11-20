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
dac = cFMC204_ETH();
adc = cFMC116_ETH();
% plt = cPilotBarker('BARKER11', dFp, dac.dCLKs);

switch(MODULATION)
    case MOD_OFDM
        %************************** OFDM *********************************%
        fprintf('--OFDM--\n');
        ofdmTyp = 'DCOOFDM';            % type of OFDM
        ofst = 3.2;                     % offset applied to time domain signal
        nsc = 64;                        % # total subcarriers
        msc = 16;                        % M-QAM
        cliph = 3.2;                    % clip high
        clipl = 0;                      % clip low
        syms = getQAMsyms(msc,true);    % QAM symbols
        
%         %******************* TO SET SIGNAL SCALE *******************%
%         [cdf,bins,lo,hi] = getOFDMdist(ofdmTyp,  nsc, syms, ofst, 0);
%         I = find(cdf>0.99,1,'first');
%         SIGMAX = bins(I);
%         %***********************************************************%
        SIGMAX = 10; % 5.5102 FROM ABOVE ANALYSIS
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

figure;
plot(txSig);
title('OFDM signal UpSampled');

% txFrm = [dac.setRail2Rail(plt.PILOT); txSig];
txFrm = txSig;

% CHANNEL
h = 1;
% SHIFT = 0;
txFrmCh = dac.getOutput(txFrm);
figure;
plot(txFrmCh);
title('OFDM signal DAC out');

% txFrmCh = [txFrmCh(end-SHIFT+1:end); txFrmCh(1:end-SHIFT)];
% figure;
% plot(txFrmCh);
frmCh = h*txFrmCh;
rxFrmCh = updnClock(frmCh,dac.dCLKs,adc.dCLKs,UPDNFLTTYP,false);
figure;
plot(rxFrmCh);
title('OFDM signal ADC in');

% RECEIVE
rxFrm = adc.getOutput(rxFrmCh);
figure;
plot(rxFrm);
title('OFDM signal ADC out');
% rxFrmUP = updnClock(rxFrm,adc.dCLKs,dac.dCLKs,UPDNFLTTYP,false);
% idx = plt.alignPilot(rxFrmUP, dac.dCLKs);

% rxPltUS = rxFrmUP(idx:idx+plt.LENGTH-1);
% rxSigUS = [rxFrmUP(idx+plt.LENGTH:end); rxFrmUP(1:idx-1)];
% rxSig = updnClock(rxSigUS,dac.dCLKs,adc.dCLKs,UPDNFLTTYP,false);
% rxSig = updnClock(rxFrmUP,dac.dCLKs,adc.dCLKs,UPDNFLTTYP,false);

scl = 1/(h*mod.SCALE*dac.dGAIN*adc.dGAIN); % this seems correct

rxFrm = rxFrm*scl;
figure;
plot(rxFrm);
title('OFDM signal ADC out and scaled');

rxFrm = rxFrm - min(rxFrm);
demod.write(rxFrm);
rxBits = demod.read(demod.COUNTOUT);

% BERR
berr = biterr2(txBits,rxBits);
fprintf('BERRs = %d\n',berr);
























































