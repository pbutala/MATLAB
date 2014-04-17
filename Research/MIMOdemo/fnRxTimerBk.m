function fnRxTimer(obj,event)
global S Stx RXCLT RXRUN RXFRAME RXDATASCL datab rxDatb;
persistent FIGSIG FIGTITLE;
persistent ITER;
switch(lower(event.Type))
    case{'startfcn'}
        fwrite(RXCLT,'a');
        ITER = 1;
        FIGTITLE = 'Off';
    case{'stopfcn'}
        fclose(RXCLT);
        delete(obj);
        delete(RXCLT);
        clear RXCLT;
    case{'timerfcn'}
        BYTECOUNT = RXCLT.BytesAvailable;
%         fprintf('BYTECOUNT = %d\n',BYTECOUNT);
        if (BYTECOUNT >= S.frmLEN16SF)
            % READ DATA FROM RECEIVER
            frameLong = typecast(uint8(fread(RXCLT,BYTECOUNT,'uchar')),'int16');
            frame = frameLong(1:S.frmLEN16SF);
            frame = double(frame*-1);
            fwrite(RXCLT,'a'); % Start another acquisition
            
%             % Take FFT of received signal
%             frmL = length(frame);
%             frmN = 2^nextpow2(frmL);
%             frmF = fft(frame,frmN)/frmL;
%             frmf = S.dCLKs/2*linspace(0,1,frmN/2+1);
%             frmf2 = S.dCLKs/2*linspace(0,2,frmN);
%             frmFH = abs(frmF);
%             
%             % Create filter
%             flt = zeros(frmN,1);
%             flt(frmf2<S.dFs/2) = 1;
%             flt(frmf2>S.dCLKs-S.dFs/2) = 1;
%             
%             % low pass filter the signal to prevent aliasing
%             frmI = frmF(flt==1);
%             frmIM = abs(frmI);
%             iL = length(frmI);
%             frf = S.dFs/2*linspace(0,1,iL/2+1);
%             frmIP = ifft(frmI,iL,1,'symmetric');
            frmIP = updnClock(frame,S.dCLKs,S.dFs,true);
            RXFRMSHFT = frmIP(1:S.frmLEN16);
            % append 
            frmcyc = [RXFRMSHFT; RXFRMSHFT(1:floor(S.frmPILOTLEN16))];
            frmcyc = frmcyc - min(frmcyc);
            Pilot = S.getPILOT01();
            [acor, lag] = xcorr(Pilot,frmcyc);
            [~,I] = max(acor(lag<=0));
            pSTART = rem(abs(lag(I)),S.frmLEN16)+1;
            pSTOP = rem(pSTART+S.frmPILOTLEN16-1,S.frmLEN16);
            
            % FRAME RECONSTRUCT
            RXFRAME = [RXFRMSHFT(pSTART:end);RXFRMSHFT(1:pSTART-1)];
            RXPILOT = RXFRAME(1:S.frmPILOTLEN16);
            RXDATA = (RXFRAME(S.frmPILOTLEN16+1:end));
            
            % CALCULATE OFFSET AND GAIN
            CHNLOFST = min(RXPILOT);
            CHNLGAIN = (max(RXPILOT)-min(RXPILOT))/Stx.dSIGMAX;
            
            % SCALE RECEIVED SIGNAL FOR DETECTION
            RXFRAMESCL = (RXFRAME-CHNLOFST)/CHNLGAIN;
            RXDATASCL = S.modCLIPL + ((RXDATA-CHNLOFST)*(S.modCLIPH-S.modCLIPL))/(CHNLGAIN*Stx.dSIGMAX);
            
            % DECODE OFDM
            rxDat = decodeOFDMsignal(RXDATASCL,...
                'OFDMtype',S.modTYPE,...    
                'N',S.modN,...
                'Symbols',getQAMsyms(S.modM));
            rxDatb = dec2binMat(rxDat-1,log2(S.modM));
            if ~isequal(size(rxDatb),[S.frmDATALEN16 log2(S.modM)])
                error('binary data size error');
            end
            % CALCULATE BIT-ERRORS
            err = biterr2(datab,rxDatb);
            ber = err/numel(rxDatb);
            fprintf('%d. BER = %0.2e\n',ITER,ber);
            
            % FIGURES
            % Signal Analysis
            FIGSIG(end+1) = figure('Name',sprintf('%d. Signal Downsample',ITER),'NumberTitle',FIGTITLE);
            figure(FIGSIG(end));
            
%             subplot(2,3,1);
%             plot(frame);
%             axis tight;
%             title(sprintf('Received frame\nADC sample rate = %0.2f Msps',S.dCLKs*1e-6));
%             xlabel('Samples');
%             
%             subplot(2,3,4);
%             ax = plotyy(frmf*1e-6,2*frmFH(1:frmN/2+1),frmf*1e-6,flt(1:frmN/2+1));
%             axis(ax,'tight');
%             title('Single sided amplitude spectrum of received frame');
%             xlabel('Frequency (MHz)');
%             ylabel('|frmFT|');
%             
            subplot(2,3,2);
            plot(RXFRMSHFT);
            axis tight;
            hold on;
            plot([pSTART,pSTART],get(gca,'YLim'),'g');
            plot([pSTOP,pSTOP],get(gca,'YLim'),'r');
            title(sprintf('Downsampled frame\nSignal sample rate = %0.2f Msps',S.dFs*1e-6));
            xlabel('Samples');
            
%             subplot(2,3,5);
%             plot(frf*1e-6,2*frmIM(1:iL/2+1));
%             axis tight;
%             title('Single sided amplitude spectrum of downsampled frame');
%             xlabel('Frequency (MHz)');
%             ylabel('|frmI|');
            
            subplot(2,3,3);
            plot(RXFRAMESCL);
            title(sprintf('Signal frame reconstructed after pilot detection.\nBit error fraction = %0.2e',ber));
            axis tight;
            xlabel('Samples');
            
            subplot(2,3,6);
            plot(lag,acor);
            axis tight;
            title(sprintf('Correlation between downsampled received signal and pilot.\nPilot start index = %d',pSTART));
            
%             % Signal Analysis
%             FIGRXSIG(end+1) = figure('Name',sprintf('%d. Received Signal',ITER),'NumberTitle',FIGTITLE);
%             plot(frame);
%             axis tight;
            % Increment ITER
            ITER = ITER + 1;
        else
            obj.TasksToExecute = obj.TasksToExecute + 1;
            fprintf('Incremented timer tasks to execute\n');
        end
        
        if ~RXRUN
            stop(obj);
        end
    otherwise
        error('Unknown timer event: %s',char(event.Type));
end
end








