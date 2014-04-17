function fnRxTimer(obj,event)
global S Stx RXCLT RXRUN RXFRAME RXDATASCL rxDatb;
persistent FIGSIG FIGTITLE AXSYM;
persistent ITER;
switch(lower(event.Type))
    case{'startfcn'}
        fwrite(RXCLT,'a');
        ITER = 1;
        FIGTITLE = 'Off';
        figure('Name',sprintf('Received Symbols'),'NumberTitle',FIGTITLE);
        AXSYM = gca;
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
            %**************************************************************
            rxFrUP = updnClock(frame,S.dCLKs,Stx.dCLKs,false);
            RXFRMSHFT = rxFrUP(1:Stx.frmLEN16SF);
            frmcyc = [RXFRMSHFT; RXFRMSHFT(1:floor(Stx.frmPILOTLEN16SF))];
            frmcyc = frmcyc - min(frmcyc);
            Pilot = updnClock(S.getPILOT01(),S.dFs,Stx.dCLKs,false);
            [acor, lag] = xcorr(Pilot,frmcyc);
            [~,I] = max(acor(lag<=0));
            pSTART = rem(abs(lag(I)),Stx.frmLEN16SF)+1;
            pSTOP = rem(pSTART+Stx.frmPILOTLEN16SF-1,Stx.frmLEN16SF);
            % FRAME RECONSTRUCT
            rxFrUP0 = [RXFRMSHFT(pSTART:end);RXFRMSHFT(1:pSTART-1)];
            %**************************************************************
            RXFRAME = updnClock(rxFrUP0,Stx.dCLKs,S.dFs,false);
            RXPILOT = RXFRAME(1:S.frmPILOTLEN16);
            RXDATA = (RXFRAME(S.frmPILOTLEN16+1:end));
            
            % CALCULATE OFFSET AND GAIN
            CHNLOFST = min(RXPILOT);
            CHNLGAIN = (max(RXPILOT)-min(RXPILOT))/Stx.dSIGMAX;
            
            % SCALE RECEIVED SIGNAL FOR DETECTION
            RXFRAMESCL = (RXFRAME-CHNLOFST)/(CHNLGAIN*Stx.dSIGMAX);
            RXDATASCL = S.modCLIPL + ((RXDATA-CHNLOFST)*(S.modCLIPH-S.modCLIPL))/(CHNLGAIN*Stx.dSIGMAX);
            
            % DECODE OFDM
            rxDat = decodeOFDMsignal(RXDATASCL,...
                'OFDMtype',S.modTYPE,...    
                'N',S.modN,...
                'Symbols',getQAMsyms(S.modM,true),...
                'ShowRcv',true,...
                'ShowRcvAx',AXSYM);
            rxDatb = dec2binMat(rxDat-1,log2(S.modM));
            rng('default'); % TODO DELETE WHEN DIFFERENT DATA IN DIFF FRAMES
            data = randi(S.modM,[S.frmDATALEN16 1]);            % transmitted data
            datab = dec2binMat(data-1);
            % CALCULATE BIT-ERRORS
            err = biterr2(datab,rxDatb);
            ber = err/numel(rxDatb);
            fprintf('%d. BER = %0.2e\n',ITER,ber);
            
%             % FIGURES
%             % Signal Analysis
%             FIGSIG(end+1) = figure('Name',sprintf('%d. Signal Downsample',ITER),'NumberTitle',FIGTITLE);
%             figure(FIGSIG(end));
%            
%             subplot(2,3,2);
%             plot(RXFRMSHFT);
%             axis tight;
%             hold on;
%             plot([pSTART,pSTART],get(gca,'YLim'),'g');
%             plot([pSTOP,pSTOP],get(gca,'YLim'),'r');
%             title(sprintf('Downsampled frame\nSignal sample rate = %0.2f Msps',S.dFs*1e-6));
%             xlabel('Samples');
%             
%             subplot(2,3,3);
%             plot(RXFRAMESCL);
%             title(sprintf('Signal frame reconstructed after pilot detection.\nBit error fraction = %0.2e',ber));
%             axis tight;
%             xlabel('Samples');
%             
%             subplot(2,3,6);
%             plot(lag,acor);
%             axis tight;
%             title(sprintf('Correlation between downsampled received signal and pilot.\nPilot start index = %d',pSTART));
            
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








