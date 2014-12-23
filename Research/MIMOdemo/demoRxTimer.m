function demoRxTimer(obj,event)
global FIGTITLE FIGRX FIGNR FIGNC demo SYNC datBits;
global FIGWID FIGHGT RXSCR SCRSZ;
global FIGXOFF FIGYOFF FIGXDLT FIGYDLT;
persistent RXCLT BS;
persistent CHNLREAD BERs CHNLRDCNT;
persistent FIGTXTH;

switch(lower(event.Type))
    case{'startfcn'}
        % start tcp server on interface
        cmd = sprintf('start C:\\ProgramData\\_4DSP_Training\\FMC116\\Debug\\Fmc116APP.exe 1 ML605 %d %d',...
            demo.ADC.dETHID,demo.ADC.dCLKSRC);
        dos(cmd);
        pause(2);                                                                   % wait till the interface app initializes
        
        % create tcp client
        RXCLT = tcpip('localhost', demo.ADC.SKT_PORT, 'NetworkRole', 'client');
        
        % set BurstSize for data to receive
        BS = 32*ceil(demo.frmRxNSmp16/32);                                           % BurstSize for FMC116 must be an integer multiple of 32
        
        % Connect to FMC116 interface
        RXCLT.InputBufferSize = BS*2;                                               % FMC116 BurstSize*2 (for int16)
        fopen(RXCLT);
        
        % set BurstSize for data to receive
        fwrite(RXCLT,[demo.ADC.CMD_BURSTSIZE demo.ADC.CHNL_1 demo.ADC.LEN_BS_LSB demo.ADC.LEN_BS_MSB]);
        fwrite(RXCLT,typecast(int16(BS),'uint8'));
        
        % initialize persistent variables
        BERs = zeros(demo.ADC.dCHNLCNT,1);
        CHNLRDCNT = zeros(demo.ADC.dCHNLCNT,1);
        
        % Start acquisition routine for 1st channel
        CHNLREAD = demo.ADC.CHNL_4;
        
        % subplot dimensions
        FIGNR = 2;
        FIGNC = 2;
        
        % Create figures to plot received signal and processing chain
        for i=1:4
            L = FIGXOFF + RXSCR*SCRSZ(3) + rem(i-1,2)*(FIGWID+FIGXDLT);
            B = FIGYOFF + SCRSZ(4) - (floor((i-1)/2)+1)*(FIGHGT+FIGYDLT);
            FIGRX(i) = figure('Name', sprintf('Demo Receive %d',i), 'NumberTitle', FIGTITLE,...
                'Position',[L,B,FIGWID,FIGHGT]);
            
%             subplot(FIGNR,FIGNC,2);
%             axis off;
%             FIGTXTH(i) = text(0,0,sprintf('BER = %0.1g',0));
        end
        
        
    case{'stopfcn'}
        fclose(RXCLT);
        delete(obj);
        delete(RXCLT);
        clear RXCLT;
        
    case{'timerfcn'}
        try
            if (SYNC == 0)  % if sync state == receive
                % Set channel to read
                switch(CHNLREAD)
                    case{demo.ADC.CHNL_1}
                        CHNLREAD = demo.ADC.CHNL_2;
                    case{demo.ADC.CHNL_2}
                        CHNLREAD = demo.ADC.CHNL_3;
                    case{demo.ADC.CHNL_3}
                        CHNLREAD = demo.ADC.CHNL_4;
                    case{demo.ADC.CHNL_4}
                        CHNLREAD = demo.ADC.CHNL_1;
                end
                CHNLREADIDX = log2(double(CHNLREAD))+1;
                fwrite(RXCLT,[demo.ADC.CMD_DATA CHNLREAD demo.ADC.ZERO_UC demo.ADC.ZERO_UC]);      % Read data from next channel
                
                BYTECOUNT = 0;
                while (BYTECOUNT < BS*2)
                    BYTECOUNT = RXCLT.BytesAvailable;
                end
                
                % READ DATA FROM RECEIVER
                frameLong = typecast(uint8(fread(RXCLT,BS*2,'uchar')),'int16');
                frame = frameLong(1:demo.frmRxNSmp16);
                
                % RECEIVER INVERTS SIGNAL TOO. TODO: FIX FOR FINAL DEMO!!
                switch(CHNLREAD)
                    case{demo.ADC.CHNL_1, demo.ADC.CHNL_4}
                        frame = double(frame); % receiver hardware inverts signal; reverses ADC inversion
                    case{demo.ADC.CHNL_2, demo.ADC.CHNL_3}
                        frame = double(frame*-1); % compensate for ADC inversion
                end
                CHNLRDCNT(CHNLREADIDX) = CHNLRDCNT(CHNLREADIDX) + 1;
                
                % Plot receive Frame (@ADC clock)
                figure(FIGRX(CHNLREADIDX));
                subplot(FIGNR,FIGNC,1);
                plot(1:demo.frmRxNSmp16, frame);
                axis([1 demo.frmRxNSmp16 min(frame) max(frame)]);
                xlabel('Normalized time');
                ylabel('Signal value');
                title('Receive frame (@ADC clock)');
                
                % Find frame starting index by aligning pilot
                idx = demo.plt.alignPilot(frame-min(frame), demo.ADC.dCLKs);
                
                % Find pilot in received frame
                NFRM = demo.frmRxNSmp16;
                pltIs = idx:idx+demo.plt.getLength(demo.ADC.dCLKs)-1;
                while ~isempty(pltIs(pltIs<=0))
                    pltIs(pltIs<=0) = pltIs(pltIs<=0) + NFRM;
                end
                while ~isempty(pltIs(pltIs>NFRM))
                    pltIs(pltIs>NFRM) = pltIs(pltIs>NFRM) - NFRM;
                end
                rxPltUS = frame(pltIs);
                
                % Show pilot start and stop cursors on figure
                figure(FIGRX(CHNLREADIDX));
                subplot(FIGNR,FIGNC,1);
                drCursor([pltIs(1), pltIs(end)],'Vertical',['g-';'r-']);
                
                % Estimate channel gain from pilot
                hhat = demo.plt.getScale(rxPltUS, demo.ADC.dCLKs);
                scl = 1/(hhat*demo.mod.SCALE*demo.DAC.dGAIN);
                
                % Find signal in frame
                if pltIs(end) < pltIs(1) % if index roll over
                    sigIs = pltIs(end)+1:pltIs(1)-1;
                else
                    sigIs = [pltIs(end)+1:NFRM, 1:pltIs(1)-1];
                end
                rxSig = frame(sigIs);
                
                % Show signal start and stop cursors on figure
                figure(FIGRX(CHNLREADIDX));
                subplot(FIGNR,FIGNC,1);
                drCursor([sigIs(1), sigIs(end)],'Vertical',['g:';'r:']);
                
                % Plot receive signal (@fSig*2 clock)
                figure(FIGRX(CHNLREADIDX));
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
                berr = biterr2(rxBits,datBits(CHNLREADIDX).deQ(datBits(CHNLREADIDX).COUNT));
                BERs(CHNLREADIDX) = (BERs(CHNLREADIDX)*(CHNLRDCNT(CHNLREADIDX)-1) + berr)/CHNLRDCNT(CHNLREADIDX);
                fprintf('BERs(%d) = %3d\n',CHNLREADIDX,berr);
                
                % Show bit error rate on figure
%                 figure(FIGRX(CHNLREADIDX));
%                 subplot(FIGNR,FIGNC,2);
%                 set(FIGTXTH(CHNLREADIDX),'String', sprintf('BER = %0.1e',BERs(CHNLREADIDX)));
                
                if (CHNLREAD == demo.ADC.CHNL_4)
                    SYNC = 1;       % Set sync state to transmit
                end
            end
        catch
        end
    otherwise
        error('Unknown timer event: %s',char(event.Type));
end
end