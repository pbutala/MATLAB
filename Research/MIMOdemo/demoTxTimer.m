function demoTxTimer(obj,event)
global FIGTITLE FIGTX demo BPFrm SYNC datBits;
global FIGWID FIGHGT TXSCR SCRSZ;
global FIGXOFF FIGYOFF FIGXDLT FIGYDLT;
persistent TXCLT;
persistent CHNLWRITE;

switch(lower(event.Type))
    case{'startfcn'}
        % start tcp server on interface
        cmd = sprintf('start C:\\ProgramData\\_4DSP_Training\\FMC204\\Debug\\Fmc204APP.exe 1 ML605 %d %d',...
            demo.DAC.dETHID,demo.DAC.dCLKSRC);
        dos(cmd);
        pause(2); % wait till the interface app initializes
        
        % start tcp client and Connect to FMC204 interface
        TXCLT=tcpip('localhost', demo.DAC.SKT_PORT, 'NetworkRole', 'client');
        TXCLT.OutputBufferSize = demo.frmTxNSmp8;
        fopen(TXCLT);
        
        % set BurstSize for data to transmit
        fwrite(TXCLT,[demo.DAC.CMD_BURSTSIZE demo.DAC.CHNL_ALL demo.DAC.LEN_BS_LSB demo.DAC.LEN_BS_MSB]);
        fwrite(TXCLT,typecast(int16(demo.frmTxNSmp16),'uint8'));
        
        % Initialize persistent variables
        CHNLWRITE = demo.DAC.CHNL_4; % so in timer, channel 1 is selected first
        
        % Create figures to plot transmitted signal and processing chain
        for i=1:4
            L = FIGXOFF + TXSCR*SCRSZ(3) + rem(i-1,2)*(FIGWID+FIGXDLT);
            B = FIGYOFF + SCRSZ(4) - (floor((i-1)/2)+1)*(FIGHGT+FIGYDLT);
            FIGTX(i) = figure('Name', sprintf('Frame - Transmit %d (%d Msps)', i, demo.DAC.dCLKs/1e6), 'NumberTitle', FIGTITLE,...
                            'Position',[L,B,FIGWID,FIGHGT]);
        end
        
    case{'stopfcn'}
        fclose(TXCLT);
        delete(obj);
        delete(TXCLT);
        clear TXCLT;
        
    case{'timerfcn'}
        if (SYNC == 1)  % if sync state == transmit
            % Set channel to transmit on
            switch(CHNLWRITE)
                case{demo.DAC.CHNL_1}
                    CHNLWRITE = demo.DAC.CHNL_2;
                case{demo.DAC.CHNL_2}
                    CHNLWRITE = demo.DAC.CHNL_3;
                case{demo.DAC.CHNL_3}
                    CHNLWRITE = demo.DAC.CHNL_4;
                case{demo.DAC.CHNL_4}
                    CHNLWRITE = demo.DAC.CHNL_1;
            end
            CHNLWRITEIDX = log2(double(CHNLWRITE))+1;
            
            % Queue bits to transmit in modulator
            txBits = randi([0 1],BPFrm,1);
            demo.mod.write(txBits);
            datBits(CHNLWRITEIDX).enQ(txBits);
            
            % Read modulated signal
            txSig = demo.mod.read(demo.mod.COUNTOUT);
            txSig(txSig>demo.DAC.dSIGMAX) = demo.DAC.dSIGMAX;
            txSig(txSig<demo.DAC.dSIGMIN) = demo.DAC.dSIGMIN;
            
            % Generate transmit frame
            txFrm = [demo.DAC.setRail2Rail(demo.plt.PILOT); txSig];
            
            % Plot transmit frame
            figure(FIGTX(CHNLWRITEIDX));
            plot(1:demo.frmTxNSmp16, txFrm);
            axis([1 demo.frmTxNSmp16  min(txFrm) max(txFrm)]);
            xlabel('Normalized time');
            ylabel('Signal value');
            title('Transmit frame');
            
            % Write data to channel 1
            fwrite(TXCLT,[demo.DAC.CMD_DATA CHNLWRITE typecast(int16(demo.frmTxNSmp8),'uint8')]);
            fwrite(TXCLT,typecast(int16(txFrm),'uint8'));
            
            if (CHNLWRITE == demo.ADC.CHNL_4)
                % Enable channels
                fwrite(TXCLT,[demo.DAC.CMD_ENCHNL demo.DAC.CHNL_ALL demo.DAC.ZERO_UC demo.DAC.ZERO_UC]);    % All channels are enabled irrespective of CHNL selection
                
                % Arm DAC
                fwrite(TXCLT,[demo.DAC.CMD_ARMDAC demo.DAC.CHNL_ALL demo.DAC.ZERO_UC demo.DAC.ZERO_UC]);
                SYNC = 0;       % Set sync state to receive
            end
        end
    otherwise
        error('Unknown timer event: %s',char(event.Type));
end