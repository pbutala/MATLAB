function demoTxTimer(obj,event)
global FIGTITLE demo txBits;
persistent TXCLT;

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
        
    case{'stopfcn'}
        fclose(TXCLT);
        delete(obj);
        delete(TXCLT);
        clear TXCLT;
        
    case{'timerfcn'}
        % Queue bits to transmit in modulator
        demo.mod.write(txBits);
        
        % Read modulated signal
        txSig = demo.mod.read(demo.mod.COUNTOUT);
        txSig(txSig>demo.DAC.dSIGMAX) = demo.DAC.dSIGMAX;
        txSig(txSig<demo.DAC.dSIGMIN) = demo.DAC.dSIGMIN;
        
        % Generate transmit frame
        txFrm = [demo.DAC.setRail2Rail(demo.plt.PILOT); txSig];
        
        % Plot transmit frame
        figure('Name', sprintf('Frame - Transmit (%d Msps)',demo.DAC.dCLKs/1e6), 'NumberTitle', FIGTITLE);
        plot(1:demo.frmTxNSmp16, txFrm);
        axis([1 demo.frmTxNSmp16  min(txFrm) max(txFrm)]);
        xlabel('Normalized time');
        ylabel('Signal value');
        title('Transmit frame');
        
        % Write data to channel 1
        fwrite(TXCLT,[demo.DAC.CMD_DATA demo.DAC.CHNL_1 typecast(int16(demo.frmTxNSmp8),'uint8')]);
        fwrite(TXCLT,typecast(int16(txFrm),'uint8'));
        
        % Write data to channel 2
        fwrite(TXCLT,[demo.DAC.CMD_DATA demo.DAC.CHNL_2 typecast(int16(demo.frmTxNSmp8),'uint8')]);
        fwrite(TXCLT,typecast(int16(txFrm),'uint8'));
        
        % Write data to channel 3
        fwrite(TXCLT,[demo.DAC.CMD_DATA demo.DAC.CHNL_3 typecast(int16(demo.frmTxNSmp8),'uint8')]);
        fwrite(TXCLT,typecast(int16(txFrm),'uint8'));
        
        % Write data to channel 4
        fwrite(TXCLT,[demo.DAC.CMD_DATA demo.DAC.CHNL_4 typecast(int16(demo.frmTxNSmp8),'uint8')]);
        fwrite(TXCLT,typecast(int16(txFrm),'uint8'));
        
        % Enable channels
        fwrite(TXCLT,[demo.DAC.CMD_ENCHNL demo.DAC.CHNL_4 demo.DAC.ZERO_UC demo.DAC.ZERO_UC]);    % All channels are enabled irrespective of CHNL selection
        
        % Arm DAC
        fwrite(TXCLT,[demo.DAC.CMD_ARMDAC demo.DAC.CHNL_ALL demo.DAC.ZERO_UC demo.DAC.ZERO_UC]);
        
    otherwise
        error('Unknown timer event: %s',char(event.Type));
end