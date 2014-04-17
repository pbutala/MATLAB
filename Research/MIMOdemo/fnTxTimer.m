function fnTxTimer(obj,event)
global S TXCLT TXRUN TXFRAME pilot;
% global FIGTXDPL FIGTITLE;
% persistent ITER;
switch(lower(event.Type))
    case{'startfcn'}
        %         ITER = 1;
    case{'stopfcn'}
        fclose(TXCLT);
        delete(obj);
        delete(TXCLT);
        clear TXCLT;
    case{'timerfcn'}
        % Generate Data
        data = randi(S.modM,[S.frmDATALEN16 1]);
        
        % Generate OFDM symbol
        txSig = genOFDMsignal(... % Variable Arguments to the function
            'data',data,...
            'OFDMtype',S.modTYPE,...
            'N',S.modN,...
            'Symbols',getQAMsyms(S.modM,true),...
            'ClipLow',S.modCLIPL,...
            'ClipHigh',S.modCLIPH,...
            'OffsetDcoStddev', S.modOFST,...
            'ShowConst',false);
        
        % Scale OFDM symbol
        signal = txSig/S.modCLIPH;  % signal max is 1
        
        % Generate frame
        frame = [pilot(:);signal(:)];
        frIP = updnClock(frame,S.dFs,S.dCLKs,true);
        
        % Ensure Pilot is between 0 and 1
        frIP = frIP/(max(frIP(1:S.frmPILOTLEN16SF))- min(frIP(1:S.frmPILOTLEN16SF)));
        frIP = frIP - min(frIP(1:S.frmPILOTLEN16SF));
        
        % Shift, Scale
        TXFRAME = frIP*S.dSIGMAX;
        TXFRAME(TXFRAME > S.dSIGMAX) = S.dSIGMAX;
        TXFRAME(TXFRAME < S.dSIGMIN) = S.dSIGMIN;
        
        fwrite(TXCLT,typecast(int16(TXFRAME),'uint8'));
%         axis tight;
        if ~TXRUN
            stop(obj);
        end
    otherwise
        error('Unknown timer event: %s',char(event.Type));
end
end