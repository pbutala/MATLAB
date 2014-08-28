classdef cChnlState < handle
    properties
        SYMS;
        TxIdx;
        RxIdx;
        RxSymEst;
        ChH;
        SNRdB;
        BER;
    end % properties
    
    methods
        function obj = cChnlState(arrLen,szRxSyms,syms,chh)
            if exist('StSz','var')
                obj.TxIdx = nan(1,arrLen);
                obj.RxIdx = nan(1,arrLen);
            else
                obj.TxIdx = [];
                obj.RxIdx = [];
            end
            
            if exist('szRxSyms','var')
                obj.RxSymEst = nan(szRxSyms);
            else
                obj.RxSymEst = [];
            end
            
            if exist('syms','var')
                obj.SYMS = syms;
            else
                obj.SYMS = [];
            end
            
            if exist('chh','var')
                obj.ChH = chh;
            else
                obj.ChH = [];
            end
            obj.SNRdB = [];
            obj.BER = [];
        end % cChnlState
    end % methods
end