classdef cFMCETHID
    methods(Static)
        function val = GetNetPCIcards(strType)
            [~,str] = system('"C:\Program Files (x86)\Windows Kits\8.1\Tools\x64\devcon" find =Net PCI\*');
            iNA = strfind(str,'4&11050A08&0&00E5');
            iTX = strfind(str,'4&18BDCD5B&0&48F0');
            iRX = strfind(str,'4&18BDCD5B&0&20F0');
            idxs = sort([iNA;iTX;iRX],1,'descend');
            switch(lower(strType))
                case 'tx'
                val = find(idxs == iTX)-1;
                case 'rx'
                val = find(idxs == iRX)-1;
            end
        end
    end
end