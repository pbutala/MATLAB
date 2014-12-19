classdef cCSK < handle
    % class to handle Color Shift Keying
    properties(Constant)
        CB0 = struct('LMin',380,'LMax',478,'Code','000','Center',429,'x',0.169,'y',0.007); % Color Band 0
        CB1 = struct('LMin',478,'LMax',540,'Code','001','Center',509,'x',0.011,'y',0.733); % Color Band 1
        CB2 = struct('LMin',540,'LMax',588,'Code','010','Center',564,'x',0.402,'y',0.597); % Color Band 2
        CB3 = struct('LMin',588,'LMax',633,'Code','011','Center',611,'x',0.669,'y',0.331); % Color Band 3
        CB4 = struct('LMin',633,'LMax',679,'Code','100','Center',656,'x',0.729,'y',0.271); % Color Band 4
        CB5 = struct('LMin',679,'LMax',726,'Code','101','Center',703,'x',0.734,'y',0.265); % Color Band 5
        CB6 = struct('LMin',726,'LMax',780,'Code','110','Center',753,'x',0.734,'y',0.265); % Color Band 6
    end % properties - constant
    
    properties(SetAccess = immutable)
        CBCid;
        CBC;
    end % properties - immutable
    
    methods
        function obj = cCSK(iCBC)
            switch iCBC
                case 1 % CBC 1
                    obj.CBC = [obj.CB6 obj.CB2 obj.CB0]; % 6;2;0
                case 2 % CBC 2
                    obj.CBC = [obj.CB6 obj.CB1 obj.CB0]; % 6;1;0
                case 3 % CBC 3
                    obj.CBC = [obj.CB5 obj.CB2 obj.CB0]; % 5;2;0
                case 4 % CBC 4
                    obj.CBC = [obj.CB5 obj.CB1 obj.CB0]; % 5;1;0
                case 5 % CBC 5
                    obj.CBC = [obj.CB4 obj.CB2 obj.CB0]; % 4;2;0
                case 6 % CBC 6
                    obj.CBC = [obj.CB4 obj.CB1 obj.CB0]; % 4;1;0
                case 7 % CBC 7
                    obj.CBC = [obj.CB3 obj.CB2 obj.CB0]; % 3;2;0
                case 8 % CBC 8
                    obj.CBC = [obj.CB3 obj.CB1 obj.CB0]; % 3;1;0
                case 9 % CBC 9
                    obj.CBC = [obj.CB2 obj.CB1 obj.CB0]; % 2;1;0
                otherwise
                    error('CBC-%d not supported',iCBC);
            end
            obj.CBCid = iCBC;
        end % cCSK
    end % methods
end % cCSK