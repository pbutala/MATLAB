classdef cCBC < handle
    % class to handle Color Band Combinations according to IEEE 802.15.7
    properties(Constant, GetAccess=private)
        CB0 = struct('id',0,'LMin',380,'LMax',478,'Code','000','Center',429,'x',0.169,'y',0.007); % Color Band 0
        CB1 = struct('id',1,'LMin',478,'LMax',540,'Code','001','Center',509,'x',0.011,'y',0.733); % Color Band 1
        CB2 = struct('id',2,'LMin',540,'LMax',588,'Code','010','Center',564,'x',0.402,'y',0.597); % Color Band 2
        CB3 = struct('id',3,'LMin',588,'LMax',633,'Code','011','Center',611,'x',0.669,'y',0.331); % Color Band 3
        CB4 = struct('id',4,'LMin',633,'LMax',679,'Code','100','Center',656,'x',0.729,'y',0.271); % Color Band 4
        CB5 = struct('id',5,'LMin',679,'LMax',726,'Code','101','Center',703,'x',0.734,'y',0.265); % Color Band 5
        CB6 = struct('id',6,'LMin',726,'LMax',780,'Code','110','Center',753,'x',0.734,'y',0.265); % Color Band 6
    end % properties - constant, get=private
    
    properties(Constant)
        CBCNT = 7;
    end
    properties(SetAccess = private)
        CBs;
    end
    methods
        function obj = cCBC()
            obj.CBs = [obj.CB0,obj.CB1,obj.CB2,obj.CB3,obj.CB4,obj.CB5,obj.CB6];
        end
        function CBs = getColorBands(obj,iCBC)
            switch iCBC
                case 1 % CBC 1
                    CBs = [obj.CB6 obj.CB2 obj.CB0]; % 6;2;0
                case 2 % CBC 2
                    CBs = [obj.CB6 obj.CB1 obj.CB0]; % 6;1;0
                case 3 % CBC 3
                    CBs = [obj.CB5 obj.CB2 obj.CB0]; % 5;2;0
                case 4 % CBC 4
                    CBs = [obj.CB5 obj.CB1 obj.CB0]; % 5;1;0
                case 5 % CBC 5
                    CBs = [obj.CB4 obj.CB2 obj.CB0]; % 4;2;0
                case 6 % CBC 6
                    CBs = [obj.CB4 obj.CB1 obj.CB0]; % 4;1;0
                case 7 % CBC 7
                    CBs = [obj.CB3 obj.CB2 obj.CB0]; % 3;2;0
                case 8 % CBC 8
                    CBs = [obj.CB3 obj.CB1 obj.CB0]; % 3;1;0
                case 9 % CBC 9
                    CBs = [obj.CB2 obj.CB1 obj.CB0]; % 2;1;0
                otherwise
                    error('CBC-%d not supported',iCBC);
            end
        end % getColorBands
    end % methods
end % cCBC