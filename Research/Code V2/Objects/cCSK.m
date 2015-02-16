classdef cCSK < handle
    % class to handle Color Shift Keying
    properties(Constant, GetAccess=private)
        CB0 = struct('LMin',380,'LMax',478,'Code','000','Center',429,'x',0.169,'y',0.007); % Color Band 0
        CB1 = struct('LMin',478,'LMax',540,'Code','001','Center',509,'x',0.011,'y',0.733); % Color Band 1
        CB2 = struct('LMin',540,'LMax',588,'Code','010','Center',564,'x',0.402,'y',0.597); % Color Band 2
        CB3 = struct('LMin',588,'LMax',633,'Code','011','Center',611,'x',0.669,'y',0.331); % Color Band 3
        CB4 = struct('LMin',633,'LMax',679,'Code','100','Center',656,'x',0.729,'y',0.271); % Color Band 4
        CB5 = struct('LMin',679,'LMax',726,'Code','101','Center',703,'x',0.734,'y',0.265); % Color Band 5
        CB6 = struct('LMin',726,'LMax',780,'Code','110','Center',753,'x',0.734,'y',0.265); % Color Band 6
    end % properties - constant, get=private
    
    properties(SetAccess = immutable)
        M;
        CBCid;
        CBC;
    end % properties - immutable
    
    methods
        function obj = cCSK(iCBC, M)
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
            if ~exist('M','var')
                M = 4;
            end
            switch M
                case{4;8;16} 
                    obj.M = M;
                otherwise
                    error('%d-CSK is not supported',M);
            end
        end % cCSK
        
        function [x,y] = getSyms(obj)
            switch obj.M
                case 4
                    [x,y]=obj.getSyms4(obj.CBC(1).x,obj.CBC(1).y,...
                        obj.CBC(2).x,obj.CBC(2).y,...
                        obj.CBC(3).x,obj.CBC(3).y);
                case 8
                    [x,y]=obj.getSyms8(obj.CBC(1).x,obj.CBC(1).y,...
                        obj.CBC(2).x,obj.CBC(2).y,...
                        obj.CBC(3).x,obj.CBC(3).y);
                case 16
                    [x,y]=obj.getSyms16(obj.CBC(1).x,obj.CBC(1).y,...
                        obj.CBC(2).x,obj.CBC(2).y,...
                        obj.CBC(3).x,obj.CBC(3).y);
                otherwise
                    error('%d-CSK is not supported',obj.M);
            end
        end % getSyms
        
        function hax = drCenters(obj, obs,fsz)
            hax = gca;
            obs.showGamut(false);
            hold on;
            if ~exist('fsz','var')
                fsz = 16;
            end
            %
            [x,y,~] = obs.getCoordinates(obj.CB0.Center);
            scatter(x, y, 80, 'o', 'k','LineWidth',2);
            str = sprintf('CB_{0}(%dnm)',obj.CB0.Center);
            text(x+0.035,y,str,'verticalalign','middle','FontSize',fsz);
            %
            [x,y,~] = obs.getCoordinates(obj.CB1.Center);
            scatter(x, y, 80, 'o', 'k','LineWidth',2);
            str = sprintf('CB_{1}(%dnm)',obj.CB1.Center);
            text(x+0.015,y,str,'verticalalign','middle','FontSize',fsz);
            %
            [x,y,~] = obs.getCoordinates(obj.CB2.Center);
            scatter(x, y, 80, 'o', 'k','LineWidth',2);
            str = sprintf('CB_{2}(%dnm)',obj.CB2.Center);
            text(x+0.015,y,str,'verticalalign','middle','FontSize',fsz);
            %
            [x,y,~] = obs.getCoordinates(obj.CB3.Center);
            scatter(x, y, 80, 'o', 'k','LineWidth',2);
            str = sprintf('CB_{3}(%dnm)',obj.CB3.Center);
            text(x-0.3,y,str,'verticalalign','middle','FontSize',fsz);
            %
            [x,y,~] = obs.getCoordinates(obj.CB4.Center);
            scatter(x, y, 80, 'o', 'k','LineWidth',2);
            str = sprintf('CB_{4}(%dnm)',obj.CB4.Center);
            text(x+0.020,y+0.045,str,'verticalalign','middle','FontSize',fsz);
            %
            [x,y,~] = obs.getCoordinates(obj.CB5.Center);
            scatter(x, y, 80, 'o', 'k','LineWidth',2);
            str = sprintf('CB_{5}(%dnm)',obj.CB5.Center);
            text(x+0.015,y,str,'verticalalign','middle','FontSize',fsz);
            %
            [x,y,~] = obs.getCoordinates(obj.CB6.Center);
            scatter(x, y, 80, 'o', 'k','LineWidth',2);
            str = sprintf('CB_{6}(%dnm)',obj.CB6.Center);
            text(x+0.01725,y-0.0525,str,'verticalalign','middle','FontSize',fsz);
        end % drCenters
        
        function hax = drConstRules(obj,fsz,hax)
            if ~exist('fsz','var')
                fsz = 16;
            end
            if ~exist('hax','var')
                figure;
                hax = gca;
            end
            x1 = 1; y1 = 0;
            x2 = 0.5; y2 = 1;
            x3 = 0; y3 = 0;
            switch obj.M
                case 4
                    [x,y] = obj.getSyms4(x1,y1,x2,y2,x3,y3);
                case 8
                    [x,y] = obj.getSyms8(x1,y1,x2,y2,x3,y3);
                case 16
                    [x,y] = obj.getSyms16(x1,y1,x2,y2,x3,y3);
                otherwise
                    error('%d-CSK is not supported',obj.M);
            end
            plot([x1,x2,x3,x1],[y1,y2,y3,y1],'b-','LineWidth',2);
            hold on;
            scatter(hax,x,y,80, 'o', 'r','LineWidth',2);
            
            % label points
            for i=1:obj.M
                str = sprintf('(%d)',i-1);
                text(x(i)+0.015,y(i),str,'verticalalign','middle','FontSize',fsz,'color','r');
            end
            text(x1,y1+0.05,'I','verticalalign','middle','FontSize',fsz,'color','b');
            text(x2-0.04,y2,'J','verticalalign','middle','FontSize',fsz,'color','b');
            text(x3-0.01,y3+0.04,'K','verticalalign','middle','FontSize',fsz,'color','b');
            % label axis ticks
            set(hax,'XTickMode','Manual','XTick',unique(sort(x)),...
                'XTickLabelMode','Manual','XTickLabel',num2str(unique(sort(x)'),'%0.2f'),...
                'YTickMode','Manual','YTick',unique(sort(y)),...
                'YTickLabelMode','Manual','YTickLabel',num2str(unique(sort(y)'),'%0.2f'));
            grid on;
            xlabel('{\Delta}x (normalized)')
            ylabel('{\Delta}y (normalized)')
        end
    end % methods
    
    methods(Access=private)
        function [x,y] = getSyms4(obj,x1,y1,x2,y2,x3,y3)
            x = [x2 0 x3 x1];
            y = [y2 0 y3 y1];
            x(2) = mean(x([1 3 4]));
            y(2) = mean(y([1 3 4]));
        end % getSyms4
        
        function [x,y] = getSyms8(obj,x1,y1,x2,y2,x3,y3)
            x = [0 0 0 0 x2 x3 0 x1];
            y = [0 0 0 0 y2 y3 0 y1];
            x(1)=(2*x(5)+x(6))/3; y(1)=(2*y(5)+y(6))/3;
            x(7)=(2*x(5)+x(8))/3; y(7)=(2*y(5)+y(8))/3;
            x(4)=(x(6)+x(8))/2; y(4)=(y(6)+y(8))/2;
            xc=(x(6)+x(5))/2; yc=(y(6)+y(5))/2;
            xd=mean([xc,x([6,4])]); yd=mean([yc,y([6,4])]);
            x(2)=(2*xd+xc)/3; y(2)=(2*yd+yc)/3;
            
            xb=(x(8)+x(5))/2; yb=(y(8)+y(5))/2;
            xa=mean([xb,x([8,4])]); ya=mean([yb,y([8,4])]);
            x(3)=(2*xa+xb)/3; y(3)=(2*ya+yb)/3;
        end % getSyms8
        
        function [x,y] = getSyms16(obj,x1,y1,x2,y2,x3,y3)
            x = zeros(1,16); y = zeros(1,16);
            x(1) = x2; y(1) = y2;
            x(9) = x1; y(9) = y1;
            x(10) = x3; y(10) = y3;
            x(7) = mean(x([1 9 10])); y(7) = mean(y([1 9 10]));
            
            x(4)=(2*x(1)+x(10))/3; y(4)=(2*y(1)+y(10))/3;
            x(11)=(x(1)+2*x(10))/3; y(11)=(y(1)+2*y(10))/3;
            
            x(6)=(2*x(1)+x(9))/3; y(6)=(2*y(1)+y(9))/3;
            x(5)=(x(1)+2*x(9))/3; y(5)=(y(1)+2*y(9))/3;
            
            x(13)=(2*x(9)+x(10))/3; y(13)=(2*y(9)+y(10))/3;
            x(16)=(x(9)+2*x(10))/3; y(16)=(y(9)+2*y(10))/3;
            
            x(2) = mean(x([1 4 6])); y(2) = mean(y([1 4 6]));
            x(3) = mean(x([4 7 11])); y(3) = mean(y([4 7 11]));
            x(8) = mean(x([5 6 7])); y(8) = mean(y([5 6 7]));
            x(12) = mean(x([11 16 10])); y(12) = mean(y([11 16 10]));
            x(14) = mean(x([5 9 13])); y(14) = mean(y([5 9 13]));
            x(15) = mean(x([7 13 16])); y(15) = mean(y([7 13 16]));
        end % getSyms16
    end % private methods   
end % cCSK







































