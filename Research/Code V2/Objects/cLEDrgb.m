classdef cLEDrgb < handle
    % class to handle RGB LED
    properties(SetAccess = immutable)
        R;      % Red PSD
        G;      % Green PSD
        B;      % Blue PSD
        Txyz;   % RGB2XYZ transformation matrix
        Trgb;   % XYZ2RGB transformation matrix
        obs;    % Standard observer (CIE XYZ 1978)
        flRes;  % color flux normalized increment resolution
    end
    
    properties(SetAccess = private)
        tn;     % tristimulus values
        xyz;    % SPD XY coordinate on CIE 1978
        hWB;    % handle waitbar
%         cct;    % correlated color temperature for each xyz(:,i) combination. using mccamy's formula. More accurate 
    end
    
    methods
        % destructor
        function delete(obj)
            if ishandle(obj.hWB)
                delete(obj.hWB);
            end
        end
        
        % constructor
        function obj = cLEDrgb(flRes,Rpsd,Gpsd,Bpsd,flCIE)
             % normalized optical flux resolution (between 0 and 1)
            if ~exist('flRes','var')
                flRes = 0.1;
            end
            if (flRes>0) && (flRes<=1)
                obj.flRes = flRes;
            else
                error('''flRes'' must be a scalar value between 0 and 1');
            end
            
            LAMBDAMIN = 200;LAMBDADELTA = 1;LAMBDAMAX = 1100;
            lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;
            % Red PSD
            if exist('Rpsd','var')
                if isa(Rpsd,'cPSD')
                    obj.R = Rpsd;
                else
                    error('''Rpsd'' must be of type cPSD');
                end
            else
                Rpsd = getSOG(627,10,1,lambdas);
                Rpsd = Rpsd/(sum(Rpsd)*LAMBDADELTA);
                obj.R = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rpsd);
            end
            % Green PSD
            if exist('Gpsd','var')
                if isa(Gpsd,'cPSD')
                    obj.G = Gpsd;
                else
                    error('''Gpsd'' must be of type cPSD');
                end
            else
                Gpsd = getSOG(530,10,1,lambdas);
                Gpsd = Gpsd/(sum(Gpsd)*LAMBDADELTA);
                obj.G = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gpsd);
            end
            % Blue PSD
            if exist('Bpsd','var')
                if isa(Bpsd,'cPSD')
                    obj.B = Bpsd;
                else
                    error('''Bpsd'' must be of type cPSD');
                end
            else
                Bpsd = getSOG(470,10,1,lambdas);
                Bpsd = Bpsd/(sum(Bpsd)*LAMBDADELTA);
                obj.B = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bpsd);
            end
            % Generate CIE 1978 object
            if ~exist('flCIE','var')
                flCIE = 'CIE1978.csv';
            end
            obj.obs = cCIE(flCIE);
           % initialize variables
            obj.tn = [nan;nan;nan]; obj.xyz = [nan;nan;nan];
           % initialize Txyz
           obj.Txyz = obj.getRGB2XYZ();
           % initialize Trgb
           obj.Trgb = obj.getXYZ2RGB();
        end
        
        % initialize the class
        function initialize(obj)
            resa = 0:obj.flRes:1;
            resLEN = numel(resa);
            s = [0;0;0];
            txyz = [0;0;0];
            h = waitbar(0,'0.00%','Name','Characterizing LED...',...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
            try
            setappdata(h,'canceling',0);
            TOTALLOOPS = resLEN^3;
            LOOPCOUNT = 0;
            for ix = 1:resLEN
                s(1) = resa(ix);
                for iy = 1:resLEN
                    s(2) = resa(iy);
                    for iz = 1:resLEN
                        s(3) = resa(iz);
                        ttn = s/sum(s,1);
                        if ~ismember(obj.tn',ttn','rows')
%                             S = (ttn(1)*obj.R) + (ttn(2)*obj.G) + (ttn(3)*obj.B);
%                             [txyz(1),txyz(2),txyz(3)] = obj.obs.getCoordinates(S.npsd);
                            XYZ = obj.Txyz*ttn;
                            obj.tn(:,end+1) = ttn;
                            obj.xyz(:,end+1) = XYZ/sum(XYZ,1);
                        end
                        LOOPCOUNT = LOOPCOUNT + 1;
                        PROG = LOOPCOUNT/TOTALLOOPS;
                        waitbar(PROG,h,sprintf('%0.2f%% done...',PROG*100));
                        if(getappdata(h,'canceling'))
                            delete(h);
                            error('Characterization aborted');
                        end
                    end
                end
            end
            waitbar(1);
            catch ex
            delete(h);
            rethrow(ex);
            end
            delete(h);
        end
%         function initialize(obj)
% %             resa = obj.flRes:obj.flRes:1;
%             resa = 0:obj.flRes:1;
%             resLEN = numel(resa);
%             s = [0;0;0];
%             txyz = [0;0;0];
%             h = waitbar(0,'0.00%','Name','Characterizing LED...',...
%                 'CreateCancelBtn',...
%                 'setappdata(gcbf,''canceling'',1)');
%             try
%             setappdata(h,'canceling',0);
%             TOTALLOOPS = resLEN^3;
%             LOOPCOUNT = 0;
%             for ix = 1:resLEN
%                 s(1) = resa(ix);
%                 for iy = 1:resLEN
%                     s(2) = resa(iy);
%                     for iz = 1:resLEN
%                         s(3) = resa(iz);
%                         ttn = s/sum(s,1);
%                         if ~ismember(obj.tn',ttn','rows')
%                             S = (ttn(1)*obj.R) + (ttn(2)*obj.G) + (ttn(3)*obj.B);
%                             [txyz(1),txyz(2),txyz(3)] = obj.obs.getCoordinates(S.npsd);
%                             obj.tn(:,end+1) = ttn;
%                             obj.xyz(:,end+1) = txyz;
%                         end
%                         LOOPCOUNT = LOOPCOUNT + 1;
%                         PROG = LOOPCOUNT/TOTALLOOPS;
%                         waitbar(PROG,h,sprintf('%0.2f%% done...',PROG*100));
%                         if(getappdata(h,'canceling'))
%                             delete(h);
%                             error('Characterization aborted');
%                         end
%                     end
%                 end
%             end
%             waitbar(1);
%             catch ex
%             delete(h);
%             rethrow(ex);
%             end
%             delete(h);
%         end
        
        % get PSD generated by the RGB led closest to x,y CIE point
        function [S,R,G,B,tr,tg,tb] = getPSD(obj,x,y)
            P = repmat([x;y;1-x-y],1,size(obj.xyz,2));
            vP = obj.xyz-P;
            dP = sum(vP.*vP,1);
            T = obj.tn(:,find(dP==min(dP),1,'first'));
            tr=T(1);tg=T(2);tb=T(3);
            R = tr*obj.R;G = tg*obj.G;B = tb*obj.B;
            S = R+G+B;
        end
        
    end % methods
    
    methods(Access = private)
        % compute RGB2XYZ transformation matrix
        function T = getRGB2XYZ(obj)
            [Xs,Ys,Zs] = obj.obs.getTristimulusValues([obj.R.npsd,obj.G.npsd,obj.B.npsd]);
            T = [Xs;Ys;Zs];
        end % getRGB2XYZ
        
        % compute RGB2XYZ transformation matrix
        function T = getXYZ2RGB(obj)
            T = zeros(3,3);
            ob1 = cCurve(obj.R.npsd.Xmin,obj.R.npsd.dX,obj.R.npsd.Xmax,obj.R.npsd.Y);
            ob2 = cCurve(obj.G.npsd.Xmin,obj.G.npsd.dX,obj.G.npsd.Xmax,obj.G.npsd.Y);
            ob3 = cCurve(obj.B.npsd.Xmin,obj.B.npsd.dX,obj.B.npsd.Xmax,obj.B.npsd.Y);
            
            v = (obj.obs.ob1).*ob1;
            T(1,1) = v.getIntegral(380,780);
            v = (obj.obs.ob1).*ob2;
            T(2,1) = v.getIntegral(380,780);
            v = (obj.obs.ob1).*ob3;
            T(3,1) = v.getIntegral(380,780);
            
            v = (obj.obs.ob2).*ob1;
            T(1,2) = v.getIntegral(380,780);
            v = (obj.obs.ob2).*ob2;
            T(2,2) = v.getIntegral(380,780);
            v = (obj.obs.ob2).*ob3;
            T(3,2) = v.getIntegral(380,780);
            
            v = (obj.obs.ob3).*ob1;
            T(1,3) = v.getIntegral(380,780);
            v = (obj.obs.ob3).*ob2;
            T(2,3) = v.getIntegral(380,780);
            v = (obj.obs.ob3).*ob3;
            T(3,3) = v.getIntegral(380,780);
        end % getRGB2XYZ
        
    end % methods(Access = private)
end % classdef














