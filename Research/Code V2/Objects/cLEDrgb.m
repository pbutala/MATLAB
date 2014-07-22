classdef cLEDrgb < handle
    % class to handle RGB LED
    properties(SetAccess = immutable)
        R;      % Red PSD
        G;      % Green PSD
        B;      % Blue PSD
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
        end
        
%         % initialize the class
%         function initialize(obj)
%             res = obj.flRes:obj.flRes:1;
%             N = numel(res);
%             s1 = repmat(res,N*N,1);
%             s2 = repmat(res,N,N);
%             s3 = repmat(res,1,N*N);
%             s = [s1(:)';s2(:)';s3(:)'];
%             sm = repmat(sum(s,1),3,1);
%             s = s./sm;
%             obj.tn = unique(s','rows')';
%             for j=1:size(obj.tn,2)
%                 S = (obj.tn(1,j)*obj.R) + (obj.tn(2,j)*obj.G) + (obj.tn(3,j)*obj.B);
%                 [obj.xyz(1,j),obj.xyz(2,j),obj.xyz(3,j)] = obj.obs.getCoordinates(S.npsd);
%             end
%         end
        % initialize the class
        function initialize(obj)
%             resa = obj.flRes:obj.flRes:1;
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
                            S = (ttn(1)*obj.R) + (ttn(2)*obj.G) + (ttn(3)*obj.B);
                            [txyz(1),txyz(2),txyz(3)] = obj.obs.getCoordinates(S.npsd);
                            obj.tn(:,end+1) = ttn;
                            obj.xyz(:,end+1) = txyz;
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
%             obj.tn = unique(obj.tn','rows')';
%             obj.xyz = unique(obj.xyz','rows')';
            catch ex
            delete(h);
            rethrow(ex);
            end
            delete(h);
        end
%         % initialize the class
%         function initialize(obj)
%             resa = 0:obj.flRes:1;
%             resLEN = numel(resa);
%             s = [0;0;0];
%             txyz = [0;0;0];
%             obj.hWB = waitbar(0,'0.00%','Name','Characterizing LED...',...
%                 'CreateCancelBtn',...
%                 'setappdata(gcbf,''canceling'',1)');
%             try
%             setappdata(obj.hWB,'canceling',0);
%             TOTALLOOPS = (resLEN-1)*resLEN/2;
%             LOOPCOUNT = 0;
%             for ix = 1:resLEN
%                 s(1) = resa(ix);
%                 for iy = 1:resLEN-ix
%                     s(2) = resa(iy);
%                     s(3) = 1-s(1)-s(2);
%                     ttn = s;
%                     S = (ttn(1)*obj.R) + (ttn(2)*obj.G) + (ttn(3)*obj.B);
%                     [txyz(1),txyz(2),txyz(3)] = obj.obs.getCoordinates(S.npsd);
%                     obj.tn(:,end+1) = ttn;
%                     obj.xyz(:,end+1) = txyz;
%                     
%                     LOOPCOUNT = LOOPCOUNT+1;
%                     PROG = LOOPCOUNT/TOTALLOOPS;
%                     waitbar(PROG,obj.hWB,sprintf('%0.2f%% done...',PROG*100));
%                     
%                     if(getappdata(obj.hWB,'canceling'))
%                         delete(obj.hWB);
%                         error('Characterization aborted');
%                     end
%                 end
%             end
% %             waitbar(1);
%             obj.tn = unique(obj.tn','rows')';
%             obj.xyz = unique(obj.xyz','rows')';
%             catch ex
%             delete(obj.hWB);
%             rethrow(ex);
%             end
%             delete(obj.hWB);
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
end % classdef














