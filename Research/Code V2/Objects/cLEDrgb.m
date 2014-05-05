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
%         cct;    % correlated color temperature for each xyz(:,i) combination. using mccamy's formula. More accurate 
    end
    
    methods
        % constructor
        function obj = cLEDrgb(flRes,Rpsd,Gpsd,Bpsd)
             % normalized optical flux resolution (between 0 and 1)
            if ~exist('flRes','var')
                flRes = 0.2;
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
            obj.obs = cCIE;
           % initialize variables
            obj.tn = []; obj.xyz = [];
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
            resa = obj.flRes:obj.flRes:1;
            resLEN = numel(resa);
            IDX = 1;
            s = [0;0;0];
            txyz = [0;0;0];
            h = waitbar(0,'0%','Name','Characterizing LED...',...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
            setappdata(h,'canceling',0);
            for ix = 1:resLEN
                s(1) = resa(ix);
                for iy = 1:resLEN
                    s(2) = resa(iy);
                    for iz = 1:resLEN
                        s(3) = resa(iz);
                        ttn = s/sum(s,1);
                        S = (ttn(1)*obj.R) + (ttn(2)*obj.G) + (ttn(3)*obj.B);
                        [txyz(1),txyz(2),txyz(3)] = obj.obs.getCoordinates(S.npsd);
                        obj.tn(:,end+1) = ttn;
                        obj.xyz(:,end+1) = txyz;
                        PROG = IDX/(resLEN^3);
                        waitbar(PROG,h,sprintf('%0.2f%% done...',PROG*100));
                        IDX = IDX + 1;
                        if(getappdata(h,'canceling'))
                            delete(h);
                            error('Characterization aborted');
                        end
                    end
                end
            end
            waitbar(1);
            obj.tn = unique(obj.tn','rows')';
            obj.xyz = unique(obj.xyz','rows')';
            close(h);
%             NLEN = 8;
%             NSTART = 1;
%             h = waitbar(0,'1','Name','Characterizing LED...',...
%                 'CreateCancelBtn',...
%                 'setappdata(gcbf,''canceling'',1)');
%             setappdata(h,'canceling',0);
%             while NSTART <= resLEN
%                 waitbar(NSTART/resLEN);
%                 NSTOP = min(NSTART+NLEN-1,resLEN);
%                 res = resa(NSTART:NSTOP);
%                 N = numel(res);
%                 s1 = repmat(res,N*N,1);
%                 s2 = repmat(res,N,N);
%                 s3 = repmat(res,1,N*N);
%                 s = [s1(:)';s2(:)';s3(:)'];
%                 sm = repmat(sum(s,1),3,1);
%                 s = s./sm;
% %                 ttn = unique(s','rows')';
%                 ttn = s;
%                 txyz = zeros(3,size(ttn,2));
%                 for j=1:size(ttn,2)
%                     S = (ttn(1,j)*obj.R) + (ttn(2,j)*obj.G) + (ttn(3,j)*obj.B);
%                     [txyz(1,j),txyz(2,j),txyz(3,j)] = obj.obs.getCoordinates(S.npsd);
%                 end
%                 obj.tn(1:size(ttn,1),end+1:end+size(ttn,2)) = ttn;
%                 obj.xyz(1:size(txyz,1),end+1:end+size(txyz,2))  = txyz;
%                 clear s1 s2 s3 s sm ttn txyz;
%                 NSTART = NSTOP + 1;
%                 
%             end
%             waitbar(1);
%             obj.tn = unique(obj.tn','rows')';
%             obj.xyz = unique(obj.xyz','rows')';
%             close(h);
        end
        
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














