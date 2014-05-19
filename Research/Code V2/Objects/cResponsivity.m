classdef cResponsivity < handle
% cResponsivity(PDMTRL,PDTYPE)
    properties(SetAccess = immutable)
        Material;
        Type;
    end
    
    properties(SetAccess = private)
        QTMEFFnm;
        QTMEFFval;
    end
    
    methods
        function obj = cResponsivity(varargin)
            narginchk(0,2);
            switch nargin
                case 1
                    MTRL = argin{1};
                case 2
                    MTRL = argin{1};
                    TYPE = argin{2};
            end 

            if exist('MTRL','var')
                if isa(MTRL,'PDMTRL')
                    obj.Material = MTRL;
                else
                    error('MTRL must be of enumeration type ''PDMTRL''');
                end
            else
                obj.Material = PDMTRL.SILICON;
            end
            
            if exist('TYPE','var')
                if isa(TYPE,'PDTYPE')
                    obj.Type = TYPE;
                else
                    error('TYPE must be of enumeration type ''PDTYPE''');
                end
            else
                obj.Type = PDTYPE.PIN;
            end
            obj.Initialize();
        end
        
        function Initialize(obj)
            switch obj.Material
                case PDMTRL.SILICON
                    fl = load('QEFF_SI.csv');
                    obj.QTMEFFnm = fl(:,1);
                    obj.QTMEFFval = fl(:,2);
%                     obj.QTMEFFnm = [400 800 1100];
%                     obj.QTMEFFval = [0.65 0.90 0];
                otherwise
                    error('PDMTRL.%s supported. Others will be incorporated in future versions',char(PDMTRL.SILICON));
            end
            
            switch obj.Type
                case PDTYPE.PIN
                otherwise
                    error('PDMTRL.%s supported. Others will be incorporated in future versions',char(PDTYPE.PIN));
            end
        end
        
        function R = getResponsivity(obj, lambdas)
            QTMEFF = interp1(obj.QTMEFFnm,obj.QTMEFFval,lambdas,'linear','extrap');
            QTMEFF(QTMEFF < 0) = 0;
            R = (lambdas.*QTMEFF)/1240;
        end
    end
end