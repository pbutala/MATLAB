classdef cTIA
    % Class to handle Transimpedance Amplifier
    % Very basic. May want ot modify to include various other params
    properties(SetAccess = private)
        Inhz = 5e-12; % Input referred noise current density
    end
    
    %% Class Methods
    methods
        function obj = cTIA(Inhz)
            % cTIA(In): class constructor
            % -INPUT-
            % Inhz: Input referred noise current density (A.Hz^{-1/2}). (Default = 5pA/sqrt(Hz))
            if nargin > 0
                obj.Inhz = Inhz;
            else
                obj.Inhz = 5e-12;
            end
        end
        
        function In = getNoise(obj,B)
            % In = getNoise(B)
            % -INPUT-
            % B: Receiver bandwidth
            %
            % -OUTPUT-
            % In: Input referred noise current
            In = (obj.Inhz.^2)*B;
        end
    end
end