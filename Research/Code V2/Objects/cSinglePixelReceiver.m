classdef cSinglePixelReceiver
    % cSinglePixelReceiver: Abstract Class to handle a single pixel optical receiver
    properties
        location; % location of receiver in global co-ordinates
        orientation = cOrientation(0,0,0); % orientation of receiver
    end
    properties (SetAccess = private)
        optics; % cpc
        sensor; % single pixel sensor
        tia;    % transimpedance amplifier
    end
    properties(Dependent = true, SetAccess = private)
        rxCount; % total number of receivers of this type
        rxFOV;
    end
    
    %% Class methods
    methods
        function obj = cSinglePixelReceiver(o1,o2,o3)
            % cSinglePixelReceiver(Xr,Yr,Zr): Class constructor
            % cSinglePixelReceiver(cLocation): Class constructor
            % -INPUT-
            % Xr,Yr,Zr: X,Y,Z coord of receiver locations. Default (0,0,0)
            if nargin == 1
                obj.location = o1;
            elseif nargin > 1
                obj.location = cLocation(o1,o2,o3);
            else
                obj.location = cLocation(0,0,0);
            end
            obj.optics = cCPC;
            obj.sensor = cSensor(0,0,0);
            obj.tia = cTIA;
        end
    end
    methods(Abstract)
        varargout = getSignal(varargin)
            % varargout = getSignal(varargin)
            % abstract function
    end
    
    %% Get/Set accessors
    methods
        function n = get.rxCount(obj)
            n = obj.location.count;
        end
        
        function fov = get.rxFOV(obj)
            fov = obj.optics.fov;
        end
    end
end