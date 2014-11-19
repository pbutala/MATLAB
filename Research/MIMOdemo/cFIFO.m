classdef cFIFO < handle
% cFIFO: class to handle a FIFO.

    properties(SetAccess = immutable, GetAccess = private)
        BUFSIZE;             % FIFO Buffer Size
    end % properties - immutable
    
    properties(SetAccess = private, GetAccess = private)
        BUFFER;                 % FIFO buffer
        IDXTAIL;                % Push Start Offset
        IDXHEAD;                % Pull Start Offset
    end % properties - private
    
    properties(SetAccess = private, Dependent = true)
        COUNT;                  % Number of queued values in FIFO
    end % properties - private
    
    methods
        % CONSTRUCTOR
        function obj = cFIFO(SZ)
        % cFIFO class constructor
        % obj = cFIFO(SZ) creates an instance of cFIFO with buffer size SZ.
            SZ = SZ+ 1;         % when head==tail, buffer is empty. Thus SZ+1 ensures SZ elements can be queued
            obj.BUFSIZE = SZ;
            obj.BUFFER = zeros(SZ,1);
            obj.IDXTAIL = 1;
            obj.IDXHEAD = 1;
        end % contructor
        
        % ENQUEUE
        function ovFlag = enQ(obj, data)
        % ovFlag = cFIFO.enQ(data)
        % Adds data to top of buffer and returns ovFlag = true if buffer
        % overflows
            CNT = numel(data);
            ovFlag = (CNT > (obj.BUFSIZE - obj.COUNT - 1));
            if ovFlag
                warning('Enqueue created an overflow. Data in buffer might be corrupted.');
            end
            I = obj.getIndexes(obj.IDXTAIL, CNT+1);
            obj.BUFFER(I(1:end-1)) = data;
            obj.IDXTAIL = I(end);
        end % push
        
        % DEQUEUE
        function [data,CNT] = deQ(obj, MAXCOUNT)
        % [data CNT] = cFIFO.deQ(MAXCOUNT)
        % Returns upto 'MAXCOUNT' data values from bottom of buffer.
        % Return value CNT gives number of values returned.
            if MAXCOUNT <= obj.COUNT
                CNT = MAXCOUNT;
            else
                CNT = obj.COUNT;
            end
            if CNT == 0
                warning('No data in queue.');
            end
            I = obj.getIndexes(obj.IDXHEAD, CNT+1);
            data = obj.BUFFER(I(1:end-1));
            obj.IDXHEAD = I(end);
        end % pull
    end % methods
    
    methods(Access = private)
        % GET INDEXES
        function I = getIndexes(obj,IDX,CNT)
        % I = getIndexes(obj,IDX,CNT) returns CNT number of indexes
        % starting at index IDX
            I = IDX : (IDX + CNT - 1);
            while ~isempty(I(I>obj.BUFSIZE))
                I(I>obj.BUFSIZE) = I(I>obj.BUFSIZE) - obj.BUFSIZE;
            end
        end % getIndexes
    end % methods - private
    
    % GETTERS/SETTERS
    methods
        % COUNT
        function val = get.COUNT(obj)
            val = obj.IDXTAIL - obj.IDXHEAD;
            if val < 0
                val = obj.BUFSIZE + val;
            end
        end % COUNT
        
%         % BUFSIZE
%         function val = get.BUFSIZE(obj)
%             % bufsize is set to SZ+1 in constructor to avoid queue empty
%             % conflict. Hence, this getter is implemented to return number
%             % of values that queue can store without overwrite.
%             val = obj.BUFSIZE - 1;
%         end % COUNT
        
    end % methods - getters/setters
end % classdef cFIFO

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    