classdef cFIFO < handle
    properties(SetAccess = immutable)
        BUFSIZEMAX;             % FIFO Buffer Size
    end % properties - immutable
    
    properties(SetAccess = private)
        BUFFER;                 % FIFO buffer
        IDXPUSH = 1;            % Push Start Offset
        IDXPULL = 0;            % Pull Start Offset
    end % properties - private
    
    properties(SetAccess = private, Dependent = true)
        COUNT;                % Number of pending values in FIFO
        isEmpty;                % flag indicating if buffer is empty
    end % properties - private
    
    methods
        function obj = cFIFO(SZ)
        % cFIFO class constructor
        % obj = cFIFO(SZ) creates an instance of cFIFO with buffer size SZ.
            obj.BUFSIZEMAX = SZ;
            obj.BUFFER = zeros(SZ,1);
        end % contructor
        
        function ovFlag = push(obj, data)
        % ovFlag = cFIFO.push(data)
        % Adds data to top of buffer and returns ovFlag = true if buffer
        % overflows
            CNT = numel(data);
            ovFlag = (CNT > (obj.BUFSIZEMAX - obj.COUNT));
            if ovFlag
                warning('Push created an overflow. Data in buffer might be corrupted.');
            end
            I = obj.getIndexes(obj.IDXPUSH, CNT);
            if obj.isEmpty
                obj.IDXPULL = obj.IDXPUSH;
            end
            obj.BUFFER(I(1:end-1)) = data;
            obj.IDXPUSH = I(end);
        end % push
        
        function [data,CNT] = pull(obj, MAXCOUNT)
        % [data CNT] = cFIFO.pull(MAXCOUNT)
        % Returns upto 'MAXCOUNT' data values from bottom of buffer.
        % Return value CNT gives number of values returned.
            if MAXCOUNT <= obj.COUNT
                CNT = MAXCOUNT;
            else
                CNT = obj.COUNT;
            end
            I = obj.getIndexes(obj.IDXPULL, CNT);
            data = obj.BUFFER(I(1:end-1));
            obj.IDXPULL = I(end);
            if obj.IDXPULL == obj.IDXPUSH
                obj.setEmpty();
            end
        end % pull
        
        function [data,CNT] = pullPad(obj, MAXCOUNT, vPad)
        % [data CNT] = cFIFO.pull(MAXCOUNT, vPad)
        % Returns upto 'MAXCOUNT' data values from bottom of buffer.
        % Return value CNT gives number of values returned. 
        % if MAXCOUNT > cFIFO.COUNT, pads returned data with scalar vPad
            [data,CNT] = obj.pull(MAXCOUNT);
            dCNT = MAXCOUNT-CNT;
            if dCNT>0
                if ~exist('vPad','var');
                    vPad = 0;
                end
                data(end+1:end+dCNT) = vPad;
            end
        end % pullPad
    end % methods
    
    methods(Access = private)
        function I = getIndexes(obj,IDX,CNT)
        % I = getIndexes(obj,IDX,CNT) returns CNT number of indexes
        % starting at index IDX
            I = IDX : (IDX + CNT);
            while ~isempty(I(I>obj.BUFSIZEMAX))
                I(I>obj.BUFSIZEMAX) = I(I>obj.BUFSIZEMAX) - obj.BUFSIZEMAX;
            end
        end % getIndexes
        
        function setEmpty(obj)
        % setEmpty() sets the buffer to empty
            obj.IDXPULL = 0;
            % obj.IDXPUSH = 1;
        end
    end % methods - private
    
    % GETTERS/SETTERS
    methods
        function val = get.COUNT(obj)
            if obj.isEmpty
                val = 0;
            else
                val = obj.IDXPUSH - obj.IDXPULL;
                if val <= 0
                    val = obj.BUFSIZEMAX + val;
                end
            end
        end % COUNT
        
        function val = get.isEmpty(obj)
            val = (obj.IDXPULL == 0);
        end % isEmpty
        
    end % methods - getters/setters
end % classdef cFIFO

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    