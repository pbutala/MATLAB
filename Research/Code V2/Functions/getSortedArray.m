function A = getSortedArray(A,Mode)
% Y = getSortedArray(A)
% Uses insert-sort -> O(n^2)
% Returns A sorted in 'Ascending' (default mode) or 'Descending' order

N = numel(A);
if ~exist('Mode','var')
    Mode = 'Ascending';
end
switch(lower(Mode))
    case 'ascending'
        for i=2:N
            val = A(i);
            for j=(i-1):-1:1
                if(A(j) > val)          % swap if greater
                    A(j+1) = A(j);
                    A(j) = val;
                else
                    break;
                end
            end
        end
    case 'descending'
        for i=2:N
            val = A(i);
            for j=(i-1):-1:1
                if(A(j) < val)          % swap if smaller
                    A(j+1) = A(j);
                    A(j) = val;
                else
                    break;
                end
            end
        end
    otherwise
        error('Mode should be either ''Ascending'' or ''Descending''');
end

