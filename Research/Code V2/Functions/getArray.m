function A = getArray(N,p)
    % A = getArray(N,p)
    % gets centered array about 0.
    %
    % -INPUT-
    % N: number of elements, 
    % p: pitch
    %
    % -OUTPUT-
    % A: centered array about 0.
    
    D = (N-1)*p;
    A = -D/2:p:D/2;
    if isempty(A)
        A = zeros(1,N);
    end
end