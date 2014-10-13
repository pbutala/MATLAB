function BITS = getGrayBits(M)
    bps = log2(M);
    BITS = zeros(M,bps);
    if bps>0
        BITS(1,end) = 0;
        BITS(2,end) = 1;
    end
    for i=1:bps-1
        R = power(2,i);
        BITS(R+1:2*R,end-i+1:end) = BITS(R:-1:1,end-i+1:end);
        BITS(R+1:2*R,end-i) = 1;
    end
end