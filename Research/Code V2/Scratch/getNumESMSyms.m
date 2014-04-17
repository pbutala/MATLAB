function Cnt = getNumESMSyms(Nt,M)
    Cnt = 0;
    for k = 0:Nt
        Cnt = Cnt + nchoosek(Nt,k)*power(M-1,k);
    end
end