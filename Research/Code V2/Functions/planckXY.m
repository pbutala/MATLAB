function [x,y] = planckXY(CCT)
%planckXY  Converts CCT to CIE coordinates
% source: http://en.wikipedia.org/wiki/Planckian_locus#Approximation

T0 = 1667; T1 = 2222; T2 = 4000; T3 = 25000;
C3 = CCT^3;C2 = CCT^2;C1 = CCT;
if (T0 <= CCT) && (CCT <=T1)        % 1667 -> 2222
    x = -(0.2661239*1e9)/(C3) -(0.2343580*1e6)/C2 + (0.8776956*1e3)/C1 + 0.179910;
    y = -1.1063814*(x^3) - 1.34811020*(x^2) + 2.18555832*x - 0.20219683;
elseif (T1 < CCT) && (CCT <=T2)     % 2222 -> 4000
    x = -(0.2661239*1e9)/(C3) -(0.2343580*1e6)/C2 + (0.8776956*1e3)/C1 + 0.179910;
    y = -0.9549476*(x^3) - 1.37418593*(x^2) + 2.09137015*x - 0.16748867;
elseif (T2 < CCT) && (CCT <=T3)     % 4000 -> 25000
    x = -(3.0258469*1e9)/(C3) +(2.1070379*1e6)/C2 + (0.2226347*1e3)/C1 + 0.240390;
    y = 3.0817580*(x^3) - 5.87338670*(x^2) + 3.75112997*x - 0.37001483;
else
    error('CCT must be between %dK and %dK',T0,T3);
end
    
end