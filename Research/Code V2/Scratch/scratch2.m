close all;
% clearvars;
clc;

for iTsd = 1:LENLEDWID                                                       % LOOP START LED SD
    for iT = 1:LENCCT                                                           % LOOP START CCT
        [x,y] = planckXY(RNGCCT(iT));
        [S,R,G,B,tr,tg,tb] = RGBLED(iTsd).getPSD(x,y);                               % Get PSDs at CCT
        tbB = tb*B; tgG = tg*G;
        BGor = tbB + tgG;
        bg = (tbB.npsd).*(tgG.npsd);
        BGand = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,bg.Y);
%         GRor = tg*G + tr*R;
    end
end