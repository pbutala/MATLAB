function STR = getTimeString(TIME)
    SECINMIN = 60; MININHR = 60; HRINDAY = 24; DAYINYR = 365;
    if TIME < SECINMIN
        STR = sprintf('%d sec',ceil(TIME));
    elseif TIME < SECINMIN*MININHR
        STR = sprintf('%d mins',ceil(TIME/SECINMIN));
    elseif TIME < SECINMIN*MININHR*HRINDAY
        STR = sprintf('%d hr',ceil(TIME/(SECINMIN*MININHR)));
    elseif TIME < SECINMIN*MININHR*HRINDAY*DAYINYR
        STR = sprintf('%d days',ceil(TIME/(SECINMIN*MININHR*HRINDAY)));
    else
        STR = sprintf('%d yr',ceil(TIME/(SECINMIN*MININHR*HRINDAY*DAYINYR)));
    end
end