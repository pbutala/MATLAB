function CCT = mccamy(x,y)
%MCCAMY  Converts CIE coordinates to CCT
%   CCT = 449n3 + 3525n2 + 6823.3n + 5520.33
%   where n = (x - 0.3320)/(0.1858 - y)
%   The calculated CCT becomes less meaningful as the source moves further
%   away from the plankian locus. CCT is only meant to characterize near
%   white lights.
n = (x - 0.3320)/(0.1858 - y);
CCT = 449*(n^3) + 3525*(n^2) + 6823.3*n + 5520.33;

end

