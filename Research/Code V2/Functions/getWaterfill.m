% basic waterfilling script

function [u,Pp,Pj] = getWaterfill(Pavg, Noise, Gains, tol)
narginchk(3,4);
if nargin == 3
    tol = Pavg/1000;
end


G2 = Gains.*Gains;
b = Noise./G2;
Nc = numel(find(Gains > 0)); % rank of H

% Gavg = mean(Gains(:));
% Gavg2 = Gavg^2;

u = Pavg/Nc;

flag = true;
while(flag)
    Pj = u-b;
    Pj(Pj<0) = 0;
    
    Pp = sum(Pj(:))/Nc;
    dP = Pp - Pavg;
    
    if (dP > 0) || (dP < -tol)
        du = dP/Nc;
        u = u - du;
    else
        flag = false;
    end
end
    
end