%% function SOG = getSOG(mn, sd, scl, X)
% function calculates and returns 'Sum of Gaussians' (SOG)
% mean: matrix containing mean for each gaussian distribution
% sd: matrix containing standard deviation for each gaussian distribution
% scl: matrix containing scale factor for each gaussian distribution

function SOG = getSOG(mean, sd, scl, X)
if(~isequal(size(mean),size(mean),size(mean)))
    error('''Mean'', ''Variance'' and ''Scale'' should be of same size');
end
SOG = 0;
for idx = 1:numel(mean)
    if isequal(sd(idx),0)
        y = zeros(size(X));
        y(X == mean(idx)) = 1;
        SOG = SOG + scl(idx)*y;
    else
        SOG = SOG + (scl(idx)/(sqrt(2*pi)*sd(idx)))*exp(-((X-mean(idx)).^2)/(2*sd(idx)*sd(idx)));
    end
end