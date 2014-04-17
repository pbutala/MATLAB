% function [frSpOnPx nGdOfSpOnPx] = getPixelvSpot(...)
% This function calculates the fraction of the image of transmitters (spot)
% that lands on each pixel. X Y coordinates must be in the sensor plane.
% Thus, while spot co-ordinates can be matrices where each column
% corresponds to a different receiver location, the pixel co-ordinates are
% just vectors, not matrices.
% 
% The algorithm approximates the values by creating a fine grid along the
% sensor plane and approximating which grid polygons are common to the spot
% and the pixel.
% 
% Inputs: 
% spX1,spY1,spX2,spY2,spX3,spY3,spX4,spY4: XY coordinates of the four
% corners of the spots. 
% pxX1,pxY1,pxX2,pxY2,pxX3,pxY3,pxX4,pxY4: XY coordinates of the four
% corners of the pixels. 
% gdDelX, gdDelY: Grid resolution along the X and Y dimensions
%
% Outputs:
% frSpOnPx: Fraction of the spot that lands on the pixel
% frPxOnSp: Fraction of the signal received by pixel that belongs to a spot
% nGdOfSpOnPx: grid polygon count of spot that lands within pixel

function [frSpOnPx frPxOnSp nGdOfSpOnPx] = getPixelvSpot(spX1,spY1,spX2,spY2,spX3,spY3,spX4,spY4,...
                                 pxX1,pxY1,pxX2,pxY2,pxX3,pxY3,pxX4,pxY4,...
                                 gdDelX, gdDelY)

nGdOfSpOnPx = zeros(numel(spX1),numel(pxX1));   % grid count common to spot and pixel
frSpOnPx = zeros(numel(spX1),numel(pxX1));
spAREATH = realmin('double');

for idSp = 1:numel(spX1)    % for each spot
    % get the minimum and maximum x and y values for spot
    Xs = [spX1(idSp);spX2(idSp);spX3(idSp);spX4(idSp);spX1(idSp)];
    gdXmin = min(Xs);
    gdXmax = max(Xs);
    gdPX = gdXmax-gdXmin;
    resX = gdPX/gdDelX;
    Ys = [spY1(idSp);spY2(idSp);spY3(idSp);spY4(idSp);spY1(idSp)];
    gdYmin = min(Ys);
    gdYmax = max(Ys);
    gdPY = gdYmax-gdYmin;
    resY = gdPY/gdDelY;
    % Create the fine grid based on input resolution
    gdXs = gdXmin+resX/2:resX:gdXmax-resX/2;
    gdYs = gdYmin+resY/2:resY:gdYmax-resY/2;
    [gdxv gdyv] = meshgrid(gdXs,gdYs);
    gdXv = gdxv(:);
    gdYv = gdyv(:);
    clearvars gdXs gdYs gdxv gdyv gdPX gdPY resX resY;
    
    fGdInPx = zeros(numel(gdXv),numel(pxX1));       % flag: 1=grid lies inside pixel
    fGdInSp = zeros(numel(gdXv),1);       % flag: 1=grid lies inside spot
    
    % find grid points which lie inside the spot
    % not all grid points are necessarily in spot based on the shape and
    % orientation of it. grid is a bounding rectangle around it.
    if(polyarea(Xs,Ys)> spAREATH)
        IN = inpolygon(gdXv,gdYv,Xs,Ys);
        fGdInSp = IN;
    end
    
    for idPx = 1:numel(pxX1)    % for each pixel
        % compute the pixel polygon x and y vectors
        xv = [pxX1(idPx);pxX2(idPx);pxX3(idPx);pxX4(idPx);pxX1(idPx)];
        yv = [pxY1(idPx);pxY2(idPx);pxY3(idPx);pxY4(idPx);pxY1(idPx)];
        % find grid points which lie inside pixel
        IN = inpolygon(gdXv,gdYv,xv,yv);
        fGdInPx(:,idPx) = IN;
    end
    
    for idPx = 1:numel(pxX1)    % for each pixel
        % find number of grid points common to pixel and spot
        nGdOfSpOnPx(idSp,idPx) = sum(fGdInSp(find(fGdInPx(:,idPx)==1)),1);
        frSpOnPx(idSp,idPx) = nGdOfSpOnPx(idSp,idPx)/numel(gdXv);
    end
end

frSpOnPx(isnan(frSpOnPx)) = 0;  % 0/0 gives NaN. set those values to 0

% find fraction of signal received by pixel that belongs to a spot
frPxOnSp = frSpOnPx./repmat(sum(frSpOnPx,1),numel(spX1),1);
frPxOnSp(isnan(frPxOnSp)) = 0;  % 0/0 gives NaN. set those values to 0


























