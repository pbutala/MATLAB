% imPng2Eps(filename)
% converts png file to eps type
function imPng2Eps(filename)
[I,map] = imread(filename,'bmp');
imshow(I,map,'InitialMagnification','fit');
end