%% Comment

%% Function
function theta = getVtxAng(V1,V2)

if(~isequal(size(V1),size(V2)))
    error('Size of V1 and V2 must be the same');
end

numVect = numel(V1(1,:));
szV = size(V1);
szTh = szV(2:end);
if numel(szTh) == 1
theta = zeros([1 szTh]);    
else
theta = zeros(szTh);
end
for i=1:numVect
    cosTh = dot(V1(:,i),V2(:,i))/(norm(V1(:,i))*norm(V2(:,i)));
    theta(i) = acos(cosTh);
end

%% Space






























