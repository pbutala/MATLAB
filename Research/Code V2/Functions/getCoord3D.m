%% Usage
% function varargout = getCoord3D(varargin)
% Vr = getVectRot(V,angle,strType);
% Vr = getVectRot(X,Y,Z,angle,strType);
% [Xr Yr Zr] = getVectRot(V,angle,strType);
% [Xr Yr Zr] = getVectRot(X,Y,Z,angle,strType);
% V must have 3 rows, each rows corresponding to a dimention
% Vr is same size as V

%% Function
% function V = getVectRot(P,angle,strType)
function [I J K] = getCoord3D(X,Y,Z,zenxyz,azixyz,tiltxyz)

% if(~isequal(size(Ni),size(Nj),size(Nk),[3 1]))
%     error('Ni Nj Nk must be 3x1 vectors');
% end
%
% if(~isequal(dot(Ni,Nj),dot(Nj,Nk),dot(Nk,Ni),0))
%     error('Ni Nj Nk must be orthogonal');
% end

if(numel(zenxyz)==1)
    zenxyz = zenxyz*ones(size(X));
end

if(numel(azixyz)==1)
    azixyz = azixyz*ones(size(X));
end

if(numel(tiltxyz)==1)
    tiltxyz = tiltxyz*ones(size(X));
end

if(~isequal(size(X),size(Y),size(Z),size(zenxyz),size(azixyz)))
    error('X Y Z must be same size. Elevation and Azimuth can be scalars or must have the same size as X');
end

szP = size(X);
I=zeros(szP);
J=zeros(szP);
K=zeros(szP);
Nxyz = eye(3);
Nel = numel(X);
ip = 1;
while(ip <= Nel)
%     Nijk = getVectRot(getVectRot(Nxyz,zenxyz(ip),'elevation'),azixyz(ip),'azimuth');
    Nijk = getRotation('global',Nxyz,zenxyz(ip),azixyz(ip),tiltxyz(ip),eye(3));
    Vxyz = [X(ip);Y(ip);Z(ip)];
    I(ip) = Nijk(:,1)'*Vxyz;
    J(ip) = Nijk(:,2)'*Vxyz;
    K(ip) = Nijk(:,3)'*Vxyz;
    ip = ip+1;
end

end

%% Space





































