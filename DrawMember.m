function DrawMember(XY,u,NSegs,Color,Width)
%I rewrite the whole function to match with frame element (need
%transformation)
L = sqrt((XY(2,1)-XY(1,1))^2+(XY(2,2)-XY(1,2))^2);
X = zeros(NSegs+1,1);%Global coordinate
Y = zeros(NSegs+1,1);%Global coordinate
x = zeros(NSegs+1,1);%x in local coordinate
dX = (XY(2,1)-XY(1,1))/NSegs;
dY = (XY(2,2)-XY(1,2))/NSegs;
dL = L/NSegs;
for i=1:(NSegs+1)
    X(i) = XY(1,1)+(i-1)*dX;
    Y(i) = XY(1,2)+(i-1)*dY;
    x(i) = (i-1)*dL;
end;
N1 = 2*x.^3/L^3-3*x.^2/L^2+1;%all in local coordinate
N2 = x.*(1-x/L).^2;%all in local coordinate
N3 = 3*(x/L).^2-2*(x/L).^3;%all in local coordinate
N4 = (x.^2/L).*(-1+x/L);%all in local coordinate
NA1 = 1-x/L;%Linear shape function in axial direction (see truss section)
NA2 = x/L;%Linear shape function in axial direction (see truss section)
v = N1*u(2)+N2*u(3)+N3*u(5)+N4*u(6);
u = NA1*u(1)+NA2*u(4);
C = (XY(2,1)-XY(1,1))/L;
S = (XY(2,2)-XY(1,2))/L;
T = [C,S;-S,C];%Need to transform displacement into global coordinates
dXY = T'*[u';v'];
Xd = X+dXY(1,:)';%Displaced X coordinate of each point
Yd = Y+dXY(2,:)';%Displaced Y coordinate of each point
line(Xd,Yd,'color',Color,'linewidth',Width);