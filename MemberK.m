function Ke = MemberK(XY,E,A,I,RLS)
RLS2 = reshape(RLS',1,6);%Put all in one row
L = sqrt((XY(2,1)-XY(1,1))^2+(XY(2,2)-XY(1,2))^2);
EI12 = 12*E*I/L^3;
EI6 = 6*E*I/L^2;
EI4 = 4*E*I/L;
EI2 = 2*E*I/L;
EA  = E*A/L;
%----------New matrix for frame element--------------------
K =    [EA      0       0       -EA     0       0
        0       EI12    EI6     0       -EI12   EI6
        0       EI6     EI4     0       -EI6    EI2
        -EA     0       0       EA      0       0
        0       -EI12   -EI6    0       EI12    -EI6
        0       EI6     EI2     0       -EI6    EI4];
%----------Need transformation matrix for frame element---------
C = find(RLS2==0);
if size(C,2)>0
    U = find(RLS2==1);
    Kuu = K(U,U);
    Kuc = K(U,C);
    Kcu = K(C,U);
    Kcc = K(C,C);
    Kcd = Kuu-Kuc*(Kcc\Kcu);
    K = zeros(6,6);
    K(U,U)=Kcd;
end;
T = T2DFrame(XY);%New function
%-----------Transform to global coordinate----------------------
Ke = T'*K*T;
