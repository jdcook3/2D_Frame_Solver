function Qf = FixEndForce(XY,LoadType,Para,E,A,I,RLS)
L = sqrt((XY(2,1)-XY(1,1))^2+(XY(2,2)-XY(1,2))^2);
Qf = zeros(4,1);
if LoadType==1
    l1 = Para(1);
    l2 = Para(2);
    W  = Para(3);
    FSb = -W*l2^2*(3*l1+l2)/L^3;
    FMb = -W*l1*l2^2/L^2;
    FSe = -W*l1^2*(l1+3*l2)/L^3;
    FMe = W*l1^2*l2/L^2;
    Qf = [0;FSb;FMb;0;FSe;FMe];%New!!!:Add axial components
elseif LoadType==2
elseif LoadType==3
    l1 = Para(1);
    l2 = Para(2);
    w  = Para(3);
    FSb = (-w*L/2)*(1-(l1/L^4)*(2*L^3-2*l1^2*L+l1^3)-(l2^3/L^4)*(2*L-l2));
    FMb = (-w*L^2/12)*(1-(l1^2/L^4)*(6*L^2-8*l1*L+3*l1^2)-(l2^3/L^4)*(4*L-3*l2));
    FSe = (-w*L/2)*(1-(l1^3/L^4)*(2*L-l1)-(l2/L^4)*(2*L^3-2*l2^2*L+l2^3));
    FMe = (w*L^2/12)*(1-(l1^3/L^4)*(4*L-3*l1)-(l2^2/L^4)*(6*L^2-8*l2*L+3*l2^2));
    Qf = [0;FSb;FMb;0;FSe;FMe];%New!!!Add axial components
elseif LoadType==4%New!!:We use trangle only, so w1=0 and w2 = w
    l1 = Para(1);
    l2 = Para(2);
    w  = Para(3);
    FSb = ((-w*(L-l1)^3)/(20*L^3))*((3*L+2*l1)*(1+l2/(L-l1)+l2^2/(L-l1)^2)-(l2^3/(L-l1)^2)*(2+(15*L-8*l2)/(L-l1)));
    FMb = (-w*(L-l1)^3/(60*L^2))*((2*L+3*l1)*(1+l2/(L-l1)+l2^2/(L-l1)^2)-(3*l2^3/(L-l1)^2*(1+(5*L-4*l2)/(L-l1))));
    FSe = (-w/2)*(L-l1-l2)-FSb;
    FMe = ((L-l1-l2)/6)*(w*(L-l1+2*l2))+FSb*L-FMb;
    Qf = [0;FSb;FMb;0;FSe;FMe];%New!!!
elseif LoadType==5%New!!!
    l1 = Para(1);
    l2 = Para(2);
    w  = Para(3);
    FAb = -w*l2/L;
    FAe = -w*l1/L;
    Qf = [FAb;0;0;FAe;0;0];
elseif LoadType==6
elseif LoadType==7
end;
%------------------Fix-end force with release-----------------------------
RLS2 = reshape(RLS',1,6);%Put all in one row
C = find(RLS2==0);
if size(C,2)>0
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
    U = find(RLS2==1);
    Kuc = K(U,C);
    Kcc = K(C,C);
    Qfu = Qf(U,1);%Partioned fixed-end force vector.
    Qfc = Qf(C,1);%Partioned fixed-end foce vector.
    Qfcd = Qfu-Kuc*(Kcc\Qfc);% Qf vector in condensed form by equation derived in class.
    Qf = zeros(6,1);
    Qf(U,1)=Qfcd;
end;
%----------New!!!Need transformation matrix for frame element---------
T = T2DFrame(XY);%New function
%-----------%New!!!Transform to global coordinate----------------------
Qf = T'*Qf;