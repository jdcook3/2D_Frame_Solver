clear;
clc;
load Inputs\NODES -ASCII;
load Inputs\MEMBERS -ASCII;
load Inputs\MATERIALS -ASCII;
load Inputs\SECTIONS -ASCII;
load Inputs\MLOADS -ASCII;
load Inputs\NLOADS -ASCII;
load Inputs\BOUNDS -ASCII;
load Inputs\RELEASES -ASCII;
NOM = size(MEMBERS,1);
NON = size(NODES,1);
K = zeros(3*NON,3*NON);%----Three DOFs at each node (New!!!)
for i=1:NOM
    NID = MEMBERS(i,1:2);%Nodal IDs
    MID = MEMBERS(i,3);%Material ID
    SID = MEMBERS(i,4);%Section ID
    RID = MEMBERS(i,5:6);%Release ID
    RLS = RELEASES(RID,:);%Release at each end
    XY = NODES(NID,:);%Nodal coordinates
    E = MATERIALS(MID,1);%Elasticity
    A = SECTIONS(SID,1);%Section area %%%-------New!!!------------
    I = SECTIONS(SID,2);%Moment of inertial
    Ke = MemberK(XY,E,A,I,RLS);%Member stiffness matrix. New!!!: A is an input
    for j=1:2
        for k=1:2
            K(3*NID(j)-2:3*NID(j),3*NID(k)-2:3*NID(k))=K(3*NID(j)-2:3*NID(j),3*NID(k)-2:3*NID(k))+...
                Ke(3*j-2:3*j,3*k-2:3*k);%-----Three DOFs at each node--(New!!!)
        end;
    end;
end;
%-------------Calculate member fixed-end-force----------------------------
NMBLs = size(MLOADS,1);
Pf = zeros(3*NON,1);%New!!!3 DOFs at each node
Qfs = zeros(NMBLs,6);%New!!!6 DOFs for each member
for i=1:NMBLs
    MID = MLOADS(i,1);
    NID = MEMBERS(MID,1:2);
    mID = MEMBERS(MID,3);%Material ID
    SID = MEMBERS(MID,4);%Section ID
    RID = MEMBERS(MID,5:6);%Release ID
    RLS = RELEASES(RID,:);%Release at each end
    E = MATERIALS(mID,1);%Elasticity
    A = SECTIONS(SID,1);%Section area %%%-------New!!!------------
    I = SECTIONS(SID,2);%Moment of inertial
    XY = NODES(NID,:);
    LoadType = MLOADS(i,2);
    Para = MLOADS(i,3:5);
    Qf = FixEndForce(XY,LoadType,Para,E,A,I,RLS);
    T = T2DFrame(XY);%New function!!!
    Qfs(i,:) = (T*Qf)';%Save fixed-end-force in local coordinates.
    Pf(3*NID(1)-2:3*NID(1),1) = Pf(3*NID(1)-2:3*NID(1),1)+Qf(1:3,1);%New!!!3 DOFs at each node
    Pf(3*NID(2)-2:3*NID(2),1) = Pf(3*NID(2)-2:3*NID(2),1)+Qf(4:6,1);%New!!!3 DOFs at each node
end;
%-----Calculate Nodal Force Vector--------------
NNL = size(NLOADS,1);%number of nodal loads
P = zeros(3*NON,1);%New!!!3 DOFs at each node
for i=1:NNL
    NID = NLOADS(i,1);%nodal number
    Loads = NLOADS(i,2:4);%New!!!Revised to do all components at once
    P(3*NID-2:3*NID,1) = P(3*NID-2:3*NID,1)+ Loads';%New!!!Revised to do all components at once
end;
%Calculate Total Force vector
Pt = P-Pf;
%------Assign Boundary Conditions--------------
NOB = size(BOUNDS,1);%Number of Nodal boundary conditions
for i=1:NOB
    NID = BOUNDS(i,1);%Nodal ID
    FixUx = BOUNDS(i,2);%Fixed in X displacement or not?New!!!
    FixUy = BOUNDS(i,3);%Fixed in Y displacement or not?%New!!!: change column
    FixRT = BOUNDS(i,4);%Fixed in Rotation or not?New!!!: change column
    if FixUx==1
        K(3*NID-2,:)=0;
        K(3*NID-2,3*NID-2)=1;
        Pt(3*NID-2,1) = 0;
    end;
    if FixUy ==1
        K(3*NID-1,:)=0;
        K(3*NID-1,3*NID-1)=1;
        Pt(3*NID-1,1) = 0;
    end;
    if FixRT ==1
        K(3*NID,:)=0;
        K(3*NID,3*NID)=1;
        Pt(3*NID,1) = 0;
    end;
end;
%-----Solve for displacement---------------------
d = K\Pt;
%-----Calculate the member force-----------------
Qe = zeros(NOM,6);
for i=1:NOM
    NID = MEMBERS(i,1:2);%Nodal IDs
    MID = MEMBERS(i,3);%Material ID
    SID = MEMBERS(i,4);%Section ID
    RID = MEMBERS(i,5:6);%Release ID
    RLS = RELEASES(RID,:);%Release at each end
    XY = NODES(NID,:);%Nodal coordinates
    E = MATERIALS(MID,1);%Elasticity
    A = SECTIONS(SID,1);%New!!!:SEction area
    I = SECTIONS(SID,2);%Moment of inertial
    Ke = MemberK(XY,E,A,I,RLS);%Member stiffness matrix. New!!:A is an input
    de = zeros(6,1);%New: chang to 6 DOFs for each member
    de(1:3,1) = d(3*NID(1)-2:3*NID(1),1);%New: 3 DOFs for first node
    de(4:6,1) = d(3*NID(2)-2:3*NID(2),1);%New: 3 DOFs for second node
    Qe(i,:) = (Ke*de)';%Global coordinate
    T = T2DFrame(XY);%New function!!!
    Qe(i,:) = T*Qe(i,:)';
end;
for i=1:NMBLs
    MID = MLOADS(i,1);
    Qe(MID,:) = Qe(MID,:)+Qfs(i,:);
end;
%------------------------------------------------
grid on;hold on;
DrawResult(NODES,MEMBERS,MATERIALS,SECTIONS,RELEASES,Qe,d,20,20,'b',1);
DrawResult(NODES,MEMBERS,MATERIALS,SECTIONS,RELEASES,Qe,d,0,20,'k',2);