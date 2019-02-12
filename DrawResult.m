function DrawResult(NODES,MEMBERS,MATERIALS,SECTIONS,RELEASES,Qe,d,Scale,NSegs,Color,Width)
NMBs = size(MEMBERS,1);
NON = size(NODES,1);
dXY = reshape(d,3,NON)';
for i=1:NMBs
    NID = MEMBERS(i,1:2);
    MID = MEMBERS(i,3);
    SID = MEMBERS(i,4);
    RID = MEMBERS(i,5:6);
    RLS = RELEASES(RID,:);
    RL2 = reshape(RLS',1,6);
    XY = NODES(NID,:);
    E = MATERIALS(MID,1);%Elasticity
    A = SECTIONS(SID,1);%Section area %%%-------New!!!------------
    I = SECTIONS(SID,2);%Moment of inertial
    Ke = MemberK(XY,E,A,I,[1,1,1;1,1,1]);%Member stiffness matrix. New!!!: A is an input
    Q = Qe(i,:)';
    de = reshape(dXY(NID,:)',6,1);
    T = T2DFrame(XY);%New!!!
    u = T*de;
    uc = find(RL2==1);%Uncondensed degrees of freedoms
    Q(uc,1) = u(uc,1);
    Ke(uc,:)=0;
    for j=1:size(uc,2)
        Ke(uc(j),uc(j))=1;
    end;
    u = Ke\Q;
    DrawMember(XY,Scale*u,NSegs,Color,Width);
end;