function DrawBeam(NODES,MEMBERS,CL,LW)
%CL:Color
%LW: Linewidth
NOM = size(MEMBERS,1);%Number of members
for i=1:NOM
    NID = MEMBERS(i,1:2);
    XY  = NODES(NID,:);
    line(XY(:,1),XY(:,2),'color',CL,'linewidth',LW);
end;