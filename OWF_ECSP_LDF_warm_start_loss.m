    %% Offshore wind farm collector system planning with radial topology/without "N-1" criterion based on LDF
% This function solve an OWF-ECS planning problem, with wind turbines (WTs)
% turbines (WTs) location fixed and 1 substation located in their centre
% Radial topology of OWF farm is considered.
% Power flow model is based on Linear DistFlow
%  For basic notations, see more at:
%  Shen, Xinwei, S. Li, and H. Li. "Large-scale Offshore Wind Farm Electrical Collector System Planning: A Mixed-Integer Linear Programming Approach." arXiv preprint arXiv:2108.08569 (2021).

%% Variable statement
L_c=1;

n_cab=size(LineCap,2);
x=binvar(L,n_cab,'full');     %Vars for line construction, x(i,1)==1 dentoes that line i is constructed.
y=binvar(L,n_cab,'full');   %Vars for line operation flag in different contigencies, y(line operation,Cont_l)==0
y_ij=binvar(L,n_cab,2,'full');  %Vars for parent-child relationship of line l (i→j)
Pij=sdpvar(L,1,'full');    %Vars for active power flow in each line, unit: p.u. 
Pij_all=sdpvar(L,n_cab,'full'); %Vars for active power flow in different candidate line, unit: p.u. 
Pw_shed=sdpvar(WTs,1, 'full'); % shedded active wind power 
g_Sub_P=sdpvar(Ns,1,'full');    %Vars for generated power of Subs
% topology
F=sdpvar(L,1,'full'); % Fictitious flow in each line, to complete SCF constraints
D=ones(N,1); % Fictitious demand in each load
be=binvar(L,2,'full'); % Vars for line direction flag in ST cons
s=ConsInf(:,2);
t=ConsInf(:,3);
%% ***********Objective Function***********
Obj_inv=sum(sum(Cost.*x))*1e4;
Ur_Hours=2000/8760;  % annual hours of OWFs
Pr_ele=200;  % 0.85 ￥/kWh = 850 ￥/MWh
Obj_loss=Ur_Hours*8760*20*sum(sum(r.*(Pij_all.^2)))*Pr_ele;
Obj_WindCurt=M*sum(sum(Pw_shed));
Obj=Obj_inv+Obj_WindCurt+Obj_loss;
% Obj=Obj_inv+Obj_WindCurt; %HKY
% Obj=sum(Cost.*x);

%% *********** Constraints ***********
Cons=[];

%% Cons1: Operation logic y<=x in any contigency C_l
Cons_Op=[];
Cons_Op=[Cons_Op, y<=x];
Cons_Op=[Cons_Op, sum(sum(y,2),1)==N-length(N_Subs)];
Cons_Op=[Cons_Op, y==sum(y_ij,3)];

Cons=[Cons,Cons_Op];
display('****** Cons. on Construction Logic Completed!******');
size(Cons_Op);
size(Cons);
%% Cons2: Crossing-avoidance constraints (CAC)
Cons_CAC=[];
Cr_Cab=Find_Cr_Cab(In,L,N_WT,N_Subs,I,J,Coord_WT,Coord_OS);
for l=1:length(Cr_Cab)
    Cons_CAC=[Cons_CAC,sum(x(Cr_Cab(l,1),:),2)+sum(x(Cr_Cab(l,2),:))<=1];
end
Cons=[Cons,Cons_CAC];
display('****** Crossing-avoidance Cons.(CAC) Completed!******');


%% Cons3: Spanning Tree Constr.
Cons_ST=[];
for l=1:L
    Os=In(s(l),:);
    Os1=find(Os==-1); % the set of lines end at node s(l)
    Os2=find(Os==1); % the set of lines start from node s(l)
    Ot=In(t(l),:);
    Ot1=find(Ot==-1); % the set of lines end at node t(l)
    Ot2=find(Ot==1); % the set of lines start from node t(l)
    % Cons_ST=[Cons_ST,be(l,1)+be(l,2)==sum(zij(l,:),2)];
    Cons_ST=[Cons_ST,be(l,1)+be(l,2)==x(l,1)];
    if ismember(s(l),Ns)
        Cons_ST=[Cons_ST,be(l,2)==0]; % be(l,2)==0 means node t(l) cannot be the parent node of s(l)
    end
    if ismember(t(l),Ns)
        Cons_ST=[Cons_ST,be(l,1)==0]; % be(l,1)==0 means node s(l) cannot be the parent node of t(l)
    end
    if ismember(s(l),WTs)
        Cons_ST=[Cons_ST,sum([be(Os1,1);be(Os2,2)])==1]; % sum(be(Os,2))==1 means only one node in Os can be the parent node of s(l)
    end
    if ismember(t(l),WTs)
        Cons_ST=[Cons_ST,sum([be(Ot1,1);be(Ot2,2)])==1]; % sum(be(Ot,1)==1 means only one node in Ot can be the parent node of t(l)
    end
end
Cons=[Cons,Cons_ST];
display('******Cons. on Spanning Tree Completed!******');
size(Cons_ST);
size(Cons);

% % Cons: Single Commodity Flow Constr.
% Cons_SCF=[];
% for k=1:length(WTs)
%     li=In(k,:);
%     lij=find(li==1); % the set of lines start from node i
%     lki=find(li==-1); % the set of lines end at node i
%     Cons_SCF=[Cons_SCF,sum(F(lij,:))-D(k,:)==sum(F(lki,:))];
%     % Cons_SCF=[Cons_SCF,abs(F(k,:))<=sum(zij(k,:),2)*M]
%     Cons_SCF=[Cons_SCF,abs(F(k,:))<=x(k,:)*M];
% end
% Cons=[Cons,Cons_SCF];
%% Cons4: Power balance
Cons_S=[];
Cons_S=[Cons_S,In*Pij(:,1)==[-(Pw-Pw_shed(:,1));g_Sub_P(:,1)]];
Cons_S=[Cons_S,Pij==sum(Pij_all,2)];
Cons_S=[Cons_S,Pw>=Pw_shed(:,1)>=0];
Cons=[Cons,Cons_S];
display('******Cons. on Power Balance Completed!******')
size(Cons_S);
size(Cons);
%% Cons5: Power flow limitation in each line
Cons_Line=[];
% for C_l=1:L_c % for Contigency C_l and Normal condition, the power limits must maintain
%     Cons_Line=[Cons_Line,-y(:,C_l).*LineCap<=Pij(:,C_l)<=y(:,C_l).*LineCap];
% end

for i=1:L
    for j=1:n_cab
        Cons_Line=[Cons_Line,-y(i,j)*LineCap(j)<=Pij_all(i,j)<=y(i,j)*LineCap(j)];
    end
end
Cons=[Cons,Cons_Line];
display('******Cons. on Power Limits of Lines Completed!******')
size(Cons_Line);
size(Cons);

%% Cons6: Power limits of Subs
Cons_Sub=[0<=g_Sub_P<=g_Sub_Max];
Cab_Sub=find(ConsInf(:,3)==N_Subs);
Cons_Sub=[Cons_Sub,sum(sum(x(Cab_Sub,:),2),1)<=n_feeders];
Cons=[Cons,Cons_Sub];
display('******Cons. on Power Limits of Subs Completed!******')
size(Cons_Sub);
size(Cons);

%% Solving Options: B&B or B&C, Heuristics percentage
% load ('60WT_1Sub_LDF_66kV_1218_1435No_Loss&MST.mat','s_x');
% load ('42WT_1Sub_LDF_66kV_1220_1150_No_loss&ST.mat','s_x');
% load ('60WT_1Sub_LDF_66kV_1220_2243_WarmStart_8feeders.mat','s_x');
% load ('HKY_60WT_8feeders.mat','s_x');
load('60WT_1Sub_LDF_66kV_0109_1557No_Loss&MST_8feeders.mat','s_x');
% load('42WT_1Sub_LDF_66kV_0110_1550No_Loss&MST_10feeders.mat','s_x');
assign(x,s_x);
% ops=sdpsettings('solver','COPT');
% ops=sdpsettings('solver','cplex','verbose',2,'usex0',0,'Cplex.Benders.Strategy',0); %,'gurobi.MIPGap',5e-2,,'Gurobi.MIPFocus',3,'gurobi.MIPGap',1e-3,'gurobi.TimeLimit',30000,
% ops=sdpsettings('solver','gurobi','usex0',1, 'gurobi.MIPGap',1e-4,'verbose',2);%,'gurobi.Cuts',0,'usex0',1,
ops=sdpsettings('solver', 'gurobi', 'usex0', 1, 'gurobi.MIPGap', 1e-4, 'verbose', 2); % , 'gurobi.MIPFocus', 1
% [model,recoverymodel] = export(Cons,Obj,ops);
% ops.cplex.exportmodel='abc.lp';
%% solve the problem
sol1=optimize(Cons,Obj,ops);
%% Save the solution with "s_" as start
s_x=round(value(x));
s_y=round(value(y));
s_y_ij=round(value(y_ij));
s_Pij=value(Pij);
s_Pij_all=value(Pij_all);
s_Pw_shed=value(Pw_shed);
s_g_Sub=value(g_Sub_P);
s_Obj=value(Obj);
s_Obj_inv=value(Obj_inv);
s_Obj_WindCurt=value(Obj_WindCurt);
s_Obj_loss=value(Obj_loss);

%% display the results' infor
display(['********* 风电场集电系统规划结束！']);
display(['********* 最优拓扑结构如图所示！']);
display(['********* 1 线路规划建设成本：￥ ',num2str(s_Obj_inv)]);
display(['********* 2 弃风成本：￥ ',num2str(s_Obj_WindCurt)]);
display(['********* 3 集电系统网损成本：￥ ',num2str(s_Obj_loss)]);
Ploss=sum(sum(r.*(s_Pij.^2)));
Ploss_rate=Ploss/sum(Pw);
C_EENG=abs(Pr_ele*2000*20*0.0045*(len_l*(sum(s_x,2).*s_Pij)));
display(['********* 4 集电系统可靠性成本：￥ ',num2str(C_EENG)]);
display(['********* （集电系统网损率为： ',num2str(Ploss_rate*100), ' %']);
display(['********* 总成本：￥ ',num2str(s_Obj+C_EENG)]);
