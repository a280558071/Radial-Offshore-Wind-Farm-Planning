%% This function define the condition of OWF-ECSP for 30 WTs example

%% Produce the coordinates for WTs
% WTs=60;
% Coord_WT=xlsread('WT2.xlsx','B2:C61');
WTs=42;
Coord_WT=xlsread('HKY_WT.xlsx','B2:C43');

G=graph;
G=addnode(G,WTs);
figure;
p=plot(G,'Layout','force');
p.XData=[Coord_WT(:,1)'];
p.YData=[Coord_WT(:,2)'];

N_WT=1:WTs;

%% Step 1: to find out the cluster centres to be Offshore substations locations
Ns=1; % number of substations
Coord_OS=zeros(Ns,2); % Coordinates of Offshore Substations
% [Coord_OS,U] = fcm(Coord_WT,Ns); % Option1: Use Fuzzy C-means 
% Coord_OS = xlsread('WT2.xlsx','B62:C62'); % Option2: Put the substation node in the centre
Coord_OS = xlsread('HKY_WT.xlsx','B44:C44'); 
N_Subs=(length(N_WT)+1):(length(N_WT)+Ns);
G=addnode(G,length(N_Subs));

%% Step 2: Find all the potential connections between neighbourhood WTs and OSs
% calculate the length based on coordinates
% put them into the Graph G
% MAX RANGE of 60 WTs
% Max_Range= norm([2 3]);% define a max range to limit the number of candidate lines
% I=[];
% J=[];
% len_l=[];
% for i=2:length(N_WT) % Find all the potential connections between neighbourhood WTs
%     for j=1:i-1
%         l=norm(Coord_WT(i,:)-Coord_WT(j,:));
%         if l<=Max_Range
%             G=addedge(G,i,j);
%             I=[I,i];
%             J=[J,j];
%             len_l=[len_l,l];
%         end
%     end
% end
% 
% for i=N_WT % Find all the potential connections between WTs and OSs
%     for j=1:length(N_Subs)
%         l=norm(Coord_WT(i,:)-Coord_OS(j,:));
%         if l<=11 %Max_Range
%             G=addedge(G,i,N_Subs(j));
%             I=[I,i];
%             J=[J,N_Subs(j)];
%             len_l=[len_l,l];
%         end
%     end
% end

% MAX RANGE of 42 WTs
I=[];
J=[];
len_l=[];
for i=2:length(N_WT) % Find all the potential connections between neighbourhood WTs
    for j=1:i-1
        l=norm(Coord_WT(i,:)-Coord_WT(j,:));
        if l<=5
            G=addedge(G,i,j);
            I=[I,i];
            J=[J,j];
            len_l=[len_l,l];
        end
    end
end

for i=N_WT % Find all the potential connections between WTs and OSs
    for j=1:length(N_Subs)
        l=norm(Coord_WT(i,:)-Coord_OS(j,:));
        if l<=5 %Max_Range
            G=addedge(G,i,N_Subs(j));
            I=[I,i];
            J=[J,N_Subs(j)];
            len_l=[len_l,l];
        end
    end
end

% Get some infor about G
In=myincidence(I,J);  % 节支关联矩阵, node-branch incidence matrix, where 
L=length(I); % the total number of candidate lines

%% Step 3: define the Critical Cable Set, i.e. CCS
% now the CCS is defined as the cables connecting offshore substations and
% WTs
CCS=[];
for l=1:L
    if sum(abs(In(N_Subs,l)))~=0
        CCS=[CCS,l];
    end
end

%% Plot the topology of WTs planning scheme
figure;
p=plot(G,'Layout','force');
p.XData=[Coord_WT(:,1)',Coord_OS(:,1)'];
p.YData=[Coord_WT(:,2)',Coord_OS(:,2)'];
hold on;
labelnode(p,N_Subs,{'Sub'});
N_WT=1:WTs;
highlight(p,N_Subs,'Marker','s','NodeColor','c','Markersize',30);
highlight(p,N_WT,'NodeColor','y','Markersize',20);
highlight(p,I,J,'EdgeColor','k','LineStyle','-.','LineWidth',2);
% text(p.XData, p.YData, p.NodeLabel,'HorizontalAlignment', 'center','FontSize', 15); % put nodes' label in right position.
p.NodeLabel={};
hold off

%% Save current line measurements into xlsx file
ConsInf=[[1:length(I)]',I',J',len_l'];
% xlswrite('HKY_lines.xlsx',ConsInf,['F3:I',num2str(3+length(I)-1)]);
n=1;  % power factor of wind turbine
% ConsInf=xlsread('HKY_lines.xlsx',['A1:D',num2str(1+length(I)-1)]);
% Pw=xlsread('HKY_lines.xlsx',['D3:D',num2str(3+WTs-1)]);  % unit: MW
Pw=ones(WTs,1)*10;


%% Base Values
Sbase=1e6;  % unit:VA
% --------------------------------------------------% 66kV--------------------------------------------------
Ubase=66e3;  % unit:V
U=66;
Ibase=Sbase/Ubase/1.732;  %unit: A
Zbase=Ubase/Ibase/1.732;  %unit: Ω
LineCap=[2,3,4,5,6]*12; % Submarine cables: 3×95/185/300/500/800 %42 WTs
C_lines=[185.4,229.2,280.3,355.5,482.7]; %42 WTs
z_b=[0.246+0.135i 0.1328+0.122i 0.0819+0.113i 0.0491+0.105i 0.03+0.094i]; %42 WTs 
% LineCap=[3,4,5,6,7,8]*10; %60 WTs
% C_lines=[180.3,221.8,267.8,306.9,404.2,474.7]; %60 WTs
% z_b=[0.246+0.135i 0.1328+0.122i 0.0819+0.113i 0.0614+0.11i 0.039+0.101i 0.03+0.094i]; %60 WTs
% z=ConsInf(:,4)*(0.028+0.325i)/Zbase; % line impedance unit:p.u., (0.039+0.101i) is the impedance of 3×630 mm2 submarine cables (ohm/km), According to ABB's Cable Catalogue-2016
z=repmat(ConsInf(:,4),1,length(C_lines)).*repmat(z_b,length(I),1)/Zbase; % line impedance unit:p.u., (0.039+0.101i) is the impedance of 3×630 mm2 submarine cables (ohm/km), According to ABB's Cable Catalogue-2016

% z=ConsInf(:,4)*(0.028+0.325i)/Zbase;
% LineCap=86.07; C_lines=482.7;  
% Cost of 66kV line (including construction), unit: 1e4 ￥/km 
% Line Capacity: 86.07 MVA for 3×800 38/66 XLPE cable


% --------------------------------------------------% 35kV--------------------------------------------------
% Ubase=35e3;  % unit:V
% U=35;
% Ibase=Sbase/Ubase/1.732;  %unit: A
% Zbase=Ubase/Ibase/1.732;  %unit: Ω
% LineCap=[1,2,3]*12;
% C_lines=[146.7,207.9,367.4];
% z=ConsInf(:,4)*(0.042+0.348i)/Zbase; % line impedance unit:p.u., (0.039+0.101i) is the impedance of 3×630 mm2 submarine cables (ohm/km), According to ABB's Cable Catalogue-2016
% % LineCap=37.15; % Line Capacity: 37.15 MVA for 3×630 26/35 XLPE cable
% % C_lines=367.4;  % Cost of 35kV line for 3×630 26/35 XLPE cable, unit: 1e4 ￥/km
% % LineCap=53.35; % Line Capacity: 53.35 MVA for 2*3×240 26/35 XLPE cable
% % C_lines=415.8;  % Cost of 35kV line for 2*3×240 26/35 XLPE cable, unit: 1e4 ￥/km


P_Max=ones(length(ConsInf(:,4)),1)*LineCap*1e6/Sbase;   % unit: p.u.

r=real(z);
y_line=1./z;   % The mutual admittance between node i and node j
Bij=abs(imag(y_line));  % The mutual susceptantce between node i and node j
v_min=0.95;
v_max=1.05;
n_feeders=10; % #OSS 出线数
g_Sub_Max=600;
M=1e9;
M2=1e4;
Cost=(1.05*ConsInf(:,4)+0.13*2)*C_lines; % more conservative than engineering
N=WTs+Ns;
L_c=length(CCS);