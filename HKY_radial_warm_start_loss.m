% HKY_radial
%% Offshore wind farm collector system planning with radial topology/without "N-1" criterion based on LDF
% This function solve aN OWF-ECS planning problem, with 30 wind turbines (WTs)
% turbines (WTs) location fixed and 1 substation located in their centre
% Radial topology of OWF farm is considered.
% Power flow model is based on Linear DistFlow
%  For basic notations, see more at:
%  Shen, Xinwei, S. Li, and H. Li. "Large-scale Offshore Wind Farm Electrical Collector System Planning: A Mixed-Integer Linear Programming Approach." arXiv preprint arXiv:2108.08569 (2021).
clear all
close all

%% Conditions initialization for 30 WTs cases
% HKY30_conditions_warm_start_loss_60WT;
HKY30_conditions_warm_start_loss_42WT;

%% Model formulation and solve the problem
% OWF_ECSP_LDF_warm_start_loss;
OWF_ECSP_LDF_warm_start_loss_noCurtail;

%% save the data
% 获取当前时间
currentTime = datestr(now, 'mmdd_HHMM');
% 构造文件名
filename_1 = [num2str(WTs),'WT_',num2str(Ns),'Sub_LDF_',num2str(U), 'kV_',currentTime,'_ST_WarmStart_loss_10feeders'];
% filename_1 = [num2str(WTs),'WT_',num2str(Ns),'Sub_LDF_',num2str(U), 'kV_',currentTime,'No_Loss&MST_10feeders'];
% 保存文件
save(filename_1);
%% Highlight the lines to be bulit and plot all the operation conditions
% pi_final=plot_ECSP_Colored(2,L,N_WT,N_Subs,I,J,s_x,s_y,s_y_ij,s_Pij,Coord_WT,Coord_OS,[],Inp,Inn,[]);
% pi_final=plot_ECSP_DCPF(1,L,N_WT,N_Subs,I,J,s_x,s_y,s_Pij,[],Coord_WT,Coord_OS,[],n_cab);
pi_final=plot_purple(1,L,N_WT,N_Subs,I,J,s_x,s_y,s_Pij,[],Coord_WT,Coord_OS,[],n_cab);
% print(3,'-dpng',[num2str(WTs),'WT_',num2str(Ns),'Sub_LDF_35kV_3_real_New_24']);
% print(1,'-dpng',[num2str(WTs),'WT_',num2str(Ns),'Sub_LDF_66kV_5_real_New_8']);
for figNum = 1:4
    filename_2 = [num2str(WTs),'WT_',num2str(Ns),'Sub_LDF',num2str(U), 'kV_',currentTime, '_Fig', num2str(figNum)];
    print(figNum, '-dpng', filename_2);
end
cab_len=len_l*s_x; 
%% 3. 修改ops调用BD
