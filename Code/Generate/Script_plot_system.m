% Script to plot the system
% change noUsers and noSBS
noUsers = 10;
noSBS   = 6;
flag_plot = 1;

[UE_BS, UEs, BS] = location_voronoi(noUsers, noSBS, flag_plot); % function in ..\voronoi
xlim([-500 500]);
ylim([-500 500]);
% save('pos_BS_UEs.mat', 'UEs', 'BS', "UE_BS")