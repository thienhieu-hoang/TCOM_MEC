% Function to initial locations of UEs and SBSs
% Input:
%   noUsers   == double == number of active UEs
%   noSBS     == double == number of SBSs
%   flag_plot == bollean, = 1 to plot to test, = 0 to skip
function [UE_BS, UEs, BS] = location_voronoi(noUsers, noSBS, flag_plot)
% Output:
% UEs == 1x1 struct
%       UEs.active   == N_active x 2 matrix
%                               (1st col == x-coordinate
%                                2nd col == y-coordinate)
%       UEs.inactive == N_inactive x 2 matrix
%                               (1st col == x-coordinate
%                                2nd col == y-coordinate)
%       UEs.inBS     == 1 x N_active  : SBS that covers the active UEs
%                              example: UEs.inBS(2) = 4 means...
%                                         ...UE 2 in coverage of SBS 4
% BS  == 1x1 struct
%       BS.positions == N_sbs x 2 matrix
%                               (1st col == x-coordinate
%                                2nd col == y-coordinate)
%       BS.SBS       == N_sbs x 1 cell : save the positions of UEs that the SBS covers
%                       example: BS.SBS{1} == [150 100;
%                                              120 200;
%                                             -125 100]
%                                     --> SBS1 covers the UEs at (150,100),
%                                                     (120,200), (-125,100)
% UE_BS  == N x M matrix % matrix of relation of UEs and SBSs

% Initialization
x0 = 0;
y0 = 0;
mbs_center = [0,0];          % the origin is the center of the simulation region

% Square Region 2D
Sqr.frameSize = [1000, 1000];        % [m] width - height
Sqr.area = Sqr.frameSize(1) * Sqr.frameSize(2);     % [m^2]
%%

lambda_active_ue = 10*1e-6;             % [point/m] mean density
intensity_active_ue = lambda_active_ue*Sqr.area;    % average number of point inside circular region

N_active_ue = noUsers;
UEs.active(:,1) = x0 + Sqr.frameSize(1)/2*(-1+2*rand(N_active_ue,1));  % [m] x-coordinate of points
UEs.active(:,2) = y0 + Sqr.frameSize(2)/2*(-1+2*rand(N_active_ue,1));  % [m] y-coordinate of points

%%

lambda_inactive_ue = 5*1e-6;             % [point/m] mean density
% random setting
% intensity_inactive_ue = poissrnd(lambda_inactive_ue*Sqr.area)
% fixed setting
intensity_inactive_ue = 5;

N_inactive_ue = 2;

UEs.inactive(:,1) = x0 + Sqr.frameSize(1)/2*(-1+2*rand(N_inactive_ue,1));  % [m] x-coordinate of points
UEs.inactive(:,2) = y0 + Sqr.frameSize(2)/2*(-1+2*rand(N_inactive_ue,1));  % [m] y-coordinate of points

%%

lambda_sbs = 10*1e-6;             % [point/m] mean density
intensity_sbs = lambda_sbs*Sqr.area;

N_sbs = noSBS;

% Locations of SBSs
BS.positions(:,1) = x0 + Sqr.frameSize(1)/2*(-1+2*rand(N_sbs,1));  % [m] x-coordinate of points
BS.positions(:,2) = y0 + Sqr.frameSize(2)/2*(-1+2*rand(N_sbs,1));  % [m] y-coordinate of points


%% Plotting Location of Points Inside 2D Circular Region
if flag_plot
    figure;

    plot(mbs_center(1), mbs_center(2),'md','MarkerFaceColor','m', 'HandleVisibility','off'); hold on;                % Location of points
    text(mbs_center(:,1)+13, mbs_center(:,2), 'MBS', ...
        'HorizontalAlignment','left')

    plot(UEs.active(:,1), UEs.active(:,2), 'bo','MarkerFaceColor','b'); hold on;                % Location of points

    plot(UEs.inactive(:,1), UEs.inactive(:,2),'ro','MarkerFaceColor','r'); hold on;                % Location of points

    plot(BS.positions(:,1), BS.positions(:,2),'gs','MarkerFaceColor','g'); hold on;                % Location of points

end
%% Voronoi network

network_size = 2*Sqr.frameSize(1);
num_cells = N_sbs;
num_users = N_active_ue;

% Plot voronoi area of cells
if flag_plot
    voronoi(BS.positions(:,1), BS.positions(:,2),'b'); hold on

    % Assign labels to the points.
    nump = size(BS.positions,1);
    plabels = arrayfun(@(n) {sprintf('SBS%d', n)}, (1:nump)');
    hold on
    Hpl = text(BS.positions(:,1)+13, BS.positions(:,2), plabels, ...
        'HorizontalAlignment','left');
end
%% Extract information of the network
% 1- Positions of cluster heads
% 2- How many nodes in each clusters (which is automatically assigned by
% voronoi function of Matlab)?
% 3- Positions of nodes in each cluster

cluster_info = cell(num_cells,1);
BS.SBS = cell(num_cells,1);
UE_BS  = zeros(size(UEs.active,1), size(BS.positions,1));
% matrix of relation of UEs and SBSs

for kk=1:num_users
    min_distance = network_size;
    cluster_id = 0;
    for uu=1:num_cells
        % distance between each user and its cluster head
        dist_user_ch = sqrt((UEs.active(kk,1)-BS.positions(uu,1))^2 ...
            + (UEs.active(kk,2)-BS.positions(uu,2))^2);
        if dist_user_ch < min_distance
            min_distance = dist_user_ch;
            cluster_id = uu;
        end
    end
    if cluster_id > 0
        BS.SBS{cluster_id} = [BS.SBS{cluster_id}; UEs.active(kk,1) UEs.active(kk,2) kk];
        UEs.inBS(kk) = cluster_id;
        UE_BS(kk,cluster_id) = 1;
    end
end

%% Save positions of nodes in each cluster to separated files

% test
if flag_plot
    for jj = 1:num_cells
        if (~isempty(BS.SBS{jj}))
            text(BS.SBS{jj,1}(:,1)+13, ...
                BS.SBS{jj,1}(:,2), num2str(jj), ...
                'HorizontalAlignment','left'); hold on

            plot(BS.SBS{jj}(:,1), ...
                BS.SBS{jj}(:,2), ...
                '.', 'Color', 'None' ,'MarkerSize',15,'MarkerEdgeColor','none', 'HandleVisibility','off'); hold on
        end
    end

    legend('active UEs', 'inactive UEs', 'SBS');
    xlabel('$x$ [m]','Interpreter','LaTex');
    ylabel('$y$ [m]','Interpreter','LaTex');
    % save('pos_BS_UEs.mat', 'UEs', 'BS', "UE_BS")
end

end

