%------------------------------
% compare proposed algorithm with exhaustive search in term of convergence behavior and runtime
% ----------------------------- 
clear
clear all

addpath('..\WOA_voronoi\')
addpath('..\Generate\')
addpath('..\')
tic
load('..\Parameters\parameters2.mat')

rng('default')

noSearchAgents = 30;
maxIter = 1500;
NoUsers = 2; % values of N
noBSs   = 4;
noAnten = 4;
noSubcs = 5;

noRealizations = 1;

doTol = 1;

% po: percentage offloading
% su: system utility

po_MECNOMA21   = zeros(length(NoUsers), noRealizations);   % 5 x noReal matrix
su_MECNOMA21   = zeros(length(NoUsers), noRealizations);   % 5 x noReal matrix
time_MECNOMA21 = zeros(length(NoUsers), noRealizations);   % 5 x noReal matrix

po_EX   = zeros(length(NoUsers), noRealizations);   % 5 x noReal matrix
su_EX   = zeros(length(NoUsers), noRealizations);   % 5 x noReal matrix
time_EX = zeros(length(NoUsers), noRealizations);   % 5 x noReal matrix

for iN = 1:length(NoUsers)

    users_no = NoUsers(iN);
    name = sprintf('../Conver_behave/position_data/pos_BS_UEs_%dUE.mat', users_no); 
    load(name);

    xt(iN) = {num2str(users_no)};
    H2H = cell(noRealizations,1); 
    % channelGain == noRealizations x 1 cell
    %                  each cell is a N x N x M x K matrix (h2h), 
    %                   ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
    %                       h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k| 
    R_nm   = zeros(users_no, noBSs, noRealizations);
    for iReal = 1:noRealizations
        [h2h, hA, dA] = channelMod(UEs, BS, noAnten, noSubcs, logNormalMean, logNormalDeviation);
            % function in ..\
            % return hA  == N x M x K cell, 
            %               each cell is a L x 1 vector  == vector of channel gain
            %                (each SBS has L antennas)
            %        h2h == N x N x M x K matrix
            %                   ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
            %                       h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k|
        H2H{iReal, 1} = h2h;
        R_nm(:,:,iReal) = dA;
    end

    for iReal = 1:noRealizations
        fprintf('iReal:%i/%i    iN:%i/%i',iReal,noRealizations,NoUsers(iN),NoUsers(length(NoUsers)));
        
        t = randi(800, 1);
        f_l = f_user(t:t+users_no-1);
        T_l = C_n./f_l;
        E_l = kappa.*C_n.*(f_l).^2;
        
        eta = beta_t.*D_n./(T_l);
        theta = beta_e.*D_n./(zeta.*E_l);

        h2h_ = H2H{iReal, 1};
        %        h2h_ == N x N x M x K matrix
        %                   ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
        %                       h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k|
        r_nm   = R_nm(:,:,iReal); % N x M matrix
        
        Adet = 0;
        noUsers = size(UEs.active, 1);   % N
        noBSs   = size(BS.positions, 1); % M

        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails2('MEC_NOMA21', noUsers, noSubcs, noBSs, UE_BS, h2h_, ...
            p_min, p_max, P_tol, n0, W, eta, theta, nu, lambda, beta, f0, f_l, Adet);
            % function in ..\
         % MEC_NOMA21
        
        [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, ~, ~, ~, time] = BWOA2(...
              'MEC_NOMA21', doTol, noSearchAgents, noUsers, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, Adet);
            % function in ..\WOA_voronoi
        po_MECNOMA21(iN, iReal) = sum(sum(sum(leader_pos_bwoa)))/users_no;
        su_MECNOMA21(iN, iReal) = leader_score_bwoa;
        time_MECNOMA21(iN, iReal) = time;
        
        % Exhaustive search
          [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, time] = exhaustive2(...
              'MEC_NOMA21', noSearchAgents, noUsers, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
              theta, eta, W, h2h_, n0, p_min, p_max, nu);
            % function in ..\

          po_EX(iN, iReal) = sum(sum(sum(leader_pos_bwoa)))/users_no;
          su_EX(iN, iReal) = leader_score_bwoa;
          time_EX(iN, iReal) = time;
        
    end
end
po_MECNOMA21 = mean(po_MECNOMA21, 2);
su_MECNOMA21 = mean(su_MECNOMA21, 2);
time_MECNOMA21 = mean(time_MECNOMA21, 2);

po_EX = mean(po_EX, 2);
su_EX = mean(su_EX, 2);
time_EX = mean(time_EX, 2);
% NoUsers = 2:2:10;

% save('exhastive.mat', 'NoUsers', 'xt', 'po_MECNOMA21', 'su_MECNOMA21', 'time_MECNOMA21', 'po_EX', 'su_EX', 'time_EX');
% xticks = 2:2:6;
% xt = {'2','3','4', '5', '6'};
% 
% h2 = figure(2)
% hold on
% plot(NoUsers, su_MECNOMA21, 'b-o', 'linewidth', 2.0, 'markers', 13.0);
% plot(NoUsers, su_EX, 'r-s', 'linewidth', 2.0, 'markers', 13.0);
% grid on
% set(gca, 'FontSize', 17.5, 'XLim', [2 10]);
% set(gca, 'xtick', xticks);
% set(gca, 'xticklabel', xt);
% xlabel('Number of Users');
% ylabel('System Utility');
% legend({'MEC NOMA 21', 'EX'}, 'Location', 'Best', 'FontSize', 17.5);
% box on
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% hold
% x = [0.3 0.5], y = [0.6 0.5]
% dim = [0.2 0.2 0.1 0.1]
% annotation('ellipse',dim)
% a = annotation('textarrow',x,y,'String','Optimality gap = 1.22%', 'FontSize', 24);
% savefig(h2, 'ex_su');
% 
% h3 = figure(3)
% ax(1) = axes();
% l1=line('parent',ax(1),'xdata',NoUsers,'ydata',time_MECNOMA21, 'LineWidth', 2.0, 'Marker', 'o', 'MarkerSize', 13.0, 'Color', 'b');
% set(ax(1), 'XLim', [2 7], 'FontSize', 17.5);
% set(ax(1), 'xtick', xticks);
% set(ax(1), 'xticklabel',xt);
% xlabel('parent', ax(1), 'Number of Users');
% ylabel('parent', ax(1), 'Algorithm Runtime (s)');
% grid on
% box on
% 
% ax = ax(1);
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% 
% ax(2) = axes('position', [0.18 0.57 0.25 0.4]);
% l2=line('parent',ax(2),'xdata',NoUsers,'ydata',time_EX,'LineWidth', 2.0, 'Marker', 's', 'MarkerSize', 13.0, 'Color', 'r');
% set(ax(2), 'XLim', [2 7], 'xticklabel', [], 'FontSize', 17.5);
% axis tight
% box on
% grid on
% 
% lngd = legend( [l1;l2] , {'MEC NOMA 21','EX'}, 'FontSize', 20);
% savefig(h3, 'ex_time');

