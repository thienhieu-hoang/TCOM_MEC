%-------------------
% script: change the number of users
%-------------------
tic
clear all
addpath('..\WOA_voronoi\')
addpath('..\Generate\')
addpath('..\')
load('..\Parameters\parameters2.mat')
% noSearchAgents = 30;
NoUsers = 4:6:28; %4:6:28 ;
noBSs   = 4;
noSubcs = 5;
noAnten = 4;

noRealizations = 200; %200

doTol = 1; % result tolerant: 1 for early break / 0 to run all iterations

% po: offloading percentage
% su: system utiliy

po_ALCA = zeros(length(NoUsers), noRealizations);
su_ALCA = zeros(length(NoUsers), noRealizations);

po_ARJOA = zeros(length(NoUsers), noRealizations);
su_ARJOA = zeros(length(NoUsers), noRealizations);

po_IWOA_BWOA = zeros(length(NoUsers), noRealizations);
su_IWOA_BWOA = zeros(length(NoUsers), noRealizations);
ti_IWOA_BWOA = zeros(length(NoUsers), noRealizations);

po_PSO_BWOA = zeros(length(NoUsers), noRealizations);
su_PSO_BWOA = zeros(length(NoUsers), noRealizations);
ti_PSO_BWOA = zeros(length(NoUsers), noRealizations);

po_WOA_BWOA = zeros(length(NoUsers), noRealizations);
su_WOA_BWOA = zeros(length(NoUsers), noRealizations);
ti_WOA_BWOA = zeros(length(NoUsers), noRealizations);

po_IOJOA = zeros(length(NoUsers), noRealizations);
su_IOJOA = zeros(length(NoUsers), noRealizations);

po_OFDMA = zeros(length(NoUsers), noRealizations);
su_OFDMA = zeros(length(NoUsers), noRealizations); 
xt = cell(1, length(NoUsers));
for iN = 1:length(NoUsers)
    users_no = NoUsers(iN);
    xt(iN) = {num2str(users_no)};
    channelGain = zeros(users_no, noBSs, noSubcs, noRealizations);
    
    H2H = cell(noRealizations,1); 
    % channelGain == noRealizations x 1 cell
    %                  each cell is a N x N x M x K matrix (h2h), 
    %                   ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
    %                       h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k| 
    R_nm   = zeros(users_no, noBSs, noRealizations);
    UE_BS_cell = cell(noRealizations, 1); % == noRealizations x 1 cell
                                          %    to save UE_BS of each realization

    for iReal = 1:noRealizations
        [UE_BS_, UEs, BS] = location_voronoi(users_no, noBSs, 0);
        [h2h, hA, dA] = channelMod(UEs, BS, noAnten, noSubcs, logNormalMean, logNormalDeviation);
            % function in ..\
            % return hA  == N x M x K cell, 
            %               each cell is a L x 1 vector  == vector of channel gain
            %                (each SBS has L antennas)
            %        h2h == N x N x M x K matrix
            %                   ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
            %                       h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k|
        H2H{iReal, 1} = h2h;
        UE_BS_cell{iReal, 1} = UE_BS_;
        R_nm(:,:,iReal) = dA;
    end
    
    for iReal = 1:noRealizations
        fprintf('iN:%i/%i   iReal:%i/%i\n', NoUsers(iN), NoUsers(length(NoUsers)), iReal, noRealizations);
        t = randi(800, 1);
        f_l = f_user(t:t+users_no-1);
        T_l = C_n./f_l;
        E_l = kappa.*C_n.*(f_l).^2;

        UE_BS = UE_BS_cell{iReal,1};
        h2h_  = H2H{iReal, 1};
        %        h2h_ == N x N x M x K matrix
        %                   ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
        %                       h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k|
        r_nm  = R_nm(:,:,iReal); % N x M matrix
        

         eta     = lambda_t.*D_n./(T_l);
         theta   = lambda_e.*D_n./(zeta.*E_l);

        %     %%%%%%%%%%%%%%%%%%%%%
        %     %       ALCA
        %     %%%%%%%%%%%%%%%%%%%%%
         fprintf('ALCA\n');
         po_ALCA(iN, iReal) = 0;
         su_ALCA(iN, iReal) = 0;
         
         
         
             %%%%%%%%%%%%%%%%%%%%%
             %       ARJOA
             %%%%%%%%%%%%%%%%%%%%%
          fprintf('ARJOA\n')
           [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails2('ARJOA', users_no, noSubcs, noBSs, UE_BS, h2h_, ...
            p_min, p_max, P_tol, n0, W, eta, theta, nu, lambda, beta, f0, f_l, 0);
            % function in ..\
           
           [leader_score_bwoa, leader_pos_bwoa, ~, ~, ~, ~, ~] = BWOA2(...
              'ARJOA', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
           
           % offloading vector
           A_ARJOA      = sum(sum(leader_pos_bwoa, 2),3); % N x 1 matrix
           off_users_no = sum(A_ARJOA);
           
           po_ARJOA(iN, iReal) = off_users_no/users_no;
           su_ARJOA(iN, iReal) = leader_score_bwoa;
           
         %%%%%%%%%%%%%%%%%%%%%
         %      MF SIC
         %%%%%%%%%%%%%%%%%%%%%
         fprintf('MEC NOMA21\n')
         [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails2('MEC_NOMA21', users_no, noSubcs, noBSs, UE_BS, h2h_, ...
            p_min, p_max, P_tol, n0, W, eta, theta, nu, lambda, beta, f0, f_l, 0);
%             function in ..\

         [leader_score_bwoa, leader_pos_bwoa, leader_pos_WOA, conver_curve_BWOA_woa, conver_curve_woa, ~, time_woa] = BWOA2(...
              'MEC_NOMA21', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
%             function in ..\WOA_voronoi

%          offloading vector
         A_MECNOMA21 = sum(sum(leader_pos_bwoa, 2),3); % N x 1 matrix
%          off_users_no = sum(A_MECNOMA21);
         
         po_WOA_BWOA(iN, iReal) = sum(sum(sum(leader_pos_bwoa)))/users_no;
         su_WOA_BWOA(iN, iReal) = leader_score_bwoa;
% 
        %          IWOA - BWOA
         fprintf('IWOA BWOA\n')
         [leader_score_IWOA_BWOA, leader_pos_IWOA_BWOA, leader_pos_IOWA, conver_curve_BWOA_iwoa, conver_curve_iwoa, ~, time_iwoa] = IWOA_BWOA(...
              'MEC_NOMA21', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);

%             function in ..\WOA_voronoi         
         po_IWOA_BWOA(iN, iReal) = sum(sum(sum(leader_pos_IWOA_BWOA)))/users_no;
         su_IWOA_BWOA(iN, iReal) = leader_score_IWOA_BWOA;
%        
        %          PSO BWOA
         fprintf('PSO BWOA\n')
         [leader_score_PSO_BWOA, leader_pos_PSO_BWOA, leader_pos_PSO, conver_curve_BWOA_pso, conver_curve_pso, ~, time_pso] = PSO_BWOA(...
              'MEC_NOMA21', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
         po_PSO_BWOA(iN, iReal) = sum(sum(sum(leader_pos_PSO_BWOA)))/users_no;
         su_PSO_BWOA(iN, iReal) = leader_score_PSO_BWOA;


        %%%%%%%%%%%%%%%%%%%%%
        %      IOJOA
        %%%%%%%%%%%%%%%%%%%%%
        fprintf('IOJOA\n')
        Adetermined = zeros(users_no, 1);
        % all users independently make offloading decision
        for i = 1:users_no
            if rand() > 0.5
                continue
            end
            Adetermined(i) = 1;
        end
        
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails2('IOJOA', users_no, noSubcs, noBSs, UE_BS, h2h_, ...
            p_min, p_max, P_tol, n0, W, eta, theta, nu, lambda, beta, f0, f_l, Adetermined);
            % function in ..\

        [leader_score_bwoa, leader_pos_bwoa, ~, ~, ~, ~, ~] = BWOA2(...
              'IOJOA', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, Adetermined);
            % function in ..\WOA_voronoi
        leader_score_bwoa

        A_IOJOA      = sum(sum(leader_pos_bwoa, 2),3); % N x 1 matrix
        off_users_no = sum(A_IOJOA);
        
        po_IOJOA(iN, iReal) = off_users_no/users_no;
        su_IOJOA(iN, iReal) = leader_score_bwoa;
        
        
        %%%%%%%%%%%%%%%%%%%%%
        %      OFDMA
        %%%%%%%%%%%%%%%%%%%%%
%         fprintf('OFDMA\n')
% 
%         [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails2('OFDMA', users_no, noSubcs, noBSs, UE_BS, h2h_, ...
%             p_min, p_max, P_tol, n0, W, eta, theta, nu, lambda, beta, f0, f_l, 0);
% %             function in ..\
% 
%          [leader_score_bwoa, leader_pos_bwoa, ~, ~, ~, ~, ~] = BWOA2(...
%               'OFDMA', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
%                  theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
% %             function in ..\WOA_voronoi
%          leader_score_bwoa
%         
%         % offloading vector
%         A_OFDMA      = sum(sum(leader_pos_bwoa, 2),3); 
%         off_users_no = sum(A_OFDMA);
%         
%         po_OFDMA(iN, iReal) = off_users_no/users_no;
%         su_OFDMA(iN, iReal) = leader_score_bwoa;
               
    end
end
po_WOA_BWOA = mean(po_WOA_BWOA, 2);
su_WOA_BWOA = mean(su_WOA_BWOA, 2);
ti_WOA_BWOA = mean(ti_WOA_BWOA, 2); 

po_IWOA_BWOA = mean(po_IWOA_BWOA, 2);
su_IWOA_BWOA = mean(su_IWOA_BWOA, 2);
ti_IWOA_BWOA = mean(ti_IWOA_BWOA, 2);

po_PSO_BWOA = mean(po_PSO_BWOA, 2);
su_PSO_BWOA = mean(su_PSO_BWOA, 2);
ti_PSO_BWOA = mean(ti_PSO_BWOA, 2);

% % 
% po_ARJOA = mean(po_ARJOA, 2);
% su_ARJOA = mean(su_ARJOA, 2);
% 
% po_IOJOA = mean(po_IOJOA, 2);
% su_IOJOA = mean(su_IOJOA, 2);
% 
% po_OFDMA = mean(po_OFDMA, 2);
% su_OFDMA = mean(su_OFDMA, 2);
% 
% po_ALCA = mean(po_ALCA,2);
% su_ALCA = mean(su_ALCA,2);

% save('Script_N_4_6_28_1Realv2_compareAll.mat',"noRealizations" ,'NoUsers', 'po_ALCA', 'po_WOA_BWOA', 'po_ARJOA', 'po_IOJOA', 'po_OFDMA', ...
%     'su_ALCA', 'su_WOA_BWOA', 'su_ARJOA', 'su_OFDMA', 'su_IOJOA', 'po_PSO_BWOA', 'su_PSO_BWOA', 'po_IWOA_BWOA', 'su_IWOA_BWOA');
% save('Script_N_compare_WOA_IWOA_10Real.mat',"noRealizations" ,'NoUsers', 'po_WOA_BWOA', 'su_WOA_BWOA', 'ti_WOA_BWOA', 'po_IWOA_BWOA', 'su_IWOA_BWOA', 'ti_IWOA_BWOA');
save('Script_N_compare_oldChan_WOA_IWOA_PSO_10_28_5Real.mat',"noRealizations" ,'NoUsers', 'po_WOA_BWOA', 'su_WOA_BWOA', 'ti_WOA_BWOA', 'po_IWOA_BWOA', 'su_IWOA_BWOA', 'ti_IWOA_BWOA',"ti_PSO_BWOA",...,
        "ti_PSO_BWOA", 'su_PSO_BWOA', "po_PSO_BWOA", 'conver_curve_BWOA_woa', 'conver_curve_woa', 'time_woa',...
        'leader_pos_IOWA', 'conver_curve_BWOA_iwoa', 'conver_curve_iwoa', 'time_iwoa',...
        'leader_pos_PSO', 'conver_curve_BWOA_pso', 'conver_curve_pso', 'time_pso');


%NoUsers = 4:2:28; 
%xt = {'4', '8', '12', '16', '20', '24', '28'}; 

% h2 = figure(2)
% hold on
% plot(1:length(po_MECNOMA21), po_MECNOMA21(1:length(po_MECNOMA21)), 'b-o', 'linewidth', 2.0, 'markers', 13.0);
% plot(1:length(po_ARJOA), po_ARJOA(1:length(po_ARJOA)), 'g-v', 'linewidth', 2.0 , 'markers', 13.0);
% plot(1:length(po_IOJOA), po_IOJOA(1:length(po_IOJOA)), 'k-x', 'linewidth', 2.0, 'markers', 13.0);
% plot(1:length(po_OFDMA), po_OFDMA(1:length(po_OFDMA)), 'r-s', 'linewidth', 2.0, 'markers', 13.0);
% plot(1:length(po_ALCA), po_ALCA(1:length(po_ALCA)), 'k--d', 'linewidth', 2.0, 'markers', 13.0);
% grid on
% set(gca, 'FontSize', 17.5, 'XLim', [1 length(NoUsers)]);
% xticks = 1:length(NoUsers);
% set(gca, 'xtick', xticks);
% set(gca, 'xticklabel',xt)
% xlabel('Number of Users');
% ylabel('Offloading Percentage');
% lgnd = legend({'MEC NOMA21', 'ARJOA', 'IOJOA', 'OFDMA', 'ALCA'}, 'Location', 'Best')
% temp = [lgnd; lgnd.ItemText];
% set(temp, 'FontSize', 17.5); 
% box on
% % trick: save plot with minimal white space in matlab 
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% hold
% 
% savefig(h2, 'b1000_po'); 
% 
% 
% h3 = figure(3)
% hold on
% plot(1:length(su_MECNOMA21), su_ODSTCA(1:length(su_MECNOMA21)), 'b-o', 'linewidth', 2.0 , 'markers', 13.0);
% plot(1:length(su_ARJOA), su_ARJOA(1:length(su_ARJOA)), 'g-v', 'linewidth', 2.0, 'markers', 13.0);
% plot(1:length(su_IOJOA), su_IOJOA(1:length(su_IOJOA)), 'k-x', 'linewidth', 2.0, 'markers', 13.0);
% plot(1:length(su_OFDMA), su_OFDMA(1:length(su_OFDMA)), 'r-s', 'linewidth', 2.0, 'markers', 13.0);
% plot(1:length(su_ALCA), su_ALCA(1:length(su_ALCA)), 'k--d', 'linewidth', 2.0, 'markers', 13.0);
% grid on
% set(gca, 'FontSize',17.5,  'XLim', [1 length(NoUsers)], 'YLim', [-Inf Inf]);
% xticks = 1:length(NoUsers);
% set(gca, 'xtick', xticks);
% set(gca, 'xticklabel',xt);
% xlabel('Number of Users');
% ylabel('System Utility');
% lgnd = legend({'MEC NOMA21', 'ARJOA', 'IOJOA', 'OFDMA', 'ALCA'}, 'Location', 'Best')
% temp = [lgnd; lgnd.ItemText];
% set(temp, 'FontSize', 17.5); 
% hold
% box on 
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% savefig(h3, 'b1000_su')

toc

