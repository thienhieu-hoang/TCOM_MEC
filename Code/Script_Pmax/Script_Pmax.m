% ----------------
% script: change maximum power of users
% ----------------
tic
clear all
addpath('..\WOA_voronoi\')
addpath('..\Generate\')
addpath('..\')
load('..\Parameters\parameters2.mat')

% noSearchAgents = 30; 
maxIter = 1000;
users_no = 10; 
NoUsers  = users_no;
noBSs   = 4;
noSubcs = 5;
noAnten = 4;

noRealizations = 200; %200

PMax = 17:2:25; %dBm

po_ALCA   = zeros(length(PMax), noRealizations);
su_ALCA   = zeros(length(PMax), noRealizations);

po_ARJOA  = zeros(length(PMax), noRealizations);
su_ARJOA  = zeros(length(PMax), noRealizations);

po_MECNOMA21 = zeros(length(PMax), noRealizations);
su_MECNOMA21 = zeros(length(PMax), noRealizations);

po_IOJOA  = zeros(length(PMax), noRealizations);
su_IOJOA  = zeros(length(PMax), noRealizations);

po_OFDMA  = zeros(length(PMax), noRealizations);
su_OFDMA  = zeros(length(PMax), noRealizations);

po_WOA_BWOA = zeros(length(PMax), noRealizations);
su_WOA_BWOA = zeros(length(PMax), noRealizations);

po_PSO_BWOA = zeros(length(PMax), noRealizations);
su_PSO_BWOA = zeros(length(PMax), noRealizations);

po_IWOA_BWOA = zeros(length(PMax), noRealizations);
su_IWOA_BWOA = zeros(length(PMax), noRealizations);

xt = cell(1, length(PMax));

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

        UE_BS = UE_BS_cell{iReal,1};
        h2h_  = H2H{iReal, 1};
        %        h2h_ == N x N x M x K matrix
        %                   ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
        %                       h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k|
        r_nm  = R_nm(:,:,iReal); % N x M matrix
end


for iN = 1:length(PMax)
    p_max  = PMax(iN);    
    xt(iN) = {num2str(p_max)};
    p_max  = db2lin(p_max - 30); 

    for iReal = 1:noRealizations
        fprintf('Pmax = %i/%i   iReal = %i/%i\n',PMax(iN), PMax(length(PMax)), iReal, noRealizations)
        h2h_  = H2H{iReal, 1};
        %        h2h_ == N x N x M x K matrix
        %                   ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
        %                       h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k|
        r_nm  = R_nm(:,:,iReal); % N x M matrix

        t   = randi(800, 1);
        f_l = f_user(t:t+users_no-1);
        T_l = C_n./f_l;
        E_l = kappa.*C_n.*(f_l).^2;
        
        %     %%%%%%%%%%%%%%%%%%%%%
        %     %       ALCA
        %     %%%%%%%%%%%%%%%%%%%%%
        po_ALCA(iN, iReal) = 0;
        su_ALCA(iN, iReal) = 0;
        
        eta   = beta_t.*D_n./(T_l);
        theta = beta_e.*D_n./(zeta.*E_l);

         
        %     %%%%%%%%%%%%%%%%%%%%%
        %     %       ARJOA
        %     %%%%%%%%%%%%%%%%%%%%%
        
        fprintf('ARJOA\n')
        % Adetermined = 0;
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails2('ARJOA', users_no, noSubcs, noBSs, UE_BS, h2h_, ...
            p_min, p_max, P_tol, n0, W, eta, theta, nu, lambda, beta, f0, f_l, 0);
            % function in ..\

         [leader_score_bwoa, leader_pos_bwoa, ~, ~, ~, ~, ~] = BWOA2(...
              'ARJOA', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
        
        offloading vector
        A_ARJOA     = sum(sum(leader_pos_bwoa, 2),3); % N x 1 matrix
        off_NoUsers = sum(A_ARJOA);

        po_ARJOA(iN, iReal) = off_NoUsers/NoUsers;
        su_ARJOA(iN, iReal) = leader_score_bwoa;
        
        %%%%%%%%%%%%%%%%%%%%%
        %      MF SIC
        %%%%%%%%%%%%%%%%%%%%%
        
        fprintf('WOA BWOA\n')
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails2('MEC_NOMA21', users_no, noSubcs, noBSs, UE_BS, h2h_, ...
            p_min, p_max, P_tol, n0, W, eta, theta, nu, lambda, beta, f0, f_l, 0);
%             function in ..\

         [leader_score_bwoa, leader_pos_bwoa, ~, ~, ~, ~, ~] = BWOA2(...
              'MEC_NOMA21', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
%             function in ..\WOA_voronoi

        % offloading vector
        A_MECNOMA21  = sum(sum(leader_pos_bwoa, 2),3); % N x 1 matrix
        off_users_no = sum(A_MECNOMA21);
         
        po_WOA_BWOA(iN, iReal) = off_users_no/NoUsers;
        su_WOA_BWOA(iN, iReal) = leader_score_bwoa;

                 % IWOA BWOA
        fprintf('IWOA BWOA\n')
        [leader_score_iwoa, leader_pos_iwoa, ~, ~, ~, ~, ~] = IWOA_BWOA(...
              'MEC_NOMA21', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
%             function in ..\WOA_voronoi  

        % offloading vector
        A_IWOA_BWOA  = sum(sum(leader_pos_iwoa, 2),3); % N x 1 matrix
        off_users_no = sum(A_IWOA_BWOA);
         
        po_IWOA_BWOA(iN, iReal) = off_users_no/NoUsers;
        su_IWOA_BWOA(iN, iReal) = leader_score_iwoa;

        %          PSO BWOA
         fprintf('PSO BWOA\n')
         [leader_score_PSO_BWOA, leader_pos_PSO_BWOA, leader_pos_PSO, conver_curve_BWOA_pso, conver_curve_pso, ~, time_pso] = PSO_BWOA(...
              'MEC_NOMA21', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
         po_PSO_BWOA(iN, iReal) = sum(sum(sum(leader_pos_PSO_BWOA)))/users_no;
         su_PSO_BWOA(iN, iReal) = leader_score_PSO_BWOA;
%        

        %%%%%%%%%%%%%%%%%%%%%
        %      IOJOA
        %%%%%%%%%%%%%%%%%%%%%
        fprintf('IOJOA\n')
        Adetermined = zeros(users_no, 1);
        % all users independently make offloading decision
        for i = 1:NoUsers
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

        %%% ==================================
        
        % offloading vector
        
        A_IOJOA     = sum(sum(leader_pos_bwoa, 2),3); % N x 1 matrix
        off_NoUsers = sum(A_IOJOA);

        po_IOJOA(iN, iReal) = off_NoUsers/NoUsers;
        su_IOJOA(iN, iReal) = leader_score_bwoa;
        
        
        %%%%%%%%%%%%%%%%%%%%%
        %      OFDMA
        %%%%%%%%%%%%%%%%%%%%%
        
        fprintf('OFDMA\n')
        Adetermined = 0; 

        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails2('OFDMA', users_no, noSubcs, noBSs, UE_BS, h2h_, ...
            p_min, p_max, P_tol, n0, W, eta, theta, nu, lambda, beta, f0, f_l, 0);
%             function in ..\

         [leader_score_bwoa, leader_pos_bwoa, ~, ~, ~, ~, ~] = BWOA2(...
              'OFDMA', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
%             function in ..\WOA_voronoi

        % offloading vector
        A_OFDMA     = sum(sum(leader_pos_bwoa, 2),3); 
        off_NoUsers = sum(A_OFDMA);

        po_OFDMA(iN, iReal) = off_NoUsers/NoUsers;
        su_OFDMA(iN, iReal) = leader_score_bwoa;
                
    end
end
po_MECNOMA21 = mean(po_MECNOMA21, 2);
su_MECNOMA21 = mean(su_MECNOMA21, 2);

po_ARJOA  = mean(po_ARJOA, 2);
su_ARJOA  = mean(su_ARJOA, 2);

po_IOJOA  = mean(po_IOJOA, 2);
su_IOJOA  = mean(su_IOJOA, 2);

po_OFDMA  = mean(po_OFDMA, 2);
su_OFDMA  = mean(su_OFDMA, 2);

po_ALCA   = mean(po_ALCA,2);
su_ALCA   = mean(su_ALCA,2);

po_WOA_BWOA = mean(po_WOA_BWOA, 2);
su_WOA_BWOA = mean(su_WOA_BWOA, 2);

po_PSO_BWOA = mean(po_PSO_BWOA, 2);
su_PSO_BWOA = mean(su_PSO_BWOA, 2);

po_IWOA_BWOA = mean(po_IWOA_BWOA, 2);
su_IWOA_BWOA = mean(su_IWOA_BWOA, 2);

xt = cell(1, length(PMax));
for i = 1:2:length(PMax)
    xt(i) = {num2str(PMax(i))}; 
end

save('Script_Pmax.mat',"noRealizations", 'users_no' ,'PMax', 'su_PSO_BWOA', ...
    "su_ARJOA" ,'su_IOJOA', 'su_ALCA', 'su_WOA_BWOA', 'su_OFDMA', 'su_IWOA_BWOA');

toc
