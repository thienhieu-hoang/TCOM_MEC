%-------------------
% script: change the number of UL SBSs
%-------------------
tic
clear all
addpath('..\WOA_voronoi\')
addpath('..\Generate\')
addpath('..\')
load('..\Parameters\parameters2.mat')
% noSearchAgents = 30;
users_no = 10;
NoBSs   = 3:7;
noSubcs = 5;
noAnten = 4;

noRealizations = 200;
maxIter = 3000; % higher than usual, more SBSs - higher iterations

doTol = 1;

% po: offloading percentage
% su: system utiliy

po_ALCA = zeros(length(NoBSs), noRealizations);
su_ALCA = zeros(length(NoBSs), noRealizations);

po_ARJOA = zeros(length(NoBSs), noRealizations);
su_ARJOA = zeros(length(NoBSs), noRealizations);

po_IWOA_BWOA = zeros(length(NoBSs), noRealizations);
su_IWOA_BWOA = zeros(length(NoBSs), noRealizations);
ti_IWOA_BWOA = zeros(length(NoBSs), noRealizations);

po_PSO_BWOA = zeros(length(NoBSs), noRealizations);
su_PSO_BWOA = zeros(length(NoBSs), noRealizations);
ti_PSO_BWOA = zeros(length(NoBSs), noRealizations);

po_WOA_BWOA = zeros(length(NoBSs), noRealizations);
su_WOA_BWOA = zeros(length(NoBSs), noRealizations);
ti_WOA_BWOA = zeros(length(NoBSs), noRealizations);

po_IOJOA = zeros(length(NoBSs), noRealizations);
su_IOJOA = zeros(length(NoBSs), noRealizations);

po_OFDMA = zeros(length(NoBSs), noRealizations);
su_OFDMA = zeros(length(NoBSs), noRealizations); 
xt = cell(1, length(NoBSs));
for iM = 1:length(NoBSs)
    noBSs = NoBSs(iM);
    xt(iM) = {num2str(noBSs)};

    f0    = 8*10^9 *ones(1,noBSs);
    P_tol = 1.001 *ones(1,noBSs);
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
        fprintf('iM:%i/%i   iReal:%i/%i\n', NoBSs(iM), NoBSs(length(NoBSs)), iReal, noRealizations);
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
         po_ALCA(iM, iReal) = 0;
         su_ALCA(iM, iReal) = 0;
         
         
         
         %     %%%%%%%%%%%%%%%%%%%%%
         %     %       ARJOA
         %     %%%%%%%%%%%%%%%%%%%%%
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
           
           po_ARJOA(iM, iReal) = off_users_no/users_no;
           su_ARJOA(iM, iReal) = leader_score_bwoa;
           
         %%%%%%%%%%%%%%%%%%%%%
         %      MF SIC
         %%%%%%%%%%%%%%%%%%%%%
         fprintf('WOA BWOA\n')
         [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails2('MEC_NOMA21', users_no, noSubcs, noBSs, UE_BS, h2h_, ...
            p_min, p_max, P_tol, n0, W, eta, theta, nu, lambda, beta, f0, f_l, 0);
%             function in ..\

         [leader_score_bwoa, leader_pos_bwoa, leader_pos_WOA, conver_curve_BWOA_woa, conver_curve_woa, ~, time_woa] = BWOA2(...
              'MEC_NOMA21', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
%             function in ..\WOA_voronoi

%          offloading vector
         A_MECNOMA21 = sum(sum(leader_pos_bwoa, 2),3); % N x 1 matrix
         
         po_WOA_BWOA(iM, iReal) = sum(sum(sum(leader_pos_bwoa)))/users_no;
         su_WOA_BWOA(iM, iReal) = leader_score_bwoa;
         ti_WOA_BWOA(iM, iReal) = time_woa;
% 
        %          IWOA - BWOA
         fprintf('IWOA BWOA\n')
         [leader_score_IWOA_BWOA, leader_pos_IWOA_BWOA, leader_pos_IOWA, conver_curve_BWOA_iwoa, conver_curve_iwoa, ~, time_iwoa] = IWOA_BWOA(...
              'MEC_NOMA21', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
%             function in ..\WOA_voronoi         

         po_IWOA_BWOA(iM, iReal) = sum(sum(sum(leader_pos_IWOA_BWOA)))/users_no;
         su_IWOA_BWOA(iM, iReal) = leader_score_IWOA_BWOA;
         ti_IWOA_BWOA(iM, iReal) = time_iwoa;         
%        
        %          PSO BWOA
         fprintf('PSO BWOA\n')
         [leader_score_PSO_BWOA, leader_pos_PSO_BWOA, leader_pos_PSO, conver_curve_BWOA_pso, conver_curve_pso, ~, time_pso] = PSO_BWOA(...
              'MEC_NOMA21', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
%             function in ..\WOA_voronoi         

         po_PSO_BWOA(iM, iReal) = sum(sum(sum(leader_pos_PSO_BWOA)))/users_no;
         su_PSO_BWOA(iM, iReal) = leader_score_PSO_BWOA;
         ti_PSO_BWOA(iM, iReal) = time_pso;

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

        A_IOJOA      = sum(sum(leader_pos_bwoa, 2),3); % N x 1 matrix
        off_users_no = sum(A_IOJOA);
        
        po_IOJOA(iM, iReal) = off_users_no/users_no;
        su_IOJOA(iM, iReal) = leader_score_bwoa;
        
        
        %%%%%%%%%%%%%%%%%%%%%
        %      OFDMA
        %%%%%%%%%%%%%%%%%%%%%
%         fprintf('OFDMA\n')
% 
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails2('OFDMA', users_no, noSubcs, noBSs, UE_BS, h2h_, ...
            p_min, p_max, P_tol, n0, W, eta, theta, nu, lambda, beta, f0, f_l, 0);
%             function in ..\

         [leader_score_bwoa, leader_pos_bwoa, ~, ~, ~, ~, ~] = BWOA2(...
              'OFDMA', doTol, noSearchAgents, users_no, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
                 theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
%             function in ..\WOA_voronoi
        
        % offloading vector
        A_OFDMA      = sum(sum(leader_pos_bwoa, 2),3); 
        off_users_no = sum(A_OFDMA);
        
        po_OFDMA(iM, iReal) = off_users_no/users_no;
        su_OFDMA(iM, iReal) = leader_score_bwoa;
               
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
po_ARJOA = mean(po_ARJOA, 2);
su_ARJOA = mean(su_ARJOA, 2);

po_IOJOA = mean(po_IOJOA, 2);
su_IOJOA = mean(su_IOJOA, 2);

po_OFDMA = mean(po_OFDMA, 2);
su_OFDMA = mean(su_OFDMA, 2);

po_ALCA = mean(po_ALCA,2);
su_ALCA = mean(su_ALCA,2);

save('Script_M.mat',"noRealizations", 'NoBSs', 'users_no', 'po_ALCA', 'po_WOA_BWOA', 'po_ARJOA', 'po_IOJOA', 'po_OFDMA', ...
    'su_ALCA', 'su_WOA_BWOA', 'su_ARJOA', 'su_OFDMA', 'su_IOJOA', 'po_PSO_BWOA', 'su_PSO_BWOA', 'po_IWOA_BWOA', 'su_IWOA_BWOA');
%

toc

