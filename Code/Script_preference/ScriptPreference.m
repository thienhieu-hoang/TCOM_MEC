%-----------------------------------
% change user's preference in time
%-----------------------------------
tic
clear all
addpath('..\WOA_voronoi\')
addpath('..\Generate\')
addpath('..\')
load('..\Parameters\parameters2.mat')
noRealizations = 30;
users_no = 16; %25
NoUsers = users_no;
Beta_t = 0.1:0.1:0.9; %0.1

noBSs   = 4;
noSubcs = 5;
noAnten = 4;


po_MECNOMA21 = zeros(length(Beta_t), noRealizations);   % Offloading percentage 
su_MECNOMA21 = zeros(length(Beta_t), noRealizations);   % System utility
ti_MECNOMA21 = zeros(length(Beta_t), noRealizations);   % Time
en_MECNOMA21 = zeros(length(Beta_t), noRealizations);   % Energy

xt_beta = cell(1, length(Beta_t));


channelGain = zeros(NoUsers, noBSs, noSubcs, noRealizations);

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

for iK = 1:length(Beta_t)
    beta_t = Beta_t(iK);
    beta_e = 1 - beta_t;
    beta = [beta_t beta_e];
    xt_beta(iK) = {num2str(beta_t)};
    
    for iReal = 1:noRealizations
        fprintf('beta_t = %i/1   iReal = %i/%i\n',beta_t, iReal, noRealizations)
        
        UE_BS = UE_BS_cell{iReal,1};
        h2h_  = H2H{iReal, 1};
        %        h2h_ == N x N x M x K matrix
        %                   ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
        %                       h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k|
        r_nm  = R_nm(:,:,iReal); % N x M matrix
     
        t      = randi(800, 1);
        f_l    = f_user(t:t-1+ users_no);
        T_l    = C_n./f_l;
        E_l    = kappa.*C_n.*(f_l).^2;
        
        eta    = beta_t.*D_n./(T_l);
        theta  = beta_e.*D_n./(zeta.*E_l);
        
        %%%%%%%%%%%%%%%%%%%%%
        %      NOMA MEC 
        %%%%%%%%%%%%%%%%%%%%%
        
        T_r_NOMAMEC21 = zeros(users_no, 1);
        E_r_NOMAMEC21 = zeros(users_no, 1);
        Rij_NOMAMEC21 = zeros(users_no, 1);
        
        Adetermined = 0;
        [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails2('MEC_NOMA21', NoUsers, noSubcs, noBSs, UE_BS, h2h_, ...
             p_min, p_max, P_tol, n0, W, eta, theta, nu, lambda, beta, f0, f_l, Adetermined);
         
        [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, ~, ~, ~, ~] = BWOA2(...
             'MEC_NOMA21', doTol, noSearchAgents, NoUsers, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa,...
             theta, eta, W, h2h_, n0, p_min, p_max, nu, 0);
         

        % offloading vector
        A_MECNOMA21  = sum(sum(leader_pos_bwoa, 2),3); % N x 1 matrix
        off_users_no = sum(A_MECNOMA21);               % number of offloading users 
        
        po_MECNOMA21(iK, iReal) = off_users_no/users_no;
        off_idx = (A_MECNOMA21 == 1);
        loc_idx = (A_MECNOMA21 == 0);
        Z_l_MECNOMA21 = loc_idx.*(beta_t.*T_l + beta_e.*E_l);
        
        tmp1 = sqrt(A_MECNOMA21.*beta_t.*f_l);
        f_i = f0.*(tmp1./sum(tmp1));
        
        for i = 1:users_no
            if A_MECNOMA21(i) == 0
                continue;
            end
            % j: offloading subcs
            temp1 = sum(leader_pos_bwoa, 2);  % N x K matrix
            jidx = find(temp1(i, :) == 1);

            % bs: base station that the UE offloads to
            temp2 = sum(leader_pos_bwoa, 3);  % N x K matrix
            bs    = find(temp2(i, :) == 1);

            % k: offloading users to subcs j that its gain smaller than i's gain
            % h2h_
            xkj = (temp1(:, jidx) > 0) & (diag(h2h_(:, :, bs, jidx)) < h2h_(i, i, bs, jidx));
                    % xkj == N x 1 boolean matrix
            
            % data rate  %% float
            Rij_NOMAMEC21(i) = W*log2(1 + (leader_pos_woa(i)*h2h_(i, i, bs, jidx))/(n0 + ...
                sum(xkj.*leader_pos_woa(:).*diag(h2h_(:, :, bs, jidx)))));
            
            T_r_NOMAMEC21(i) = D_n/Rij_NOMAMEC21(i) + C_n/f_i(i);
            E_r_NOMAMEC21(i) = leader_pos_woa(i)*D_n/(zeta*Rij_NOMAMEC21(i));
        end
        Z_r_NOMAMEC21 = beta_t.*T_r_NOMAMEC21 + beta_e*E_r_NOMAMEC21;
        
        T_i = loc_idx.*T_l + off_idx.*T_r_NOMAMEC21;
        E_i = loc_idx.*E_l + off_idx.*E_r_NOMAMEC21;
        
        Z_i = beta_t.*(T_l - T_i)./T_l + beta_e.*(E_l - E_i)./E_l;
        
        su_MECNOMA21(iK, iReal) = sum(Z_i);
        ti_MECNOMA21(iK, iReal) = sum(T_i);
        en_MECNOMA21(iK, iReal) = sum(E_i);
    end
end

po_MECNOMA21_1 = mean(po_MECNOMA21, 2);
su_MECNOMA21_1 = mean(su_MECNOMA21, 2);
ti_MECNOMA21_1 = mean(ti_MECNOMA21, 2);
en_MECNOMA21_1 = mean(en_MECNOMA21, 2);

xt_beta = cell(1, length(Beta_t));
for i = 1:2:length(Beta_t)
    xt_beta(i) = {num2str(Beta_t(i))};
end

save('preference.mat', 'Beta_t' ,'xt_beta','po_MECNOMA21_1', 'su_MECNOMA21_1', 'ti_MECNOMA21_1', 'en_MECNOMA21_1','users_no');

hold

toc
t = toc