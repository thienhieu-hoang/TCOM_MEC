%------------------------------------------------------------------
% Obtain the objective function and search domain for WOA and BWOA. 
% Version C4 in fraction form
% matrix in the dimension of N x M+1 x K
%------------------------------------------------------------------

% Output:
% lb_woa = N x 1 matrix = lower bound of power transmission of UEs
% ub_woa = N x 1 matrix = upper bound of power transmission of UEs
% fobj_woa = 'string' = @function name of woa
% fobj_bwoa = 'string' = @function name of bwoa

function [lb_woa, ub_woa, fobj_woa, fobj_bwoa] = getFunctionDetails2(F, noUsers, noSubcs, noBSs, UE_BS, h2h, p_min, p_max, P_tol, n0, W, eta, theta, nu, lambda, beta, f0, f_l, Adet)
% F function name: F: depends on each simulation schemes {ARJOA, ODSTCA, IOJOA, OFDMA} 
% network size: noUsers x noBSs x noSubcs == N x M x K
% UE_BS   == N x M matrix   == binary matrix of relation of UEs and BSs
% h2h     == N x N x M x K matrix
%               ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
%                   h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k|
% users' power budget: [p_min, p_max]
% tolerant power for NOMA: P_tol == 1 x M matrix
% thermal noise: n0
% W: bandwidth of subchannel
% f0 == 1 x M matrix:  total computing resource of BSs
% eta, theta, f_l: parameters of local computing (in (26))  == N x 1 matrix
% penalty factor for constraint dealing: <double> nu (here, in this code; but in pdf, they are nu and lambda)
% perference in time and energy: beta = [beta_t beta_e] == N x 2 matrix
% matrix defines offloading decision of IOJOA: Adet 	== N x 1 matrix
    
%     noUsers = size(UEs.active, 1);   % N
%     noBSs   = size(BS.positions, 1); % M

	lb_woa = zeros(noUsers, 1); 
	ub_woa = zeros(noUsers, 1);

	PMin = zeros(noUsers, 1); 
	PMax = zeros(noUsers, 1); 

	PMin(:) = p_min; 
	PMax(:) = p_max; 

	switch F
	case 'MEC_NOMA21'
		fobj_woa = @FWOA; 
		lb_woa(:) = p_min; 
		ub_woa(:) = p_max;
		fobj_bwoa = @FBWOA_ODSTCA;
	case 'ARJOA'
		fobj_woa = @FWOA;
		lb_woa(:) = p_min;
		ub_woa(:) = p_max; 
		fobj_bwoa = @FBWOA_ARJOA;
	case 'IOJOA'
		fobj_woa = @FWOA; 
		lb_woa(:) = p_min;
		ub_woa(:) = p_max; 
		fobj_bwoa = @FBWOA_IOJOA;
	case 'OFDMA'
		fobj_woa = @FWOA; 
		lb_woa(:) = p_min;
		ub_woa(:) = p_max; 
		fobj_bwoa = @FBWOA_OFDMA;
    end
                    

function o = FWOA(P, X) 
		% P == N x 1 matrix _ matrix of transmission power 
		% X == N x M+1 x K binary matrix _ matrix of association
	o = 0;
     
	rs = 0; 
	% C3:
    fc3 = P - PMax; 
    flag_fc3 = fc3 >= 0; %== G(f(x)) in (31)
    pnal_fc3 = sum(nu.*flag_fc3.*(fc3.^2));  %% penalty function for C3 
    rs = rs + pnal_fc3;
    % mu = hArray/(n0*W); %% normalized channel gain

    g4 = 0;
    for k = 1:noSubcs
		Xk = X(:,:,k);		% Xk == N x M matrix == association matrix regarding to subchannel k
%		mu_k = mu(:,:,k);	%    == N x M matrix == normalized channel gain regarding to subchannel k

		%% Find the BSs that have offloading UEs (via subchannel k)
		Xk_sum = sum(Xk,1); 		%% 1 x M matrix
		BS_no  = find(Xk_sum>0);  	%  1 x ?1 matrix  
									%  == indexes of BSs that have offloading UEs (via subchannel k)
	    if (sum(BS_no)==0) 
            continue 			% no UE offload via subchannel k --> go check k+1
        end

            [UE_off,BS_off] = find(Xk>0);		%% all UEs that use the same subchannel k
												 % UE_off = row index of Xk==1 	%% ?2 x 1 matrix
												 % BS_off = col index			%% ?2 x 1 matrix

            I_nmk = zeros(noUsers,1);  	% intercell and intracell interference 
												% to UEs using subchannel k
											 	% N x 1 matrix
		
            gamma_k = zeros(noUsers, noBSs); % N x M matrix == SNR matrix
												        
		%% Find the UEs that offload tasks to BS_off(i) (via subchannel k)
            for m_idx = 1 : length(BS_no) % consider BSs that have offloading UEs, not all BSs
                BS_idx = BS_no(m_idx);
                UE_off_m = find(Xk(:,BS_idx)>0);	%% indexes of the UEs that offloading tasks to BS_off(m)
                									%% == find A_m == n's where n \in A_m
												 	%  via subchannel k
												 	%  ?3 x 1 matrix
			% calculate the interference that UEs offloading to BS_idx suffer
                for i = 1:length(UE_off_m)		   % UE_off_m(i) = n \in A_m
                    flag_less = zeros(noUsers, 1); % N x 1 matrix  		%% flag to check whether UE (n) satisfy SIC condition
                    for n_ = 1:noUsers
                        flag_less(n_, 1) = h2h(n_,n_,BS_idx,k) <= h2h(UE_off_m(i),UE_off_m(i),BS_idx,k); % consider BS number BS_idx
                                            % ||h_{nm}^k||^2 <= ||h_{UE(i)m}^k||^2
                    end
                    flag_less(UE_off_m(i)) = 0;

                    %interfer_k(UE_off_m(i)) = sum(flag_less.*P.*Hk(:,BS_idx)); % double
%                    I_nmk(UE_off_m(i),1) = flag_less'*(P.*(h2h(UE_off_m(i),:,BS_idx,k)'.^2/(h2h(UE_off_m(i),UE_off_m(i),BS_idx,k))));
                            %  P == N x 1 matrix
                            %  h2h(i,:,BS_idx,k) == 1 x N matrix
                            %  I_nmk == N x 1 matrix
                    I_nmk(UE_off_m(i),1) = (flag_less'*(P.*diag(h2h(:,:,BS_idx,k))));
                     
                    g4		= P_tol(BS_idx)*(flag_less'*(P.*diag(h2h(:,:,BS_idx,k)))) - P(UE_off_m(i))*h2h(UE_off_m(i), UE_off_m(i), BS_idx,k); %double
                    flag_g4 = g4 >0;
                    pnal_g4 = sum(10^16*nu*flag_g4*(g4 ^2)); 	% double
                    rs = rs + pnal_g4;
                end
                % I_nmk == N x 1 matrix
            end

            for bs = 1: noBSs
                gamma_k(:, bs) = (P.*diag(h2h(:,:,bs,k)))./(n0+I_nmk); % N x 1
            end
            % gamma_k == N x M matrix

            Rk = W*log2(1+gamma_k);  % N x M matrix

            Wk = sum((eta + theta.*P).*sum(Xk./Rk,2)); % W in (31)
			  %      (Nx1 + Nx1.* Nx1)____NxM ./ NxM

            rs = rs + Wk; 

    end
   
	o = o + rs; 
end
%%
function o = FBWOA_ODSTCA(X)
	% X == N x M x K matrix = asociation matrix
	
	% constraint dealing 
 		% Constraint 2 
 		A_nm = sum(X,3); 	% == a_nm in C6 == N x M matrix
 		A_n  = sum(A_nm,2); % == a_n  in C2 == N x 1 matrix
 		fc2  = A_n -1;		% == constraint 2
        flag_fc2 = fc2>0;
       
 		pnal_fc2 = sum(nu.*flag_fc2.*(fc2.^2));

 		% Constraint 1a
 		pnal_fc1a = 0;
  		for m1 = 1:noBSs 
%  			X_nm1 = r_nm(:,m1)<r_m(m1); % flag: the UEs that in the cover of BS m1
 			X1 	  = X.*UE_BS(:, m1); % we will determine the pnal_fc1a with X1
 							 % all the associations of UEs (offloading to BS m1) with other BSs 
 							 % __ will contribute to the pnal_fc1a
 			X1(:,m1,:) = zeros(size(X1(:,m1,:))); % does not count the UEs associating with BS m1
 			fc1a 	   = sum(sum(sum(X1)));
 			flag_fc1a  = fc1a>0;
 			pnal_fc1a  = pnal_fc1a + nu.*flag_fc1a.*(fc1a.^2);
 		end
 	% objective function
 		% get from GRA (33)
		tmp = sqrt(beta(:, 1).*f_l); % numerator in (33)
									   % == N x 1 matrix
		A  	= sum(X,3).*tmp;		 % N x M matrix
		A1 	= (sum(A,1)).^2;  		 % 1 x M matrix
		V_AF= sum(A1./f0);
        A2  = sum(sum(X,3),2);

	beta_sum = sum(sum(beta, 2).*A2); % (beta_t + beta_e) in (18) (28)

	o = beta_sum - V_AF - pnal_fc1a - pnal_fc2; % obj function has another subtrahend (value of WOA obj function)
												% but we perform it in line 48 code BWOA.m

end 

%% All Remote Joint Optimization Algorithm
function o = FBWOA_ARJOA(X)
		% X == N x M x K matrix

	% dealing constraints: all UEs have to offload
	A_nm = sum(X,3); 	% == a_nm  == N x M matrix
	A_n  = sum(A_nm,2); % == a_n   == N x 1 matrix
	fc2  = A_n -1; 
	flag_fc2 = fc2~=0;
	pnal_fc2 = sum(nu.*flag_fc2.*(fc2.^2)); 

	pnal_fc1a = 0;
 		for m1 = 1:(noBSs) 
 			%X_nm1 = r_nm(:,m1)<r_m(m1); % flag: the UEs that in the cover of BS m1
 			X1 	  = X.*UE_BS(:,m1); % we will determine the pnal_fc1a with X1
 							 % all the associations of UEs (offloading to BS m1) with other BSs... 
 							 % ... will contribute to the pnal_fc1a
 			X1(:,m1,:) = zeros(size(X1(:,m1,:))); % does not count the UEs associating with BS m1
 			fc1a 	   = sum(sum(sum(X1)));
 			flag_fc1a  = fc1a>0;
 			pnal_fc1a  = pnal_fc1a + nu.*flag_fc1a.*(fc1a.^2);
 		end

	% objective function
		% get from GRA (33)
		tmp = sqrt(beta(:, 1).*f_l); % numerator in (33)
									   % == N x 1 matrix
		A  	= sum(X,3).*tmp;		 % N x M matrix
		A1 	= (sum(A,1)).^2;  		 % 1 x M matrix
		V_AF= sum(A1./f0);

	A2  = sum(sum(X,3),2);
	beta_sum = sum(sum(beta, 2).*A2); % (beta_t + beta_e) in (18) (28)

	o = beta_sum - V_AF - pnal_fc1a - pnal_fc2; % obj function has another subtrahend (value of WOA obj function)
												% but we perform it in line 48 code BWOA.m

end 

%% Independent Offloading Joint Optimization Algorithm
function o = FBWOA_IOJOA(X) % A was determined  
	% constraint dealing 
		% Constraint 2 
 		A_nm = sum(X,3); 	% == a_nm in C6 == N x M matrix
 		A_n  = sum(A_nm,2); % == a_n  in C2 == N x 1 matrix
 			% Constraint 3 (Constraint of IOJOA) <no need C2>
		fc3 = A_n - Adet; 
		flag_fc3 = (fc3 ~= 0); 
		pnal_fc3 = sum(nu.*flag_fc3.*(fc3.^2));

 		% Constraint 1a
 		pnal_fc1a = 0;
 		for m1 = 1:(noBSs) 
 			%X_nm1 = r_nm(:,m1)<r_m(m1); % flag: the UEs that in the cover of BS m1
 			X1 	  = X.*UE_BS(:,m1); % we will determine the pnal_fc1a with X1
 							 % all the associations of UEs (offloading to BS m1) with other BSs 
 							 % __ will contribute to the pnal_fc1a
 			X1(:,m1,:) = zeros(size(X1(:,m1,:))); % does not count the UEs associating with BS m1
 			fc1a 	   = sum(sum(sum(X1)));
 			flag_fc1a  = fc1a>0;
 			pnal_fc1a  = pnal_fc1a + nu.*flag_fc1a.*(fc1a.^2);
 		end

  

	% objective function
		% get from GRA (33)
		tmp = sqrt(beta(:, 1).*f_l); % numerator in (33)
									   % == N x 1 matrix
		A  	= sum(X,3).*tmp;		 % N x M+1 matrix
		A1 	= (sum(A,1)).^2;  		 % 1 x M+1 matrix
		V_AF= sum(A1./f0);	
	A2  = sum(sum(X,3),2);
	beta_sum = sum(sum(beta, 2).*A2); % (beta_t + beta_e) in (18) (28)

	o = beta_sum - V_AF - pnal_fc1a - pnal_fc3; % obj function has another subtrahend (value of WOA obj function)
												% but we perform it in line 48 code BWOA.m

end 

%% OFDMA
function o = FBWOA_OFDMA(X)
	
	% constraint dealing 
 		% Constraint 2 
 		A_nm = sum(X,3); 	% == a_nm in C6 == N x M matrix
 		A_n  = sum(A_nm,2); % == a_n  in C2 == N x 1 matrix
 		fc2  = A_n -1;		% == constraint 2
        flag_fc2 = fc2>0;
       
 		pnal_fc2 = sum(nu.*flag_fc2.*(fc2.^2));

 		% Constraint 1a
 		pnal_fc1a = 0;
 		for m1 = 1:(noBSs-1) 
%  			X_nm1 = r_nm(:,m1)<r_m(m1); % flag: the UEs that in the cover of BS m1
 			X1 	  = X.*UE_BS(:,m1); % we will determine the pnal_fc1a with X1
 							 % all the associations of UEs (offloading to BS m1) with other BSs 
 							 % __ will contribute to the pnal_fc1a
 			X1(:,m1,:) = zeros(size(X1(:,m1,:))); % does not count the UEs associating with BS m1
 			fc1a 	   = sum(sum(sum(X1)));
 			flag_fc1a  = fc1a>0;
 			pnal_fc1a  = pnal_fc1a + nu.*flag_fc1a.*(fc1a.^2);
 		end

 		% Constraint of OFDMA
 		A_k 	  = sum(sum(X,2),1);
 		fc_ofdm   = A_k-1;
 		flag_ofdm = fc_ofdm>0;

 		pnal_ofdm = sum(nu.*flag_ofdm.*(fc_ofdm.^2));


	% objective function 
		% get from GRA (33)
		tmp = sqrt(beta(:, 1).*f_l); % numerator in (33)
									   % == N x 1 matrix
		A  	= sum(X,3).*tmp;		 % N x M+1 matrix
		A1 	= (sum(A,1)).^2;  		 % 1 x M+1 matrix
		V_AF= sum(A1./f0);	
	A2  = sum(sum(X,3),2);
	beta_sum = sum(sum(beta, 2).*A2); % (beta_t + beta_e) in (18) (28)

	o = beta_sum - V_AF - pnal_fc1a - pnal_fc2 - pnal_ofdm; % obj function has another subtrahend (value of WOA obj function)
												% but we perform it in line 48 code BWOA.m

end 
end