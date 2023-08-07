%--------------------------------------------------------------------------------
% Algorithm to solve ODSTCA pb by
% using EXHAUSTIVE SEARCH to solve SA pb, instead of using BWOA. (still using WOA to solve TPC problem)
%--------------------------------------------------------------------------------

% Output:
% leader_score_bwoa	== double		== obtained value of maximum U
% leader_pos_bwoa	== N x M x K 
% leader_pos_woa 	== N x 1
% time				== double

function [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, time] = exhaustive2(F, noSearchAgents, noUsers, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb, ub, fobj, theta, eta, W, h2h, n0, p_min, p_max, nu);
	% Input:
	% F 	  == function
	% noSearchAgents: number of whales 
	% noUsers == N
	% noSubcs == K
	% noBSs   == M
	% maxIter == maxIter of WOA
	% lb 	  == N x 1 matrix == lower bound of transmit power
	% ub	  == N x 1 matrix == upper bound of transmit power
	% beta	  == N x 2 matrix == [beta_t beta_e]
	% f_l	  == N x 1 matrix == computing capacity of UEs
	% f0	  == 1 x M        == total computing resources of the servers
	% UE_BS   = N x M matrix   == binary matrix of relation of UEs and BSs
    
%     noUsers = size(UEs.active, 1);
%     noBSs   = size(BS.positions, 1);

	tic 
	sa 	= zeros(noUsers, noBSs, noSubcs); 
	n 	= 1; 
    cnt = 0; 
	  
	leader_pos_bwoa 	= zeros(noUsers, noBSs, noSubcs); 
	leader_score_bwoa 	= -inf; 
    bwoa                = -inf;
	leader_pos_woa 		= zeros(noUsers, 1); 

	phi = @(y,a,x,eta) y*log2(1 + a*x) - (a/log(2))*(eta + y*x)/(1 + a*x);
    fmin = @(y,a,x,eta) (eta+y*x)/(W*log2(1 + a*x)); 
	varepsilon = 1e-5 ;

	TRY(n); 

	function [] = solution()
		flag = 0;
		% COMPUTING RESOURCE ALLOCATION 
		fitbwoa = fobj_bwoa(sa); 

        
		if fitbwoa > 1e2 || fitbwoa <= leader_score_bwoa
				flag = 1;
		end

		% condition 2
		if (flag == 0)
			WOA_tmp = 0; 
    	    p_tmp = zeros(noUsers,1); 
			for n = 1:noUsers
				for m = 1:noBSs
					for k = 1:noSubcs
						if sa(n, m, k) == 0
							continue; 
						end 
						if phi(theta(n), h2h(n, n, m, k)/(n0), p_max, eta(n)) <= 0
							p_tmp(n) = p_max; 
						else
							p_s = p_min; p_t = p_max; 
							while (abs(p_t - p_s) > varepsilon)
								p_l = (p_t + p_s)/2; 
								if phi(theta(n), h2h(n, n, m, k)/(n0), p_l, eta(n)) <= 0
									p_s = p_l;
								else 
									p_t = p_l; 
								end 
							end 
							p_tmp(n) = (p_s + p_t)/2; 
						end
						WOA_tmp = WOA_tmp + fmin(theta(n), h2h(n, n, m, k)/(n0), p_tmp(n), eta(n));  
					end
				end 
			end
	        if (fitbwoa - WOA_tmp) <= leader_score_bwoa
    	        flag = 1;
    	    end 

    	    if flag == 0
    	    	% TRANSMIT POWER ALLOCATION 
       			[woa_rs, woa_pos, ~] = WOA(noSearchAgents, noUsers, maxIter, lb, ub, fobj, sa); 
        		bwoa = fitbwoa - woa_rs;
            end
            
        end
												
	
		cnt = cnt + 1;
		if bwoa > leader_score_bwoa
			leader_score_bwoa = bwoa; 
			leader_pos_bwoa   = sa; 
			leader_pos_woa 	  = woa_pos; 
            fprintf('exhaustive leader: %i\n', bwoa);
        end 

	end 

    function [] = TRY(n)
    	
    	for k = 1:noSubcs
    		for m = 0:noBSs
    			if ((m==0) & (k==1))
    				if n == noUsers
    					solution();
    				  else
    					TRY(n+1);
    				end
                else 
                    if m==0
                        continue;
                    end
    				sa(n,m,k) = 1;
    				if n == noUsers
    					solution();
    				  else
    					TRY(n+1);
    				end
    				sa(n,m,k) = 0;
    			end
    		end
		end
	end
    toc;  
    time = toc; 
end 
