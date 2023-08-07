%----------------------------------------------------------------------------
% apply BWOA algorithm and condition 1 to solve UL pb .
%----------------------------------------------------------------------------

% leader_score_bwoa == double
% leader_pos_bwoa 	== N x M+1 x K matrix


function [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, conver_curve, conver_curve_woa, no_WOA_run, time] = BWOA1(functionName, doTol, noSearchAgents, noUsers, noSubcs, noBSs, r_m, r_nm, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa, theta, eta, hArray, n0, p_min, p_max, nu, Adetermined)

tic 
	% initialization ======
	[positions posi_p_nouse] = Initialization2(functionName, noUsers, noSubcs, noBSs, r_m, r_nm, p_max, p_min, noSearchAgents);
		% positions == N x M+1 x K x noSA  matrix == position of noSA binary whales
		% posi_p == N x noSA matrix

	leader_pos_bwoa = zeros(noUsers, noBSs, noSubcs); % position of the whale that makes the obj function get the best fitted value
	leader_score_bwoa = -inf; % value of the best fitted obj function
	leader_score_pre = leader_score_bwoa; 
	% leader_score_woa 
	leader_pos_woa = zeros(noUsers, noBSs, noSubcs); 


	% loop counter
	todoTol = 0; % = doTol == 0 to check all 150 iterations
	delta = 1e-6; 
	flag = 0; 
	
	conver_curve = zeros(1, maxIter);
    conver_curve_woa = zeros(1, maxIter); 
	C_Step = zeros(noUsers, noBSs, noSubcs); 
	iter = 0; 
	no_WOA_run = 0; 

	while iter < maxIter
 %       display(['iter = ' num2str(iter)]);  
         whale_no = 0; 
		for nSA = 1:noSearchAgents
%            fprintf('BWOA iter:%i/%i\n', iter, maxIter)
%            fprintf('SearchAgent no %i/%i\n', nSA, noSearchAgents)
            fitbwoa = fobj_bwoa(positions(:, :, :, nSA)); 
 %           display([k fitbwoa score_bwoa leader_score_bwoa]);

 			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% define constraint to cut down on the number of WOA runs
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			

			% condition 1
			if fitbwoa > 1e2 || fitbwoa <= leader_score_bwoa
				continue
			end 



        
            no_WOA_run = no_WOA_run + 1; 

            [WOA_rs, pos_woa, conver_woa] = WOA(noSearchAgents, ...
				noUsers, maxIter, lb_woa, ub_woa, fobj_woa, positions(:,:,:,nSA));

            fitness = fitbwoa - WOA_rs; % double

			% update the leader
			if fitness >= leader_score_bwoa
				leader_score_bwoa = fitness; 
				leader_pos_bwoa   = positions(:, :, :, nSA); 
				leader_pos_woa 	  = pos_woa; 
				leader_score_woa  = WOA_rs; 
                whale_no          = nSA;
            end 
%            fprintf('BWOA fitness Whale no (SA no): %i\n', whale_no)
%            fprintf('BWOA leaderscore: %i\n',leader_score_bwoa)
		end

		a = 2 - iter*(2/maxIter); % a decreases linearly from 2 to 0
		a2 = -1 + iter*(-1/maxIter); % a2 decreases linearly from -1 to -2 
		% update the position of each search agents 

		for nSA = 1:noSearchAgents
			r1 = rand(); 
			r2 = rand(); 
			A = 2*a*r1 - a; 
			C = 2*r2; 
			% parameters for spiral updating position
			b = 1; 
			l = (a2 - 1)*rand + 1; 
			p = rand();
			for k = 1:noSubcs
				for m = 1:noBSs
			  		for n = 1:noUsers				
 						if p < 0.5
							% search for prey (exploration phase)
							if abs(A) >= 1 
								rand_idx = floor(noSearchAgents*rand + 1); 
								X_rand   = positions(:, :, :, rand_idx); % N x M+1 x K matrix
								D_X_rand = abs(C*X_rand(n, m, k) - positions(n, m, k, nSA)); 
								C_Step(n, m, k) = X_rand(n, m, k) - A*D_X_rand;
							elseif abs(A) < 1
								% shrinking encircling mechanism (exploitation phase)
								D_leader = abs(C*leader_pos_bwoa(n, m, k) - positions(n, m, k, nSA)); 
								C_Step(n, m, k) = leader_pos_bwoa(n, m, k) - A* D_leader; 
							end
						elseif p >= 0.5
							distance2leader = abs(leader_pos_bwoa(n, m, k) - positions(n, m, k, nSA));
							C_Step(n, m, k) = distance2leader*exp(b.*l).*cos(l.*2*pi) + leader_pos_bwoa(n, m, k); 
						end 

						sigmoid = 1/(1 + exp(-10*(C_Step(n, m, k)-0.5))); 

						p_rand = rand(); 
						if p_rand < sigmoid
							positions(n, m, k, nSA) = ~positions(n, m, k, nSA); 
						end 		
					end 
				end 
			end
		end

		iter = iter + 1; 
		conver_curve(iter) = leader_score_bwoa; 
		conver_curve_woa(iter) = leader_score_woa;

		fprintf('iter:%i/%i, leader_score_woa:%i, leader_score_bwoa: %i\n', iter, maxIter, leader_score_woa, leader_score_bwoa)

		if todoTol == 1 && iter > 70 && abs(leader_score_bwoa - leader_score_pre) < delta 
                                 %40
            flag = flag + 1; 
		else 
			flag = 0; 
		end 
		leader_score_pre = leader_score_bwoa; 
		if flag == 15
			conver_curve = conver_curve(1, 1:iter); 
			conver_curve_woa = conver_curve_woa(1:iter); 
			break; 
		end 
	end
	toc
	time = toc;
end 