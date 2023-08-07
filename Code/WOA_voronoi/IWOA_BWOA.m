%----------------------------------------------------------------------------
% apply BWOA algorithm and condition 1 and 2 to solve UL pb.
%----------------------------------------------------------------------------

% leader_score_bwoa == double
% leader_pos_bwoa 	== N x M x K matrix
% leader_pos_woa    == N x 1 matrix

function [leader_score_bwoa, leader_pos_bwoa, leader_pos_woa, conver_curve, conver_curve_woa, no_WOA_run, time] = IWOA_BWOA(functionName, doTol, noSearchAgents, noUsers, noSubcs, noBSs, UE_BS, maxIter, fobj_bwoa, lb_woa, ub_woa, fobj_woa, theta, eta, W, h2h, n0, p_min, p_max, nu, Adetermined)

tic 

	% initialization ======
	[positions, ~] = Initialization2(functionName, noUsers, noSubcs, noBSs, UE_BS, p_max, p_min, noSearchAgents, Adetermined);
		% positions == N x M x K x noSA  matrix == position of noSA binary whales
		% posi_p == N x noSA matrix

	leader_pos_bwoa = zeros(noUsers, noBSs, noSubcs); % position of the whale that makes the obj function get the best fitted value
	leader_score_bwoa = -inf; % value of the best fitted obj function
	leader_score_pre = leader_score_bwoa; 
	% leader_score_woa 
	leader_pos_woa = zeros(noUsers, 1); 


	% loop counter
	doTol = 0; % = doTol == 0 to check all 150 iterations
	delta = 1e-5; 
	flag = 0; 
	
	conver_curve = zeros(1, maxIter);
    conver_curve_woa = zeros(1, maxIter); 
	C_Step = zeros(noUsers, noBSs, noSubcs); 
	iter = 0; 
	phi = @(y,a,x,eta) y*log2(1 + a*x) - (a/log(2))*(eta + y*x)/(1 + a*x);
    fmin = @(y,a,x,eta) (eta+y*x)/(W*log2(1 + a*x)); 
	

	varepsilon = 1e-5 ;
 
	no_WOA_run = 0; 

	while iter < maxIter
 %       display(['iter = ' num2str(iter)]);  
 %        whale_no = 0; 
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

			% condition 2
            %tic
			WOA_tmp = 0; 
            p_tmp = zeros(noUsers,1); 
			for n = 1:noUsers
				for m = 1:noBSs
					for k = 1:noSubcs
						if positions(n, m, k, nSA) == 0
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
                continue
            end 
            %toc
            %t1 = toc;
 %           fprintf('old_time for condition 2: %i', t1);


        
            no_WOA_run = no_WOA_run + 1; 

            [WOA_rs, pos_woa, ~] = IWOA(noSearchAgents, ...
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
        
        %tic
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
        %toc
 %       t2 = toc;
 %       fprintf('old_time for searching is: %i', t2);

        
		iter = iter + 1;
%         fprintf('iter %i\n', iter);
		conver_curve(iter) = leader_score_bwoa; 
		conver_curve_woa(iter) = leader_score_woa;

%		fprintf('iter:%i/%i, leader_score_woa:%i, leader_score_bwoa: %i\n', iter, maxIter, leader_score_woa, leader_score_bwoa)

		if doTol == 1 && iter > 600 && abs(leader_score_bwoa - leader_score_pre) < delta 
                                 %600 for general, 100 for comparing
                                 %old_new
            flag = flag + 1; 
		else 
			flag = 0; 
		end 
		leader_score_pre = leader_score_bwoa; 
		if flag == 200
                %200 for general, 20 for comparing old_new
			conver_curve = conver_curve(1, 1:iter); 
			conver_curve_woa = conver_curve_woa(1:iter); 
			break; 
		end 
	end
	toc
	time = toc;
end 