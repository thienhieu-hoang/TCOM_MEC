% ---------------------------------------------------------
%WOA for power transmission P == N x 1 matrix
%----------------------------------------------------------

% Output:
% leaderScore: value of obj function after this code == double
% leaderPos == N x 1 matrix
% convergenceCurve == 1 x maxIter  matrix = value of obj function after each iteration

function [leaderScore, leaderPos, convergenceCurve] = WOA_matrix(noSearchAgents, noUsers, maxIter, lb, ub, fobj, X) 
	% Input:
	% noSearchAgents: number of whales
	% noUsers = N  = number of UEs
	% ub == N x 1 matrix == upper bound transmit power p_i^{max}
	% lb == N x 1 matrix == lower bound transmit power p_i^{min}
	% X  == N x M+1 x K matrix == association matrix
    
%    maxIter     = 300; % 150
	leaderPos   = zeros(noUsers, 1); 
	leaderScore = inf; 

	leader_score_pre = leaderScore; 

	convergenceCurve = zeros(1, maxIter); 
 
	% ======================== Initialization =================

		% If each variable has a different lb and ub
		posi_p = 0.5*ones(noUsers, noSearchAgents).*(ub - lb) + lb; %rand(noUsers, noSearchAgents).*(ub - lb) + lb; 
    
  	% ======================== Loop ===========================
	% Loop counter 
	t = 0;
	todoTol = 0; % =0 to run all iteration
	delta = 1e-7; 
	flag = 0; 
    PS_flag = 0; % population of search agents flags
    PS_max = 40;
    PS_min = 10;
    alpha = 0.25;

	% Main loop
	while t < maxIter && flag < 10

			% Return back the search agents that go beyond the boundaries of the search space
			tmp = posi_p; 
			flag4lb = tmp < lb; 
			flag4ub = tmp > ub; 
			posi_p = tmp.*(~(flag4lb + flag4ub)) + lb.*flag4lb + ub.*flag4ub; 

			% Calculate objective function for each search agent
		for i = 1:noSearchAgents
			fitness = fobj(posi_p(:, i), X); 

			% Update the leader 
			if fitness < leaderScore
				leaderScore = fitness; 
				leaderPos = posi_p(:, i);
			end 
		end
	

		% a decreases linearly from 2 to 0
		a = 2 - t*(2/maxIter); 		
		% a2 linearly decreases from -1 to -2 to calculate t 
		a2 = -1 + t*(-1/maxIter); 

		% Update the position of each search agents
		for i = 1:noSearchAgents
			r1 = rand(); 
			r2 = rand();  
			
			A = 2*a*r1 - a; 
			C = 2*r2;

			% parameters for spiral updating position
			b = 1; 
			l = (a2 - 1)*rand + 1; 

			p = rand(); 

				% follow the shrinking encircling mechanism or prey search
				if p < 0.5
					% search for prey (exploration phase)
					if abs(A) >= 1
						randLeaderIndex = floor(noSearchAgents*rand + 1); 
						X_rand = posi_p(:, randLeaderIndex); 		% -> X_rand == N x 1 matrix
						D_X_rand = abs(C*X_rand - posi_p(:,i)); 	%  N x 1
						posi_p(:, i) = X_rand - A*D_X_rand;         %  N x 1
					elseif abs(A) < 1
						D_Leader = abs(C*leaderPos - posi_p(:, i));   	% D_Leader, leaderPos == N x 1
						posi_p(:,i) = leaderPos - A*D_Leader;    % N x 1
					end
				elseif p >= 0.5
					distance2Leader = abs(leaderPos - posi_p(:,i)); 
					posi_p(:,i) = distance2Leader*exp(b.*l).*cos(l.*2*pi) + leaderPos; 
				end 
		end

		% increase the iteration index by 1  
		t = t + 1; 
		convergenceCurve(1,t) = leaderScore; 

        leader_pos_SA{t} = leaderPos;   % cell of N x 1 
                                        %       == position of best SA in generation t
        if t >2
            if (sum(leader_pos_SA{t} ~= leader_pos_SA{t-1}) >0 ) && ...
                    (sum(leader_pos_SA{t-1} ~= leader_pos_SA{t-2}) >0) && ... % update 2 gen consecutively 
                        noSearchAgents > PS_min
                PS_flag = "decrease";
            end
            if (sum(leader_pos_SA{t} ~= leader_pos_SA{t-1}) == 0 ) && ... % not update 1 gen
                        noSearchAgents == PS_max
                PS_flag = "decrease";
            end
            if (sum(leader_pos_SA{t} ~= leader_pos_SA{t-1}) ==0 ) && ...
                        noSearchAgents < PS_max
                PS_flag = "increase";
            end
        end
                        

 		if todoTol == 1 && abs(leaderScore - leader_score_pre) < delta
			flag = flag + 1; 
			convergenceCurve = convergenceCurve(1, 1:t);
		else 
			flag = 0; 
        end 
%        fprintf('WOA iter:%i, leaderScore:%i, flag:%i\n', t, leaderScore, flag) 
		leader_score_pre = leaderScore;
    end 
       plot(1:size(convergenceCurve,2),convergenceCurve);
       hold on;
        % plot conver curve of WOA in one iter, 
        % uncomment and set breakpoint to se the figure       

