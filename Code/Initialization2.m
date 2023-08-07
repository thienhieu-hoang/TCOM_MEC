%----------------------------------------------------------------------------------------
%%% Initialze positions of whales (subchannel assignment in the case of A --> posi_a 
%%%									UEs' transmiting power in the case of P --> posi_p )
%----------------------------------------------------------------------------------------

% Output:
% posi_a = N x M x K x noSA  matrix	== positions of whales in 3-D searching for the best associations
% posi_p = N x noSA  matrix			== positions of whales in 1-D searching for the best transmission powers

function [posi_a, posi_p] = Initialization2(functionname, noUsers, noSubcs, noBSs, UE_BS, p_max, p_min, noSearchAgents, Adet)
	% functionname 	== 'string'
	% noUsers = N   		   == number of UEs
	% noSubcs = K 			   == number of subchannel
	% noBSs   = M    		   == number of BS (M+1 indicates BS_0 == MBS)
	% UE_BS   = N x M matrix   == binary matrix of relation of UEs and BSs
                                 % run 'Generate\location_voronoi.m' to get
	% Adet	  = N x 1 matrix

	posi_a = zeros(noUsers, noBSs, noSubcs, noSearchAgents); % N x M x K x noSA matrix
 	posi_p = zeros(noUsers, noSearchAgents);				 % N x noSA matrix

	switch functionname
		%%
		case 'MEC_NOMA21'
			% setup posi_p
			posi_p = rand(size(posi_p)).*(p_max-p_min)+p_min;

            
			% setup posi_a
 		    for i = 1:noSearchAgents
 				for m = noBSs:-1:1   % m = M downto 1
 					for n = 1:noUsers
 						if (UE_BS(n,m))
                            % in the case UE can offload to MBS:
%   						    if (m<noBSs) 
%                                 posi_a(n, noBSs, :, i) = zeros(size(posi_a(n, noBSs, :, i))); 
%   							    % guarantee that UEs in the coverages of SBSs don't offload to MBS
%   						   	end   
 						 	
 						 	if rand > 0.5  % probability of offloading is 50%
 								rand_idx = floor(noSubcs*rand + 1); 
 								posi_a(n, m, rand_idx, i) = 1; 
 
 							end
 						end 
 					end
 				end
            end
            
%            posi_a = zeros(noUsers, noBSs, noSubcs, noSearchAgents);
            
            
 		%% All Remote Joint Optimization Algorithm 
 		case 'ARJOA'
			for i = 1:noSearchAgents
 				for m = noBSs:-1:1   % m = M downto 1
 					for n = 1:noUsers
 						if (UE_BS(n,m))
 						 	
 						 	rand_idx = floor(noSubcs*rand + 1); 
 							posi_a(n, m, rand_idx, i) = 1; 
 
 						end 
 					end
 				end
 			end
 	 	%% Independent Offloading Joint Optimization Algorithm
		case 'IOJOA'
			for i = 1:noSearchAgents
 				for m = noBSs:-1:1   % m = M downto 1
 					for n = 1:noUsers
 						if (UE_BS(n,m))
%   						 	if (m<noBSs) posi_a(n, noBSs, :, i) = zeros(size(posi_a(n, noBSs, :, i))); 
%   							    % guarantee that UEs in the coverages of SBSs don't offload to MBS
%   						   	end   
 						 	
 						 	if Adet(n) == 1
 						 		rand_idx = floor(noSubcs*rand + 1); 
 								posi_a(n, m, rand_idx, i) = 1; 
 							end
 						end 
 					end
 				end
 			end
		%% OFDMA
		case 'OFDMA'
			for i = 1:noSearchAgents  
                [xx,yy] = find(UE_BS==1); % find UE-BS associations
                temp    = min(length(xx), noSubcs);
                subc_   = randperm(noSubcs,temp); % random occupied subcarriers 
                                                  % == 1 x temp vector
                u_     = randperm(length(xx),temp); % random associations
                                                    % UE xx(u_(j)) offloads
                                                    % to BS yy(u_(j))
 				for k = 1:length(subc_)         
                    posi_a(xx(u_(k)), yy(u_(k)), subc_(k), i) = 1;
 				end
 				
 			end
	end
end


		