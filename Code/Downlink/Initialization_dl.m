%----------------------------------------------------------------------------------------
%%% Initialze positions of whales (subchannel assignment in the case of A --> posi_a 
%%%									SBSs' transmiting power in the case of P --> posi_p )
%----------------------------------------------------------------------------------------

% Output:
% posi_a = M x K x noSA  matrix	== positions of whales in 2-D searching for the best associations
% posi_p = M x N x noSA  matrix	== positions of whales in 2-D searching for the best transmission powers

function [posi_a, posi_p] = Initialization_dl(functionname, noUsers, noSubcs, noSBS, UE_BS, P_SBS_min, P_SBS_max, noSearchAgents)
	% functionname 	== 'string'
	% noUsers == N   		   == number of UEs
	% noSubcs == K 			   == number of subchannel
	% noBSs   == M    		   == number of BS (M+1 indicates BS_0 == MBS)
	% UE_BS   == N x M matrix  == binary matrix of relation of UEs and BSs
	                             % get from 'Generate_dl\gen_location_dl.m'
    % P_SBS_max = 1 x M matrix == transmit power budget of SBSs
    % P_SBS_min = double or 1 x M
                        
	posi_a = zeros(noSBS, noSubcs, noSearchAgents); % M x K x noSA matrix
 	posi_p = zeros(noUsers, noSBS, noSearchAgents);	% N x M x noSA matrix
    
    p_min = P_SBS_min.* UE_BS;             % N x M == lower bound of power transmission of M SBSs to N UEs
	p_max = P_SBS_max.* UE_BS;    % N x M == upper bound

	switch functionname
		%%
		case 'MEC_NOMA_DL'
			% setup posi_p
            for nSA = 1:noSearchAgents
			    posi_p(:,:,nSA) = UE_BS .* (1/noUsers *rand(size(posi_p(:,:,nSA))).*(p_max-p_min)+p_min);
            end
            SBS_busy = sum(UE_BS,1)>0; % SBSs that contain offloading UEs

			% setup posi_a
            k = 0;
 		    for i = 1:noSearchAgents
 				for m = 1:noSBS
 					if (SBS_busy(1,m))
 				        posi_a(m,k+1,i) = 1; 
                        k = rem(k+1, noSubcs);
 					end
 				end
            end 
end


		