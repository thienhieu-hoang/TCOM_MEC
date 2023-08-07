% This script sets the values of input parameters of the system
% you can change their values here, then save to parameters2.mat
% ------------------------------------------------------------------

noSearchAgents = 30;
noSubcs = 5;

logNormalMean = 0;
logNormalDeviation = 8.0;

maxIter = 1500;
noRealizations = 200; %300

beta_t = 0.5;
beta_e = 1 - lambda_t;
beta = [beta_t beta_e];


n0 = db2lin(-114 - 30);
W = 1e6;

p_min = 1e-8;
p_max = 0.25;
f0 = 1e9* [8 8 8 8];
D_n = 1*420e3;
C_n = 1000e6;
kappa = 5e-27;
zeta = 1;

nu = 1e14;
P_tol = [1.001, 1.001, 1.001, 1.001];

maxIter = 100; 		%500
NoUsers = [100];	%5 15 25
noSubcs = 5;
noBSs 	= 5;


f_local = 1e9*[0.5 0.8 1];
f_user = zeros(1000, 1);
for i = 1:1000
    f_user(i) = f_local(randi(length(f_local), 1));
end
