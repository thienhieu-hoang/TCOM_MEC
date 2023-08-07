noRealizations = 200;
D_n_sample = 0.2:0.2:1; % x1e6
users_no = 10;
su_PSO_BWOA = [7.715240581597731,7.08831042874054,6.724690544317954,6.705444729279829,6.575249614956282]';
su_IWOA_BWOA = [7.553256142938676,7.307140242549743,7.269539593386082,6.678401271614942,6.66717301097138]';
su_WOA_BWOA = [7.603018924256071,7.270345270958984,7.105659668437376,6.950561921658371,6.906620963157385]';
su_OFDMA = [4.741200758872632,4.629090716184119,4.512739499402244,4.511125722085831,4.408880572879769]';
su_IOJOA = [4.223735033793815,4.06210481930834,3.905404011705798,3.624179226243571,3.077570733209216]';
su_ARJOA = [4.718534723927199,3.682967846542516,2.79892469414501,0.78522291795327,0.436521205716032]';
su_ALCA = [0,0,0,0,0]';

save('Script_Datasize.mat',"noRealizations", 'D_n_sample', 'users_no', 'su_PSO_BWOA', ...
    "su_ARJOA" ,'su_IOJOA', 'su_ALCA', 'su_WOA_BWOA', 'su_OFDMA', 'su_IWOA_BWOA');
