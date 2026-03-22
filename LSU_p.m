function [kf] = LSU_p(Y)
verbose_hysime='on';
noise_type_hysime = 'additive';        
[w_hysime Rn_hysime] = estNoise(Y,noise_type_hysime,verbose_hysime);
[kf,Ek,E_hysime]=hysime(Y,w_hysime,Rn_hysime,1e-5,verbose_hysime); 

end