
%%
%---------------------------------------------
%  If this code is helpful for you, please cite our paper as follows:
%    @article{wu2026endmember,
%    title={Endmember Selection with Adaptive Double Prior Model},
%    author={Wu, Rui and Luo, Wenfei and Gao, Lianru and Ren, Longfei and Gao, Hongmin},
%    journal={IEEE Transactions on Geoscience and Remote Sensing},
%    volume={64},
%    pages={1-22}
%    year={2026},
%    publisher={IEEE}
%    }
%%
%---------------------------------------------
%  Cuprite
%---------------------------------------------
clc
clear 
mypath = '.\\Cuprite\\';
lib_name = strcat(mypath,'Cuprite_lib.mat');
load(lib_name);
big_ed = Cuprite_lib;
big_ed([1:3,107:113,150:170,223:224],:) = [];


image_name = strcat(mypath,'\Cuprite_image.mat');
load(image_name);
image = Cuprite_image;
image = image / 10000;
image =image';
image([1:3,107:113,150:170,223:224],:) = [];

AL_iters = 1000;
tol = 0;



% ADSU
lambda = 0.01;
w = ones(size(fd));
pxl_site = [449,450,451,452,968,969,970,442,443,444,447,448,1004,35:57,895:896,898:899,965,161:162,137:138,912,528:542];
w(pxl_site,:) = 0;
w = normWeight(w);
tic
[x_ad] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);                                 
x_ad_name = strcat(mypath,'\adsu.dat');                
save(x_ad_name,'x_ad');
clear x_ad
toc


% SU
lambda = 0.001;
w = ones(size(fd));
tic
[x_su] =  sunsal_w(big_ed,image,w,'lambda',lambda,'ADDONE','no','POSITIVITY','yes', ...
                        'TOL',tol, 'AL_iters',AL_iters,'verbose','yes'); 
x_name = strcat(mypath,'\su.dat');                
save(x_name,'x_su');
toc
clear x_su

% WSU 
lambda = 0.01;
w = ones(size(fd));
pxl_site = [449,450,451,452,968,969,970,442,443,444,447,448,1004,35:57,895:896,898:899,965,161:162,137:138,912,528:542];
w(pxl_site,:) = 0;
w = normWeight(w);
tic
[x_w] = sunsal_w(big_ed,image,w,'lambda',lambda,'ADDONE','no','POSITIVITY','yes', ...
                        'TOL',tol, 'AL_iters',AL_iters,'verbose','yes');                   
x_w_name = strcat(mypath,'\wsu.dat');                
save(x_w_name,'x_w');
clear x_w
toc

% JSpBLRU
Y = image;
A = big_ed;
blocknum = 3;
numpat = blocknum *ones(1, fix(size(Y,2)/blocknum));
numpat(end) = numpat(end) + mod(size(Y,2),blocknum);
if sum(numpat)~= size(Y,2)
    error('block partition does not match!')
end
K0      = length(numpat);  % number of blocks
maxiter = 300;
parameter.verbose = 0;
parameter.K0      = K0;
parameter.MaxIter = maxiter;
parameter.numpat  = numpat;
parameter.mu = 0.001;
parameter.gamma = 0.001;
parameter.tau = 0.001;
parameter.imgsize = [200,200];
tic 
X_jsp = JSpBLRU(Y,A,parameter);
toc;
X_jsp_name = strcat(mypath,'\JSpBLRU.dat');                
save(X_jsp_name,'X_jsp');
clear X_jsp parameter


% BiJSpLRU
Y = image;
A = big_ed;
blocknum = 3;
numpat = blocknum *ones(1, fix(size(Y,2)/blocknum));
numpat(end) = numpat(end) + mod(size(Y,2),blocknum);
if sum(numpat)~= size(Y,2)
    error('block partition does not match!')
end
K0      = length(numpat);  % number of blocks
maxiter = 300;
parameter.verbose = 0;
parameter.K0      = K0;
parameter.MaxIter = maxiter;
parameter.numpat  = numpat;
parameter.mu = 0.001;
parameter.gamma = 0.001;
parameter.tau = 0.001;
parameter.imgsize = [200,200];
tic 
X_bijsp   = BiJSpLRU(Y,A,parameter);
toc;
X_bijsp_name = strcat(mypath,'\BiJSpLRU.dat');                
save(X_bijsp_name,'X_bijsp');
clear parameter

% LSU 
Y = image;
n_pixel = size(Y,2);
A = big_ed;
x_2 = X_bijsp;
verbose_hysime='on';
noise_type_hysime = 'additive';        
[w_hysime Rn_hysime] = estNoise(Y,noise_type_hysime,verbose_hysime);
[kf,Ek,E_hysime]=hysime(Y,w_hysime,Rn_hysime,1e-5,verbose_hysime); 
parameter.imgsize = [200,200]; 
parameter.p = kf;
parameter.alpha = 1e-2;
parameter.beta = 0.1;
parameter.iter = 30;
tic 
[x1,x_lsu] = runLSUCuprite(Y,A,x_2,parameter);
toc
x_lsu_name = strcat(mypath,'\LSU.dat');                
save(x_lsu_name,'x_lsu');



