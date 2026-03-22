
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
%---------------------------------------------

close all;
clear global;
clear 

mypath = 'DC1_data\';
post = '_x';
maxAb = 0.7;
EmNum = 7;
EmNum_true = EmNum;
testEd = 1;

load([mypath 'x.mat']);
load([ mypath , 'Image_30DB.mat']);
load([ mypath , 'M.mat']);
Image_30DB = Image;
A = M;
Y = Image_30DB;
m_lib = size(M,2);
N = size(Image,2);
idx = [15,80,128,200,302,450,256];
%% 30db


%test other algorihtms, you can skip them if you do not intent to compare
%with them

% % JSpBLRU
% 
% blocknum = 3;
% numpat = blocknum *ones(1, fix(size(Y,2)/blocknum));
% numpat(end) = numpat(end) + mod(size(Y,2),blocknum);
% if sum(numpat)~= size(Y,2)
%     error('block partition does not match!')
% end
% 
% K0 = length(numpat);  % number of blocks
% maxiter = 300;
% load('.\\parameter\\JSpBLRU_30db_parameter.mat');
% parameter.verbose = 0;
% parameter.K0      = K0;
% parameter.MaxIter = maxiter;
% parameter.numpat  = numpat;
% parameter.imgsize = [50,50];
% tic
% X_jsp = JSpBLRU(Y,A,parameter);
% toc;
% X_jsp_name = strcat(mypath,'\x_30db_JSpBLRU.mat');
% save(X_jsp_name,'X_jsp');
% 
% % BiJSpLRU
% load('.\\parameter\\BiJSpLRU_30db_parameter.mat');
% parameter.verbose = 0;
% parameter.K0      = K0;
% parameter.MaxIter = maxiter;
% parameter.numpat  = numpat;
% tic
% X_bijsp   = BiJSpLRU(Y,A,parameter);
% toc;
% X_bijsp_name = strcat(mypath,'\x_30db_BiJSpLRU.mat');
% save(X_bijsp_name,'X_bijsp');
% 
% % LSU
% load('.\\parameter\\lsu_30db_parameter.mat');
% parameter.imgsize = [50,50];
% tic
% [x1,x_lsu] = runLSUCuprite(Y,A,X_bijsp,parameter);
% cpu_lsu = toc;
% X_lsu_name = strcat(mypath,'\x_30db_LSU.mat');
% save(X_lsu_name,'x_lsu');


%  30db 

%% test_prior
document = mypath;
it = 1000;
idx_other = 1:m_lib;
idx_other(idx) = [];
othersz = length(idx_other);
p=0;



% su
lambda = 0.01;
tic
[x_su] = sunsal_w(A,Y,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_su_name = strcat(mypath,'\x_30db_SU.mat');
save(x_su_name,'x_su');

Image = Image_30DB;
lambda = 0.001;
p=0;
parameter.lambda_x = 0.5;
parameter.lambda_s0 = 10;

% C1
%1100000
w = ones(m_lib,N);
w(idx(1:2),:) = p;
prior = w; 
w = normWeight(w);
label = '30db_C1'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad'); 
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');

% C2
%1100000(2-)
w = ones(m_lib,N);
w(idx(1:2),:) = p;
dismis_id = [20,60];
w(idx_other(dismis_id),:) = 2;
prior = w;
w = normWeight(w);
label = '30db_C2'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad'); 
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');

% C3
% 1111000
w = ones(m_lib,N);
w(idx(1:4),:) = p;
prior = w;
w = normWeight(w);
label = '30db_C3'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad'); 
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');

% C4
%1100000(2-)
w = ones(m_lib,N);
w(idx(1:4),:) = p;
dismis_id = [20,60];
w(idx_other(dismis_id),:) = 2;
prior = w;
w = normWeight(w);
label = '30db_C4'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad'); 
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C5
%1111111
w = ones(m_lib,N);
w(idx,:) = p;
prior = w;
w = normWeight(w);
label = '30db_C5'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');

% C6
% 1111111(20)-(60)-
w = ones(m_lib,N);
w(idx,:) = p;
dismis_id = [20,60];
w(idx_other(dismis_id),:) = 2;
prior = w; 
w = normWeight(w);
label = '30db_C6'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C7
%0000000(50)(51)
w = ones(m_lib,N);
% w(idx,:) = p;
mis_id = [50,51];
w(idx_other(mis_id),:) = p;
prior = w; 
w = normWeight(w);
label = '30db_C7'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C8
% 0000000(50)(51)(100)(368)
w = ones(m_lib,N);
% w(idx,:) = p;
mis_id = [50,51,100,368];
w(idx_other(mis_id),:) = p;
prior = w; 
w = normWeight(w);
label = '30db_C8'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C9
% 1100000(50)(51)
w = ones(m_lib,N);
w(idx(1:2),:) = p;
mis_id = [50,51];
w(idx_other(mis_id),:) = p;
prior = w; 
w = normWeight(w);
label = '30db_C9'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');

% C10
% 1100000(50)(51)(100)(368)
w = ones(m_lib,N);
w(idx(1:2),:) = p;
mis_id = [50,51,100,368];
w(idx_other(mis_id),:) = p;
prior = w; 
w = normWeight(w);
label = '30db_C10'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');

% C11
% 1111000(50)(51)
w = ones(m_lib,N);
w(idx(1:4),:) = p;
mis_id = [50,51];
w(idx_other(mis_id),:) = p;
prior = w;
w = normWeight(w);
label = '30db_C11'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');

% C12
% 1111000(50)(51)(100)(368)
w = ones(m_lib,N);
w(idx(1:4),:) = p;
mis_id = [50,51,100,368];
w(idx_other(mis_id),:) = p;
prior = w;
w = normWeight(w);
label = '30db_C12'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C13
%1-1-00000
w = ones(m_lib,N);
w(idx(1:2),:) = 2;
prior = w; 
w = normWeight(w);
label = '30db_C13'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C14
% 1-1-1-1-000
w = ones(m_lib,N);
w(idx(1:4),:) = 2;
prior = w; 
w = normWeight(w);
label = '30db_C14'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C15
% 111-1-000
w = ones(m_lib,N);
w(idx(1:2),:) = p;
w(idx(3:4),:) = 2;
prior = w;
w = normWeight(w);
label = '30db_C15'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C16
% 111-1-000(20)-(60)-
w = ones(m_lib,N);
w(idx(1:2),:) = p;
w(idx(3:4),:) = 2;
dismis_id = [20,60];
w(idx_other(dismis_id),:) = 2;
prior = w; 
w = normWeight(w);
label = '30db_C16'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C17
% 111-1-000(20)-(60)-(100)-(368)-
w = ones(m_lib,N);
w(idx(1:2),:) = p;
w(idx(3:4),:) = 2;
dismis_id = [20,60,100,368];
w(idx_other(dismis_id),:) = 2;
prior = w; 
w = normWeight(w);
label = '30db_C17'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C18
% 1-1-1-1-1-1-1-
w = ones(m_lib,N);
w(idx,:) = 2;
% dismis_id = [20,60,100,368];
% w(idx_other(dismis_id),:) = 2;
prior = w; 
w = normWeight(w);
label = '30db_C18'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C19
% 1-1-00000(50)(51)
w = ones(m_lib,N);
w(idx(1:2),:) = 2;
mis_id = [50,51]
w(idx_other(mis_id),:) = 0;
prior = w; 
w = normWeight(w);
label = '30db_C19'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C20
%1-1-1-1-000(50)(51)
w = ones(m_lib,N);
w(idx(1:4),:) = 2;
mis_id = [50,51]
w(idx_other(mis_id),:) = p;
prior = w; 
w = normWeight(w);
label = '30db_C20'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C21
%111-1-000(50)(51)
w = ones(m_lib,N);
w(idx(1:2),:) = p;
w(idx(3:4),:) = 2;
% dismis_id = [20,60,100,368];
% w(idx_other(dismis_id),:) = 2;
mis_id = [50,51]
w(idx_other(mis_id),:) = p;
prior = w; 
w = normWeight(w);
label = '30db_C21'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


% C22
% 111-1-000(20)-(50)(51)(60)-
w = ones(m_lib,N);
w(idx(1:2),:) = p;
w(idx(3:4),:) = 2;
dismis_id = [20,60];
w(idx_other(dismis_id),:) = 2;
mis_id = [50,51]
w(idx_other(mis_id),:) = p;
prior = w;
w = normWeight(w);
label = '30db_C22'
tic
[x_w] = sunsal_w(M,Image,w,'POSITIVITY','yes','lambda',lambda,'AL_ITERS',it);
toc
x_w_name = strcat(document,'su_w_',label,'.mat');
save(x_w_name,'x_w');  
tic;
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc
x_ad_name = strcat(document,'adsu_',label,'.mat');
save(x_ad_name,'x_ad');
tic
x_BMSPI = BMSPI(A,Y,prior,parameter);
elapsedtime = toc;
x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
save(x_BMSPI_name,'x_BMSPI');


%%
%---------------------------------------------
%  evaluate 
%---------------------------------------------

document = mypath;
data_id = '30db';
p = 0;
C = {};
% 生成label
for i = 1:22
    c_name = strcat(data_id,'_C',num2str(i) );
    C{end+1} =  c_name;
end


C_su_w_RMSE = [];
C_su_w_SRE = [];
C_adsu_RMSE = [];
C_adsu_SRE = [];
C_BMSPI_RMSE = [];
C_BMSPI_SRE = [];
x_threshold = 0+1e-60;
abundance_threshold = 0;

for i = 1:length(C)
    label = strjoin(C(i));
    x_w_name = strcat(document,'su_w_',label,'.mat');
    x_ad_name = strcat(document,'adsu_',label,'.mat');
    x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
    load(x_w_name);
    load(x_ad_name);
    load(x_BMSPI_name); 
    C_su_w_RMSE = [C_su_w_RMSE,getRMSE(x-x_w)];
    C_su_w_SRE = [C_su_w_SRE,getSRE(x,x_w)];
    C_adsu_RMSE = [C_adsu_RMSE,getRMSE(x-x_ad)];
    C_adsu_SRE = [C_adsu_SRE,getSRE(x,x_ad)];
    C_BMSPI_RMSE = [C_BMSPI_RMSE,getRMSE(x-x_BMSPI)];
    C_BMSPI_SRE = [C_BMSPI_SRE,getSRE(x,x_BMSPI)];
end

% SU
% 20DB
% x_su_name = strcat(mypath,'\x_20db_SU.mat');
% load(x_su_name);
% c_rmse_20db = getRMSE(x-x_su);
% c_sre_20db = getSRE(x,x_su);
x_su_name = strcat(mypath,'\x_30db_SU.mat');
load(x_su_name);
c_rmse_30db = getRMSE(x-x_su);
c_sre_30db = getSRE(x,x_su);
% x_su_name = strcat(mypath,'\x_40db_SU.mat');
% load(x_su_name);
% c_rmse_40db = getRMSE(x-x_su);
% c_sre_40db = getSRE(x,x_su);


% JSpBLRU
% 20DB
% X_jsp_name = strcat(mypath,'\x_20db_JSpBLRU.mat');
% load(X_jsp_name);
% c_rmse_20db = getRMSE(x-X_jsp);
% c_sre_20db = getSRE(x,X_jsp);
X_jsp_name = strcat(mypath,'\x_30db_JSpBLRU.mat');
load(X_jsp_name);
c_rmse_30db = getRMSE(x-X_jsp);
c_sre_30db = getSRE(x,X_jsp);
% X_jsp_name = strcat(mypath,'\x_40db_JSpBLRU.mat');
% load(X_jsp_name);
% c_rmse_40db = getRMSE(x-X_jsp);
% c_sre_40db = getSRE(x,X_jsp);



%BiJSpLRU
% % 20DB
% X_bijsp_name = strcat(mypath,'\x_20db_BiJSpLRU.mat');
% load(X_bijsp_name);
% c_rmse_20db = getRMSE(x-X_bijsp);
% c_sre_20db = getSRE(x,X_bijsp);
X_bijsp_name = strcat(mypath,'\x_30db_BiJSpLRU.mat');
load(X_bijsp_name);
c_rmse_30db = getRMSE(x-X_bijsp);
c_sre_30db = getSRE(x,X_bijsp);
% X_bijsp_name = strcat(mypath,'\x_40db_BiJSpLRU.mat');
% load(X_bijsp_name);
% c_rmse_40db = getRMSE(x-X_bijsp);
% c_sre_40db = getSRE(x,X_bijsp);



%LSU
% 20DB
% X_lsu_name = strcat(mypath,'\x_20db_LSU.mat');
% load(X_lsu_name);
% c_rmse_20db = getRMSE(x-x_lsu);
% c_sre_20db = getSRE(x,x_lsu);
X_lsu_name = strcat(mypath,'\x_30db_LSU.mat');
load(X_lsu_name);
c_rmse_30db = getRMSE(x-x_lsu);
c_sre_30db = getSRE(x,x_lsu);
% X_lsu_name = strcat(mypath,'\x_40db_LSU.mat');
% load(X_lsu_name);
% c_rmse_40db = getRMSE(x-x_lsu);
% c_sre_40db = getSRE(x,x_lsu);


%SUnCNN
% 20DB
% X_SUnCNN_name = strcat(mypath,'\x_20db_SUnCNN.mat');
% load(X_SUnCNN_name);
% x_SUnCNN = reshape(out_avg_np_20db_21,50*50,498);
% x_SUnCNN = x_SUnCNN';
% c_rmse_20db = getRMSE(x-x_SUnCNN);
% c_sre_20db = getSRE(x,x_SUnCNN);
X_SUnCNN_name = strcat(mypath,'\x_30db_SUnCNN.mat');
load(X_SUnCNN_name);
x_SUnCNN = reshape(out_avg_np_30db_21,50*50,498);
x_SUnCNN = x_SUnCNN';
c_rmse_30db = getRMSE(x-x_SUnCNN);
c_sre_30db = getSRE(x,x_SUnCNN);
% X_SUnCNN_name = strcat(mypath,'\x_40db_SUnCNN.mat');
% load(X_SUnCNN_name);
% x_SUnCNN = reshape(out_avg_np_40db_21,50*50,498);
% x_SUnCNN = x_SUnCNN';
% c_rmse_40db = getRMSE(x-x_SUnCNN);
% c_sre_40db = getSRE(x,x_SUnCNN);


%%
%---------------------------------------------
%  fig_DC1
%---------------------------------------------
sum_re = sum(abs(x),2)';

% su
x_su_name = strcat(mypath,'\x_30db_SU.mat');
load(x_su_name);
sum_su = sum(abs(x_su),2)';
su_map = createstem(sum_re,sum_su,'SU');
su_name = strcat(document,'\SU');
print(su_map,'-djpeg',su_name,'-r1800');

% jsp
X_jsp_name = strcat(mypath,'\x_30db_JSpBLRU.mat');
load(X_jsp_name);
sum_jsp = sum(abs(X_jsp),2)';
jsp_map = createstem(sum_re,sum_jsp,'JSpBLRU');
jsp_name = strcat(document,'\JSpBLRU');
print(jsp_map,'-djpeg',jsp_name,'-r1800');

% BiJSpLRU
X_bijsp_name = strcat(mypath,'\x_30db_BiJSpLRU.mat');
load(X_bijsp_name);
sum_bijsp = sum(abs(X_bijsp),2)';
bijsp_map = createstem(sum_re,sum_bijsp,'BiJSpLRU');
bijsp_name = strcat(document,'\BiJSpLRU');
print(bijsp_map,'-djpeg',bijsp_name,'-r1800');
% 
% LSU
X_lsu_name = strcat(mypath,'\x_30db_LSU.mat');
load(X_lsu_name);
sum_lsu = sum(abs(x_lsu),2)';
lsu_map = createstem(sum_re,sum_lsu,'LSU');
lsu_name = strcat(document,'\LSU');
print(lsu_map,'-djpeg',lsu_name,'-r1800');
close all



% WSU、ADSU 
document = mypath;
data_id = '30db';
p = 0;
C = {};
% 生成label
for i = 1:22
    c_name = strcat(data_id,'_C',num2str(i) );
    C{end+1} =  c_name;
end



for i = 1:length(c_label)
    label = strcat(data_id,'_',strjoin(c_label(i)) );
    x_w_name = strcat(document,'su_w_',label,'.mat')
    x_ad_name = strcat(document,'adsu_',label,'.mat');
    x_BMSPI_name = strcat(document,'x_BMSPI_',label,'.mat');
    load(x_w_name);
    load(x_ad_name);
    load(x_BMSPI_name); 
    x_w_name = strcat(document1,'su_w_',label,'p_',num2str(p),'.mat');

    % WSU
    sum_w = sum(abs(x_w),2)';
    wsu_C_name = strcat('WSU(case:C',num2str(i),')');
    wsu_map = createstem(sum_re,sum_w,wsu_C_name);
    wsu_name = strcat(document,'\WSU_C',num2str(i));
    print(wsu_map,'-djpeg',wsu_name,'-r1800');
    
    % ADSU
    sum_ad = sum(abs(x_ad),2)';
    adsu_C_name = strcat('ADSU(case:C',num2str(i),')');
    adsu_map = createstem(sum_re,sum_ad,adsu_C_name);
    adsu_name = strcat(document_save_all,'\ADSU_C',num2str(i));
    print(adsu_map,'-djpeg',adsu_name,'-r1800');
    close all

    % BMSPI
    sum_BMSPI = sum(abs(x_BMSPI),2)';
    BMSPI_C_name = strcat('BMSPI(case:C',num2str(i),')');
    BMSPI_map = createstem(sum_re,sum_BMSPI,BMSPI_C_name);
    BMSPI_name = strcat(document,'\BMSPI_C',num2str(i));
    print(BMSPI_map,'-djpeg',BMSPI_name,'-r1800');
    close all
end

% 
% True
re_C_name = 'True Abundance'
sum_re1 = zeros(size(sum_re));
re_map = createstem(sum_re,sum_re1,re_C_name);
re_name = strcat(document,'\True Abundance');
print(re_map,'-djpeg',re_name,'-r1800');
close all 

% SUnCNN
X_SUnCNN_name = strcat(mypath,'\x_30db_SUnCNN.mat');
load(X_lsu_name);
x_SUnCNN = reshape(out_avg_np_20db_21,50*50,498);
x_SUnCNN = x_SUnCNN';
sum_SUnCNN = sum(abs(x_SUnCNN),2)';
SUnCNN_map = createstem(sum_re,sum_SUnCNN,'SUnCNN');
SUnCNN_name = strcat(document,'\SUnCNN');
print(SUnCNN_map,'-djpeg',SUnCNN_name,'-r1800');
close all