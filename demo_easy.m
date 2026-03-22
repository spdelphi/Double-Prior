
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



%  30db 

%% test_prior
it = 1000;
idx_other = 1:m_lib;
idx_other(idx) = [];
othersz = length(idx_other);



Image = Image_30DB;
lambda = 0.001;
p=0;
parameter.lambda_x = 0.5;
parameter.lambda_s0 = 10;

% only for C1, more cases can be seen in demo_DC1.m
%1100000
w = ones(m_lib,N);
w(idx(1:2),:) = p;
prior = w; 
w = normWeight(w);
label = '30db_C1';
tic
[x_ad,phi,P,Q] = adsu_pq(M,Image,w,ones(m_lib,N)./sqrt(m_lib),'POSITIVITY','yes','lambda1',lambda,'AL_ITERS',it);
toc


%%
%---------------------------------------------
%  evaluate 
%---------------------------------------------

document = mypath;


C_adsu_RMSE = getRMSE(x-x_ad);
C_adsu_SRE = getSRE(x,x_ad);

%%
%---------------------------------------------
%  fig_DC1
%---------------------------------------------
sum_re = sum(abs(x),2)';
   
% ADSU
sum_ad = sum(abs(x_ad),2)';
adsu_map = createstem(sum_re,sum_ad,'');



