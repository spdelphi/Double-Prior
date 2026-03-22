close all;
clear global;
clear 

mypath = 'D:\work\TGRS\data\BMSPI_lambda';
post = '_x';
maxAb = 0.7;
EmNum = 7;
EmNum_true = EmNum;
testEd = 1;

if true%~exist([ mypath 'simu',post,'.mat'])
    index = 1:EmNum;%indexall(1:EmNum);
    lib_name = strcat(mypath,'\USGS_1995_Library.mat');
    load(lib_name);
    [dummy indexM] = sort(datalib(:,1));
    M =  datalib(indexM,4:end);
    names = names(4:end,:);
    M = normSpectral(M);
    m_lib = size(M,2);
    m = 9;
    idx = [15,80,128,200,302,450,256];%randi([1,m_lib],m,1);%
%     M = M(:,randi(500,EmNum,1));%1:EmNum);
    Em = M(:,idx);
    maxAbs = maxAb;
    N = 1000;

    [Image,~,x_sq]=Simulation_data(2,Em,0.5,m,N,maxAbs,[]);%use USGS
    x = zeros(m_lib,N);
    x(idx,:) = x_sq;
    linear = Image;
    [bands,N] = size(Image);            
    oldImg = Image;          
    [Image_20DB, noise,sigma ] = SimuNoise( oldImg , 'SNR' ,20 );
    [Image_30DB, noise,sigma ] = SimuNoise( oldImg , 'SNR' ,30 );
    [Image_40DB, noise,sigma ] = SimuNoise( oldImg , 'SNR' ,40 );
end


%% 测试参数
% 参数设置 required parameters
parameter.epsilon = 1e-6;%1e-5;
parameter.maxiter = 1000;
parameter.mu = 0.1;
parameter.verbose = 1;
lambda_x_list = [1,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001,5e-5,1e-5];
lambda_s0_list = [20,10,1,1e-1,1e-2,1e-3,1e-4,1e-5];
idx_other = 1:m_lib;
idx_other(idx) = [];
othersz = length(idx_other);


% 20db
Xtrue =x;
A = M;
Y = Image_20DB;
no_l1 = length(lambda_x_list);
no_l2 = length(lambda_s0_list);
% case(1)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(3:4),:) = 0;
prior = w;
label = '20db_case(1)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');


% case(2)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(2:6),:) = 0;
prior = w;
label = '20db_case(2)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

% case(3)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(1:2),:) = 0;
w([23,45],:) = 0;
prior = w;
label = '20db_case(3)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

%case(4)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(1:4),:) = 0;
w([10,15,180,401],:) = 0;
prior = w;
label = '20db_case(4)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');


% case(5)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(1:2),:) = 0;
w([23,45,76,89,123,234,345],:) = 2;
prior = w;
label = '20db_case(5)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

% case(6)
SRE_temp_jsp = zeros(no_l1,no_l2);
w = ones(m_lib,1);
w(idx(1:2),:) = 0;
w([idx(3:end)],:) = 2;
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
prior = w;
label = '20db_case(6)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

% case(7)
SRE_temp_jsp = zeros(no_l1,no_l2);
w = ones(m_lib,1);
w(idx(1:2),:) = 0;
w([idx(3:4)],:) = 2;
w([76,89],:) = 2;
w([234,345,401],:) = 0;
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
prior = w;
label = '20db_case(7)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

% 
%% 30db
Xtrue =x;
A = M;
Y = Image_30DB;
no_l1 = length(lambda_x_list);
no_l2 = length(lambda_s0_list);
% case(1)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(3:4),:) = 0;
prior = w;
label = '30db_case(1)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');


% case(2)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(2:6),:) = 0;
prior = w;
label = '30db_case(2)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

% case(3)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(1:2),:) = 0;
w([23,45],:) = 0;
prior = w;
label = '30db_case(3)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

%case(4)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(1:4),:) = 0;
w([10,15,180,401],:) = 0;
prior = w;
label = '30db_case(4)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');


% case(5)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(1:2),:) = 0;
w([23,45,76,89,123,234,345],:) = 2;
prior = w;
label = '30db_case(5)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

% case(6)
SRE_temp_jsp = zeros(no_l1,no_l2);
w = ones(m_lib,1);
w(idx(1:2),:) = 0;
w([idx(3:end)],:) = 2;
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
prior = w;
label = '30db_case(6)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

% case(7)
SRE_temp_jsp = zeros(no_l1,no_l2);
w = ones(m_lib,1);
w(idx(1:2),:) = 0;
w([idx(3:4)],:) = 2;
w([76,89],:) = 2;
w([234,345,401],:) = 0;
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
prior = w;
label = '30db_case(7)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

%% 40db
Xtrue =x;
A = M;
Y = Image_40DB;
no_l1 = length(lambda_x_list);
no_l2 = length(lambda_s0_list);
% case(1)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(3:4),:) = 0;
prior = w;
label = '40db_case(1)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');


% case(2)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(2:6),:) = 0;
prior = w;
label = '40db_case(2)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

% case(3)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(1:2),:) = 0;
w([23,45],:) = 0;
prior = w;
label = '40db_case(3)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

%case(4)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(1:4),:) = 0;
w([10,15,180,401],:) = 0;
prior = w;
label = '40db_case(4)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');


% case(5)
SRE_temp_jsp = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
w = ones(m_lib,1);
w(idx(1:2),:) = 0;
w([23,45,76,89,123,234,345],:) = 2;
prior = w;
label = '40db_case(5)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

% case(6)
SRE_temp_jsp = zeros(no_l1,no_l2);
w = ones(m_lib,1);
w(idx(1:2),:) = 0;
w([idx(3:end)],:) = 2;
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
prior = w;
label = '40db_case(6)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');

% case(7)
SRE_temp_jsp = zeros(no_l1,no_l2);
w = ones(m_lib,1);
w(idx(1:2),:) = 0;
w([idx(3:4)],:) = 2;
w([76,89],:) = 2;
w([234,345,401],:) = 0;
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
prior = w;
label = '40db_case(7)';
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_jsp(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
save(SRE_name,'SRE_temp_jsp');



%% DC2_test_lambda
clc
close all;
clear global;
clear 
mypath = 'D:\work\TGRS\data\BMSPI_lambda';
post = '_x';
maxAb = 0.7;
EmNum = 5;
EmNum_true = EmNum;
testEd = 1;

if true%~exist([ mypath 'simu',post,'.mat'])
    index = 1:EmNum;%indexall(1:EmNum);
    lib_name = strcat(mypath,'\USGS_1995_Library.mat');
    [dummy indexM] = sort(datalib(:,1));
    M =  datalib(indexM,4:end);
    names = names(4:end,:);
    M = normSpectral(M);
    m_lib = size(M,2);
    m = 5;
    idx = [15,80,128,200,302];%randi([1,m_lib],m,1);%
    Em = M(:,idx);
    maxAbs = maxAb;
    N = 100;
        [Image,~,x_sq]=Simulation_data(2,Em,0.5,m,N,maxAbs,[]);%use USGS
    x = zeros(m_lib,N);
    x(idx,:) = x_sq;
    linear = Image;
    [bands,N] = size(Image);            
    oldImg = Image;          
    [Image_20DB, noise,sigma ] = SimuNoise( oldImg , 'SNR' ,20 );
    [Image_30DB, noise,sigma ] = SimuNoise( oldImg , 'SNR' ,30 );
    [Image_40DB, noise,sigma ] = SimuNoise( oldImg , 'SNR' ,40 );
end


% required parameters
parameter.epsilon = 1e-6;%1e-5;
parameter.maxiter = 1000;
parameter.mu = 0.1;
parameter.verbose = 1;
lambda_x_list = [1,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001,5e-5,1e-5];
lambda_s0_list = [100,50,20,10,1,1e-1,1e-2,1e-3,1e-4,1e-5];
idx_other = 1:m_lib;
idx_other(idx) = [];
othersz = length(idx_other);

Xtrue = x
A = M;
% 20db
Y = Image_20DB;
w = ones(m_lib,1);
w(idx(3:4),:) = 0;
w(idx(1),:) = 2;
w([23,45],:) = 0;
prior = w;
no_l1 = length(lambda_x_list);
no_l2 = length(lambda_s0_list);
SRE_temp_bmspi = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_bmspi(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_DC2_20DB.mat');
save(SRE_name,'SRE_temp_bmspi');

% 30db
Y = Image_20DB;
w = ones(m_lib,N);
w(idx(3:4),:) = 0;
w(idx(1),:) = 2;
w([23,45],:) = 0;
prior = w;
no_l1 = length(lambda_x_list);
no_l2 = length(lambda_s0_list);
SRE_temp_bmspi = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_bmspi(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_DC2_30DB.mat');
save(SRE_name,'SRE_temp_bmspi');

% 40db
Y = Image_40DB;
w = ones(m_lib,N);
w(idx(3:4),:) = 0;
w(idx(1),:) = 2;
w([23,45],:) = 0;
prior = w;
no_l1 = length(lambda_x_list);
no_l2 = length(lambda_s0_list);
SRE_temp_bmspi = zeros(no_l1,no_l2);
kk0 = waitbar(0, 'Searching optimal parameters for BMSPI ...');
kk = 1;
for j1 = 1:no_l1
    for j2 = 1:no_l2
        waitbar(kk/(no_l1*no_l2), kk0)
        kk = kk+1;
        parameter.lambda_x  = lambda_x_list(j1);
        parameter.lambda_s0    = lambda_s0_list(j2);
        tic 
        x_BMSPI = BMSPI(A,Y,prior,parameter);
        toc
        SRE_temp_bmspi(j1,j2) = getSRE(Xtrue,x_BMSPI);
    end
end
close(kk0)
SRE_name = strcat(mypath,'\BMSPI_SRE_DC2_40DB.mat');
save(SRE_name,'SRE_temp_bmspi');


%% result
mypath = 'D:\work\TGRS\data\BMSPI_lambda';
% 20db
label = '20db_case(1)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db1 = SRE_temp_jsp;
label = '20db_case(2)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db2 = SRE_temp_jsp;
label = '20db_case(3)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db3 = SRE_temp_jsp;
label = '20db_case(4)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db4 = SRE_temp_jsp;
label = '20db_case(5)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db5 = SRE_temp_jsp;
label = '20db_case(6)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db6 = SRE_temp_jsp;
label = '20db_case(7)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db7 = SRE_temp_jsp;

% 30db
label = '30db_case(1)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db1 = SRE_temp_jsp;
label = '30db_case(2)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db2 = SRE_temp_jsp;
label = '30db_case(3)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db3 = SRE_temp_jsp;
label = '30db_case(4)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db4 = SRE_temp_jsp;
label = '30db_case(5)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db5 = SRE_temp_jsp;
label = '30db_case(6)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db6 = SRE_temp_jsp;
label = '30db_case(7)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db7 = SRE_temp_jsp;

% 40db
label = '40db_case(1)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db1 = SRE_temp_jsp;
label = '40db_case(2)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db2 = SRE_temp_jsp;
label = '40db_case(3)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db3 = SRE_temp_jsp;
label = '40db_case(4)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db4 = SRE_temp_jsp;
label = '40db_case(5)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db5 = SRE_temp_jsp;
label = '40db_case(6)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db6 = SRE_temp_jsp;
label = '40db_case(7)';
SRE_name = strcat(mypath,'\BMSPI_SRE_',label,'.mat');
load(SRE_name);
SRE_temp_jsp_20db7 = SRE_temp_jsp;

% DC2
SRE_name = strcat(mypath,'\BMSPI_SRE_DC2_20DB.mat');
max(SRE_temp_bmspi,[],"all")
