function [X] = BMSPI(A, Y, prior, parameter)
[L,N] = size(Y); 
[L,M] = size(A);
if size(A, 1) ~= L
    error(['The sizes of hyperspectral data matrix Y and spectral ',...
        'library matrix A are inconsistent!']);
end
% if isempty(prior) == 0 && (min(size(prior)) ~= 1 || ...
%         sum(ismember(prior, 1:M)) ~= length(prior))
%     error(['The input prior should be a vector whose elements are ',... 
%         'included in the indices of the spectral library!']);
% end
% 参数设置
if isfield(parameter, 'maxiter') % 迭代次数
    MaxIter = parameter.maxiter;
else
    MaxIter = 1000;
end
if isfield(parameter, 'epsilon') % 容差
    epsilon = parameter.epsilon;
else
    epsilon = 1e-5;
end
if isfield(parameter, 'mu') 
    mu = parameter.mu;
else
    mu = 0.1;
end
if isfield(parameter, 'lambda_x') 
    lambda_x = parameter.lambda_x;
else
    error('The parameter lambda_x is missing!');
end
if isfield(parameter, 'lambda_s0') 
    lambda_s0 = parameter.lambda_s0;
else
    error('The parameter lambda_s0 is missing!');
end
if isfield(parameter, 'verbose')
    verbose = parameter.verbose;
else
    verbose = 1;
end
%%
%---------------------------------------------
%  Initializations
%
lambda_s =1e-6; % 初始化为1e-4
H = eye(M);
for i = 1:length(prior)
    H(i, i) = prior(i);
end
% 使用HySime获得D
verbose_hysime='on';
noise_type_hysime = 'additive';        

[w_hysime Rn_hysime] = estNoise(Y,noise_type_hysime,verbose_hysime);
[kf,Ek,E_hysime]=hysime(Y,w_hysime,Rn_hysime,1e-5,verbose_hysime);  
D = sqrt(Rn_hysime);
D = 1./D;
D(D==inf) =0;
I = eye(M);
IF = ((D*A)'*(D*A)+mu*(I+H))^-1;
U = IF*((D*A)'*(D*Y));
V1 = Y-A*U;
V2 = H*U;
V3 = U;
D1 = V1*0;
D2 = V2*0;
D3 = V3*0;
S = RH(V1-D1,2*lambda_s/mu);

%error tolerance
tol = sqrt((3*M + L)/2*N/2)*epsilon;
%current iteration number
i = 1;
%primal residual 
res_p = inf;
%dual residual
res_d = inf;
mu_changed = 0;
%%
%---------------------------------------------
%  ADMM iterations
%---------------------------------------------
while (i <= MaxIter) && ((abs(res_p) > tol) || (abs(res_d) > tol))
    if mod(i,10) == 1
        V10 = V1;
        V20 = V2;
        V30 = V3;
        S_L20 = sqrt(sum(S.^2,2));
        S_L20 = sum(S_L20>0);
    end
    
    % update U and V
    
    U = IF * ((D*A)'*D*(Y-V1) + mu*(V2-D2+V3-D3));
    S = RH(V1-D1,2*lambda_s/mu);
    V1 = ((D'*D+mu*eye(size(D'*D)))^-1) * (D'*D*(Y-A*U)+ mu*(S+D1));
    V2 = vector_soft_row(H*U+D2, lambda_x/mu);
    V3 = max(U+D3,0);
    
    % update D
    D1 = D1 - (V1 - S);
    D2 = D2 - (V2 - H*U);
    D3 = D3 - (V3 - U);
    
    if mod(i, 10) == 1
        %primal residual
        res_p = norm([V1; V2; V3] - [S; H*U; U], 'fro');
        %dual residual
        res_d = norm([V1; V2; V3] - [V10; V20; V30], 'fro');
        res_S = sqrt(sum(S.^2,2));
        res_S = sum(res_S>0);
        if res_S == S_L20 && lambda_s*2*N <= lambda_s0
            lambda_s = 2 * lambda_s;
        else 
            lambda_s = lambda_s;
        end
        obj = 0.5*norm(D*(Y-A*U-S)) + lambda_x.*sum(sum(abs(H*U))) + lambda_s*res_S;
        if verbose
            fprintf('i = %d, obj = %f,res_p = %f, res_d = %f, mu = %f\n',...
            i,obj, res_p, res_d, mu);
        end
        if res_p > 10 * res_d
            mu = mu*2;
            D1 = D1/2;
            D2 = D2/2;
            D3 = D3/2;
            mu_changed = 1;
        elseif res_d > 10*  res_p
            mu = mu/2;
            D1 = D1*2;
            D2 = D2*2;
            D3 = D3*2;
            mu_changed = 1;
        end
        if mu_changed
            IF = ((D*A)'*(D*A)+mu*(I+H))^-1;
            mu_changed = 0;
        end
    end
    i =i +1;
end
X = U.*(U>=0);
end