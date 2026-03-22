function [x,obj_list,res_p_list,res_d_list] = sunsal_w(M,y,w,varargin)

%%
%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
% if (nargin-length(varargin)) ~= 2
%    error('Wrong number of required parameters');
% end
% mixing matrixsize
[LM,p] = size(M);
% data set size
[L,N] = size(y);
if (LM ~= L)
    error('mixing matrix M and data set y are inconsistent');
end
% if (L<p)
%     error('Insufficient number of columns in y');
% end


%%
%--------------------------------------------------------------
% Set the defaults for the optional parameters
%--------------------------------------------------------------
% maximum number of AL iteration
AL_iters = 1000;
% regularizatio parameter
lambda = 1;
% display only sunsal warnings
verbose = 'off';
% Positivity constraint
positivity = 'no';
% Sum-to-one constraint
addone = 'no';
% tolerance for the primal and dual residues
tol = 1e-4;
% initialization
x0 = 0;

%%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------


%--------------------------------------------------------------
% 阅读可选参数
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'AL_ITERS'
                AL_iters = round(varargin{i+1});
                if (AL_iters < 0 )
                    error('AL_iters must a positive integer');
                end
            case 'LAMBDA'
                lambda = varargin{i+1};
                if (sum(sum(lambda < 0)) >  0 )
                    error('lambda must be positive');
                end
            case 'POSITIVITY'
                positivity =  varargin{i+1};
            case 'ADDONE'
                addone = varargin{i+1};
            case 'TOL'
                tol = varargin{i+1};
            case 'VERBOSE'
                verbose = varargin{i+1};
            case 'X0'
                x0 = varargin{i+1};
                if (size(x0,1) ~= p) | (size(x0,2) ~= N)
                    error('initial X is  inconsistent with M or Y');
                end
            case 'X_SOL'
                X_sol = varargin{i+1};                
            case 'CONV_THRE'
                conv_thre = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end;
    end;
end

%---------------------------------------------
% 如果lambda是标量，则将其转换为向量
%---------------------------------------------
lambda1 = lambda;
Nlambda = size(lambda);
if Nlambda == 1
    % same lambda for all pixels
    lambda = lambda*ones(p,N);
elseif Nlambda ~= N
    error('Lambda size is inconsistent with the size of the data set');
else
    %each pixel has its own lambda
    lambda = repmat(lambda(:)',p,1);
end

% compute mean norm
norm_m = sqrt(mean(M(:).^2))*(25+p)/p;
% rescale M and Y and lambda
% M = M/norm_m;
% y = y/norm_m;
% lambda = lambda/norm_m^2;

%%
%---------------------------------------------
%  常量和初始化
%---------------------------------------------

mu_AL = 0.001;
mu = 0.001;
% mu = 10*mean(lambda(:)) + mu_AL;

F = M'*M + 2*eye(p);
IF = pinv(F);



%%
%---------------------------------------------
%  初始化
%---------------------------------------------
U = IF*M'*y;
V1 = M*U;
V2 = U;
V3 = U;

D1 = V1*0;
D2 = V2*0;
D3 = V3*0;


%%
%---------------------------------------------
%  AL 迭代-主体
%---------------------------------------------
tol1 = sqrt(N*p)*tol;
tol2 = sqrt(N*p)*tol;
i=1;
res_p = inf;
res_d = inf;
maskz = ones(size(V2));
mu_changed = 0;

% c_norm=sqrt(p*N);
z_err_set=ones(1,10)*1e10;
obj_list = [];
res_p_list = [];
res_d_list = [];
%--------------------------------------------------------------------------
%%
% ADMM iterations
%--------------------------------------------------------------------------
while(i <= AL_iters) && ((abs (res_p) > tol1) || (abs (res_d) > tol2))
    if mod(i,10) == 1
        V10 = V1;
        V20 = V2;
        V30 = V3;
    end
    
    % update U and V
    U = IF*(M'*(V1+D1)+V2+D2+V3+D3);
    V1 = 1/(1+mu) .* (y + mu.* (M*U-D1));
    V2 = soft(U-D2, (lambda/mu).*w);
    V3 = max(U-D3, 0);
    
    % update D
    D1 = D1 - M*U +V1;
    D2 = D2 - U + V2;
    D3 = D3 - U + V3;
    
    if mod(i,10) == 1
        obj = 1/2*norm(y - M*U, 'fro')^2 + lambda1.*sum(sum(abs(w.*U)));
        res_p = norm([V1;V2;V3] - [M*U;U;U], 'fro');
        res_d = norm([V1;V2;V3] - [V10;V20;V30], 'fro');
        obj_list = [obj_list,obj];
        res_p_list = [obj_list,res_p];
        res_d_list = [obj_list,res_d];
        if strcmp(verbose,'yes')
            fprintf('i = %d, obj = %f, res_p = %f, res_d = %f\n',...
                i, obj, res_p, res_d);
        end
        
        if res_p > 10*res_d
            mu = mu *2;
            D1 = D1/2;
            D2 = D2/2;
            D3 = D3/2;
        elseif res_d > 10*res_p
            mu = mu/2;
            D1 = D1*2;
            D2 = D2*2;
            D3 = D3*2;
        end

    end
    i = i + 1;
end
    x = U;
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
