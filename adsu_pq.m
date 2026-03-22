function [x,obj_list,res_p_list,res_d_list,i_list] = adsu_pq(A,y,w1,w2,varargin)
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

%%
%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
% if (nargin-length(varargin)) ~= 2
%    error('Wrong number of required parameters');
% end
% mixing matrixsize
M = [A,A];
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
tol = 0;
% initialization
x0 = 0;

%%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------


%--------------------------------------------------------------
% 闃呰鍙?夊弬鏁?
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
            case 'LAMBDA1'
                lambda1 = varargin{i+1};
                if (sum(sum(lambda1 < 0)) >  0 )
                    error('lambda1 must be positive');
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
% 濡傛灉lambda鏄爣閲忥紝鍒欏皢鍏惰浆鎹负鍚戦噺
%---------------------------------------------
lambda_1 = lambda1;
Nlambda1 = size(lambda1);

if Nlambda1 == 1
    % same lambda for all pixels
    lambda1 = lambda1*ones(p,N);
elseif Nlambda1 ~= N
    error('Lambda size is inconsistent with the size of the data set');
else
    %each pixel has its own lambda
    lambda1 = repmat(lambda1(:)',p,1);
end
% compute mean norm
norm_m = sqrt(mean(M(:).^2))*(25+p)/p;
% rescale M and Y and lambda
% M = M/norm_m;
% y = y/norm_m;
% lambda1 = lambda1/norm_m^2;
% lambda2 = lambda2/norm_m^2;
%%
%---------------------------------------------
%  甯搁噺鍜屽垵濮嬪寲
%---------------------------------------------

mu_AL = 0.001;
mu = 0.001;
% mu = 10*mean(lambda1(:)) + mu_AL;
B = [eye(p/2),eye(p/2)];

F = M'*M  + 2*eye(p);
IF = pinv(F);
% 鍏?1鐭╅樀
w_1 = ones(p/2, N);
IF_A = pinv(A'*A + 2*eye(p/2));
%%
%---------------------------------------------
%  鍒濆鍖?
%---------------------------------------------
phi = 0.5*ones(p/2,N);
U = [phi.*(w1.*(IF_A*A'*y));(w_1-phi).*w2.*(IF_A*A'*y)];
V1 = M*U;
V2 = U;
V3 = U;

D1 = V1*0;
D2 = V2*0;
D3 = V3*0;
C = [w1;w2];
%%
%---------------------------------------------
%  AL 杩唬-涓讳綋
%---------------------------------------------
tol1 = sqrt(N*p)*tol;
tol2 = sqrt(N*p)*tol;
i=1;
res_p = inf;
res_d = inf;
maskz = ones(size(V2));
mu_changed = 0;

obj_list = [];
res_p_list = [];
res_d_list = [];
i_list = [];

% c_norm=sqrt(p*N);
z_err_set=ones(1,10)*1e10;
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
    V2 = soft(U-D2, (lambda1/mu).*C);
    V3 = max(U-D3, 0);
    
    % update D
    D1 = D1 - M*U +V1;
    D2 = D2 - U + V2;
    D3 = D3 - U + V3;
    
    if mod(i,10) == 1
        obj = 1/2*norm(y - M*U, 'fro')^2 + lambda_1.*sum(sum(abs(C.*U)));
        res_p = norm([V1;V2;V3] - [M*U;U;U], 'fro');
        res_d = norm([V1;V2;V3] - [V10;V20;V30], 'fro');
        obj_list(end+1) = obj;
        res_p_list(end+1) = res_p;
        res_d_list(end+1) = res_d;
        i_list(end+1) = i;
        if strcmp(verbose,'yes')
            fprintf('i = %d, obj = %f, res_p = %f, res_d = %f\n',i, obj, res_p, res_d);
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
    P = U(1:p/2,:);
    Q = U(1+p/2:p,:);
    P(P<0)=0;
    Q(Q<0)=0;
    x =  P + Q;
    idx = x==0;
    phi = P./x;
    phi(idx)= 0;
end
