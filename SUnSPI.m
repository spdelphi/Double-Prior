function X = SUnSPI(Y, A, prior, parameter)
% Sparse Unmixing of Hyperspectral Data Using Spectral a Priori 
% Information [SUnSPI]
% 
% [SUnSPI] Wei Tang, Zhenwei Shi, Ying Wu, and Changshui Zhang, "Sparse 
% Unmixing of Hyperspectral Data Using Spectral a Priori Information", 
% IEEE Transactions on Geoscience and Remote Sensing, Article in Press.
%
% -------------------------------------------------------------------
%
% Usage:
%
% X = SUnSPI(Y, A, prior, parameter)
%
% ------- Input variables -------------------------------------------
%
%  Y - hyperspectral data matrix with dimensions L(bands) x K(pixels)
%
%  A - spectral library matrix with dimentsions L(bands) x m(spectra)
%
%  prior - set of indices of spectral signatures in the library that are
%  known to exist in the hyperspectral scene
%
%  parameter.
%            * lambda_p - scalar, regularization parameter of the spectral 
%              a priori information regularizer
%            * lambda_s - scalar, regularization parameter of the l1 norm  
%              regularizer          
%            mu - scalar, initial augmented Lagrangian penalty parameter, 
%              default: 0.1
%            epsilon - scalar, scaled error tolerance, default: 1e-5
%            maxiter - scalar, maximum iteration number, default: 500
%            verbose - output the process (1) or not (0), default: 1
%
%  NOTE: PARAMETERS WITH SYMBOL * ARE NECESSARY
%
% ------- Output variables -------------------------------------------
%
% X - the estimated abundance matrix with dimensions m(spectra) x K(pixels)
%
% ---------------------------------------------------------------------
%
% Please see [SUnSPI] for more details.
%
% Please contact Wei Tang (tangwei@sa.buaa.edu.cn) to report bugs or 
% provide suggestions and discussions for the codes.
%
% ---------------------------------------------------------------------
% version: 1.0 (29-May-2014)
% ---------------------------------------------------------------------
%
% Copyright (May, 2014):       Wei Tang (tangwei@sa.buaa.edu.cn)
%                              Zhenwei Shi (shizhenwei@buaa.edu.cn)
%                              Ying Wu (yingwu@eecs.northwestern.edu)
%                              Changshui Zhang (zcs@mail.tsinghua.edu.cn)
%
% SUnSPI is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ---------------------------------------------------------------------

%%
%--------------------------------------------------------------
% load and test required parameters
%--------------------------------------------------------------
[L, K] = size(Y);
m = size(A, 2);
if size(A, 1) ~= L
    error(['The sizes of hyperspectral data matrix Y and spectral ',...
        'library matrix A are inconsistent!']);
end

if isempty(prior) == 0 && (min(size(prior)) ~= 1 || ...
        sum(ismember(prior, 1:m)) ~= length(prior))
    error(['The input prior should be a vector whose elements are ',... 
        'included in the indices of the spectral library!']);
end

if isfield(parameter, 'maxiter')
    MaxIter = parameter.maxiter;
else
    MaxIter = 500;
end

if isfield(parameter, 'epsilon')
    epsilon = parameter.epsilon;
else
    epsilon = 1e-5;
end

if isfield(parameter, 'mu')
    mu = parameter.mu;
else
    mu = 0.1;
end

if isfield(parameter, 'lambda_p')
    lambda1 = parameter.lambda_p;
else
    error('The parameter lambda_p is missing!');
end

if isfield(parameter, 'lambda_s')
    lambda2 = parameter.lambda_s;
else
    error('The parameter lambda_s is missing!');
end

if isfield(parameter, 'verbose')
    verbose = parameter.verbose;
else
    verbose = 1;
end
%%
%---------------------------------------------
%  Initializations
%---------------------------------------------
H = eye(m);
for i = 1:length(prior)
    H(prior(i), prior(i)) = 0;
end

IF = (A'*A + H + 2*eye(m))^-1;
U = IF*A'*Y;

V1 = A*U;
V2 = H*U;
V3 = U;
V4 = U;

D1 = V1*0;
D2 = V2*0;
D3 = V3*0;
D4 = V4*0;

%current iteration number
i = 1;

%primal residual 
res_p = inf;

%dual residual
res_d = inf;

%error tolerance
tol = sqrt((3*m + L)/2*K/2)*epsilon;


%%
%---------------------------------------------
%  ADMM iterations
%---------------------------------------------
while (i <= MaxIter) && ((abs(res_p) > tol) || (abs(res_d) > tol))
    if mod(i, 10) == 1
        V10 = V1;
        V20 = V2;
        V30 = V3;
        V40 = V4;
    end
    %update U and V
    U = IF*(A'*(V1 + D1) + H*(V2 + D2) + (V3 + D3) + (V4 + D4));
    V1 = 1/(1+mu)*(Y + mu*(A*U - D1));
    V2 = vecsoft(H*U - D2, lambda1/mu);
    V3 = soft(U - D3, lambda2/mu);
    V4 = max(U - D4, 0);
    %update D
    D1 = D1 - A*U + V1;
    D2 = D2 - H*U + V2;
    D3 = D3 - U + V3;
    D4 = D4 - U + V4;
    
    if mod(i, 10) == 1
        %object function
        obj = 1/2*norm(A*U - Y, 'fro')^2 + ...
            lambda1*sum(sqrt(sum((H*U).*(H*U)))) + ...
            lambda2*sum(sum(abs(U)));
        %primal residual
        res_p = norm([V1; V2; V3; V4] - [A*U; H*U; U; U], 'fro');
        %dual residual
        res_d = norm([V1; V2; V3; V4] - [V10; V20; V30; V40], 'fro');
        
        if verbose
            fprintf('i = %d, obj = %f,res_p = %f, res_d = %f, mu = %f\n',...
            i, obj, res_p, res_d, mu);
        end
        
        if res_p > 10*res_d
            mu = mu*2;
            D1 = D1/2;
            D2 = D2/2;
            D3 = D3/2;
            D4 = D4/2;
        elseif res_d > 10*res_p
            mu = mu/2;
            D1 = D1*2;
            D2 = D2*2;
            D3 = D3*2;
            D4 = D4*2;
        end
    end
    i = i + 1;
end
if i == MaxIter + 1
    display('Maximum iteration reached!');
end
X = U;
end