function X = BiJSpLRU(Y, F, parameter)

K0      = parameter.K0;
numpat  = parameter.numpat;
gamma   = parameter.gamma;
tau     = parameter.tau;
mu      = parameter.mu;
verbose = parameter.verbose;
MaxIter = parameter.MaxIter;


norm_y = sqrt(mean(Y(:).^2));
Y = Y./norm_y;
F = F./norm_y;

epsilon = 1e-5;
[L, K] = size(Y);
m = size(F, 2);

Finv = (F'*F +  4*mu*eye(m))^-1;
W = Finv*F'*Y;

%Initialization of auxiliary matrices V1,V2,V3,V4
V1 = W;
V2 = W;
V3 = W;
V4 = W;

%Initialization of Lagranfe Multipliers
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


mu_changed = 0;

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
    
    % W
    W = Finv*(F'* Y + mu*((V1 + D1) + (V2 + D2) + (V3 + D3)+(V4 + D4)));
    
    %V1
    V1_temp = W-D1;    

    V1_temp = rtc(V1_temp,sqrt(size(Y,2)),sqrt(size(Y,2)));

    col_k   = 1:numpat(1);
    V1(:,col_k) = vector_soft_row(V1_temp(:,col_k), ...
        gamma/mu* (1./( sqrt(sum(V1_temp(:,col_k).^2,2)) + eps ))+ eps );
    for k = 2:K0   % K0 = length(numpat)
        col_k = (col_k(end)+1) : (col_k(end)+ numpat(k));
        V1(:,col_k) = vector_soft_row(V1_temp(:,col_k), ...
            gamma/mu* (1./( sqrt(sum(V1_temp(:,col_k).^2,2)) + eps ))+ eps );
    end

    V1 = rtc(V1,sqrt(size(Y,2)),sqrt(size(Y,2)));
    % V2
    V2_temp = W-D2;      

    col_k   = 1:numpat(1);
    V2(:,col_k) = vector_soft_row(V2_temp(:,col_k), ...
        gamma/mu* (1./( sqrt(sum(V2_temp(:,col_k).^2,2)) + eps ))+ eps );
    for k = 2:K0   % K0 = length(numpat)
        col_k = (col_k(end)+1) : (col_k(end)+ numpat(k));
        V2(:,col_k) = vector_soft_row(V2_temp(:,col_k), ...
            gamma/mu* (1./( sqrt(sum(V2_temp(:,col_k).^2,2)) + eps ))+ eps );
    end 

    
    %V3
    [u,s,v] = svd(W - D3,'econ');
    ds = diag(s);
    V3 = u*diag( max( abs(ds) - (tau/mu)*(1./(abs(ds)+ eps)),...
        zeros(size(ds)) ) )*v';
  
    %V4
    V4 = max(W - D4, 0);
    
    %update D
    D1 = D1 - W + V1;
    D2 = D2 - W + V2;
    D3 = D3 - W + V3;
    D4 = D4 - W + V4;
    
    if mod(i, 10) == 1       
        %primal residual
        res_p = norm([V1; V2; V3; V4] - [W; W; W; W], 'fro');
        %dual residual
        res_d = mu*norm([V1; V2; V3; V4] - [V10; V20; V30; V40], 'fro');
        
        if verbose
            fprintf('i = %d, res_p = %f, res_d = %f, mu = %f\n',...
                i, res_p, res_d, mu);
        end
        
        if res_p > 10*res_d
            mu = mu*2;
            D1 = D1/2;
            D2 = D2/2;
            D3 = D3/2;
            D4 = D4/2;
            mu_changed = 1;
        elseif res_d > 10*res_p
            mu = mu/2;
            D1 = D1*2;
            D2 = D2*2;
            D3 = D3*2;
            D4 = D4*2;
            mu_changed = 1;
        end
        if  mu_changed
            % update Finv
            Finv = (F'*F +  4*mu*eye(m))^-1;
            mu_changed = 0;
            % mu
        end
        
    end
 
    i = i + 1;
end

X = W.*(W>=0);

end
