function [A_hat E_hat iter] = WNNM(D, Support, C, myeps, maxIter, tol)

% April 2015
%
% This matlab code implements the inexact augmented Lagrange multiplier 
% method for WNNM-MC in 
% "Weighted nuclear norm minimization and its applications in low level
%  vision".
% S. Gu, Q. Xie, D. Meng, W. Zuo, X. Feng, L. Zhang.
%
% The program is written based on the code provided in
%
%    "The Augmented Lagrange Multiplier Method for Exact Recovery of 
%    Corrupted Low-Rank Matrices". 
%    Z. Lin, M. Chen, L. Wu. arXiv:1009.5055, 2010
%
%
% D - m x n matrix of observations/data (required input)
%
% Support - observation data indicator (binary matrix) 
%
% C - parameter for the weight setting in the weighted nuclear norm
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
%
% maxIter - maximum number of iterations
%         - DEFAULT 1000, if omitted or -1.
% 
% Initialize A,E,Y,u
% while ~converged 
%   minimize (inexactly, update A and E only once)
%     L(A,E,Y,u) = |A|_w,*  + <Y,D-A-E> + mu/2 * |D-A-E|_F^2;
%                                  s.t. E.*Support = 0;
%   Y = Y + \mu * (D - A - E);
%   \mu = \rho * \mu;
% end


% addpath PROPACK;

UpdatingSupport = ~Support;
[m n] = size(D);

if nargin < 4
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 5
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

% initialize
Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) ;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

A_hat = zeros( m, n);
E_hat = zeros( m, n);
mu = 1/norm_two; % this one can be tuned
mu_bar = mu * 1e7;
rho = 1.05;          % this one can be tuned
d_norm = norm(D, 'fro');

iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;
sv = 10;
while ~converged       
    iter = iter + 1;
    
    temp_T = D - A_hat + (1/mu)*Y;
    E_hat = UpdatingSupport.*temp_T;

    if choosvd(n, sv) == 1
        [U S V] = lansvd(D - E_hat + (1/mu)*Y, sv, 'L');
    else
        [U S V] = svd(D - E_hat + (1/mu)*Y, 'econ');
    end

    diagS = diag(S);
    [tempDiagS,svp]=ClosedWNNM(diagS,C/mu,myeps);
    A_hat = U(:,1:svp)*diag(tempDiagS)*V(:,1:svp)';  
    %%%%%%%%%%%%
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end

%%%%%% 
    total_svd = total_svd + 1;
  
    Z = D - A_hat - E_hat;
    
    Y = Y + mu*Z;
    mu = min(mu*rho, mu_bar);
        
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / d_norm;
    if stopCriterion < tol
        converged = true;
    end    
    
    % if mod( total_svd, 50) == 0
    %     disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(rank(A_hat))...
    %         ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
    %         ' stopCriterion ' num2str(stopCriterion)]);
    % end    
    
    if ~converged && iter >= maxIter
        % disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end

end


function [SigmaX,svp]=ClosedWNNM(SigmaY,C,oureps)
temp=(SigmaY-oureps).^2-4*(C-oureps*SigmaY);
ind=find (temp>0);
svp=length(ind);
SigmaX=max(SigmaY(ind)-oureps+sqrt(temp(ind)),0)/2;

end

function y = choosvd( n, d)

if n <= 100 
    if d / n <= 0.02
        y = 1;
    else
        y = 0;
    end
elseif n <= 200
    if d / n <= 0.06
        y = 1;
    else
        y = 0;
    end
elseif n <= 300
    if d / n <= 0.26
        y = 1;
    else
        y = 0;
    end
elseif n <= 400
    if d / n <= 0.28
        y = 1;
    else
        y = 0;
    end
elseif n <= 500
    if d / n <= 0.34
        y = 1;
    else
        y = 0;
    end
else
    if d / n <= 0.38
        y = 1;
    else
        y = 0
    end
end

end
