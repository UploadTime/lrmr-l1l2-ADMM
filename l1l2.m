function X_l1l2 = l1l2(mask_image, mask, beta, mu, rho, tol, maxiter)
%
% This code extends the l1/l2 fractional algorithm (tbd) to matrix form, 
% in order to design a new algorithm for solving N/F.
%
% based on :
%   C. Wang, M. Yan, Y. Rahimi and Y. Lou, "Accelerated Schemes for 
%   the  L1/L2  Minimization," in IEEE Transactions on Signal 
%   Processing, vol. 68, pp. 2660-2669, 2020.
%
% Inputs:
%   mask_image: sampled image (the observed vector)
%   mask: sampled set (the known index set Omega)
%       note that: mask == A
%   beta: penalty parameter in original form
%   mu: regularization parameter
%   rho: penalty parameter in admm
%   tol: tolerance of convergence criterion
%   maxiter: max iter
%
% Outputs:
%   X_nf: recoverd matrix, obtained by l1l2
%
% Author: Ken Chen
%

% Initialization
X = mask_image;
[n1, n2] = size(X);
% I = ones(size(X));
M = X;
S = X;

% omega_i = ||X^i||_* / ||X^i||_F = N/F(X^i)
X_N = sum(svd(X));
X_F = norm(X, 'fro');
alpha = X_N / X_F;

% iter
for iter = 1 : maxiter

    % %check
    % if mod(iter ,10) == 1
    %     iter
    % end

    % temp 
    Xtemp = X;

    %% X-subproblem
    w = alpha / X_F;

    % -C^k in the computation
    C = (w + beta)/(beta + rho) * X + rho/(beta + rho) * M - 1/(beta + rho) * S;

    [U, Sigma, V] = svd(C);

    % D(Sigma)
    [s1, s2] = size(Sigma);
    D_Sigma = zeros(s1, s2);
    for i = 1: min(s1, s2)
        subtract = Sigma(i, i) - 1/(beta + rho);
        D_Sigma(i, i) = max(subtract, 0);
    end

    % X^k+1
    X = U * D_Sigma * V';

    %% M-subproblem
    tool_matrix = mu * mask_image + rho * X + S;

    [m1, m2] = size(M);
    for i = 1: m1
        for j = 1: m2
            M(i, j) = tool_matrix(i, j) / (rho + mask(i, j) * mu);
        end
    end

    %% S
    S = S + rho * (X - M);

    %% tol
    TOLL = norm(X - Xtemp, 'fro') / norm(X, 'fro');
    if TOLL <= tol
        % fprintf('smaller than tol')
        break;
    end

    % update
    X_N = sum(svd(X));
    X_F = norm(X, 'fro');
    alpha = X_N / X_F;
end

% return result
X_l1l2 = X;

end

