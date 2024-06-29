function X_nf = NF(mask_image, mask, lembda, mu1, mu2, k_max, t_max, tol)
%
% This code implements the N/F-ADMM algorithm
% based on :
%   Kaixin Gao, Zheng-Hai Huang, Lulu Guo,
%   "Low-rank matrix recovery problem minimizing a new ratio of two norms 
%   approximating the rank function then using an ADMM-type solver with applications",
%   Journal of Computational and Applied Mathematics, Volume 438, 2024.
%
% Inputs:
%   mask_image: sampled image (the observed vector)
%   mask: sampled set (the known index set Omega)
%       note that: A(X) is the matrix that shares same entries 
%       value with that of whose indice are in Omega
%   lembda: (>0) regularization parameter
%   mu1: penalty parameter
%   mu2: penalty parameter
%   k_max: maxiter 
%   t_max: maxiter
%   tol: tolerance of convergence criterion
%
% Outputs:
%   X_nf: recoverd matrix, obtained by N/F-ADMM
%
% Author: Ken Chen
%

% Initializaion
X = mask_image;
M = X;
S = X;
% PICKS = find(mask == 1);
N = X;
T = X;

% iter
for k = 1: k_max

    % first for X:
    B = M - S / mu1;

    % temp store X
    Xtemp_outer = X;

    for t = 1: t_max
        
        % first for X:
        C = N - T / mu2;
        [U, Sigma, V] = svd((mu1 * B + mu2 * C) / (mu1 + mu2));

        % D(Sigma)
        [nx, ny] = size(Sigma);
        D_Sigma = zeros(nx, ny);
        for i = 1: min(nx, ny)
            subtract = lembda / ((mu1 + mu2) * norm(N, 'fro'));
            D_Sigma(i, i) = max(Sigma(i, i) - subtract, 0);
        end

        % store X(before update) for further checking tol
        Xtemp = X;

        % update X
        X = U * D_Sigma * V';

        % then for N:
        F = X + T / mu2;
        d = lembda * sum(svd(X));%%%%%%%%%%%%
        q = d / (mu2 * norm(F, 'fro')^3);
        s = power((27 * q + 2 + sqrt((27 * q + 2)^2 - 4)) / 2, 1 / 3);
        gama = 1 / 3 + (s + 1 / s) / 3;

        % update N
        if F ~= 0
            N = gama * F;
        else
            % random matrix in the article
            [npx, npy] = size(F);
            N = ones(size(F)) .* (power(d / mu2, 1/3) / sqrt(npx * npy));
        end

        % last for T: update T
        T = T + mu2 * (X - N);

        % stopping criteria
        TOLL = norm(X - Xtemp, 'fro') / norm(X, 'fro');
        if TOLL <= tol
            break;
        end

    end

    % then for M:
    % not sure%%%%%%%%
    % tool matrix
    tool_mat = mask_image + mu1 * X + S;

    % Mp + mu1*M = tool_mat
    % update M
    [mx, my] = size(M);
    for i = 1: mx
        for j = 1: my
            M(i, j) = tool_mat(i, j) / (mu1 + mask(i, j));
        end
    end

    % last for S: update S
    S = S + mu1 * (X - M);

    % stopping criteria
    TOLL = norm(X - Xtemp_outer, 'fro') / norm(X, 'fro');
    if TOLL <= tol
        break;
    end

end

X_nf = X;

end

