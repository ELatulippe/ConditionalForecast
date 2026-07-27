function [pstar, aic] = selectLag(Y, pmax, Xexo)
% SELECTLAG  Choose VAR lag order by AIC on a common effective sample.
%
%   [pstar, aic] = selectLag(Y, pmax)
%   [pstar, aic] = selectLag(Y, pmax, Xexo)   % same exogenous block as the VAR
%
% AIC(p) = log det(Sigma_u(p)) + 2 * n_params / Teff, with all candidates
% estimated on the same sample (rows pmax+1..T) for comparability.

    if nargin < 3, Xexo = []; end
    [T,K] = size(Y);
    Teff = T - pmax;
    qExo = size(Xexo,2);
    aic = zeros(pmax,1);
    for p = 1:pmax
        Yr = Y(pmax+1:T, :);
        X  = ones(Teff, 1);
        for i = 1:p
            X = [X, Y(pmax+1-i:T-i, :)]; %#ok<AGROW>
        end
        if qExo > 0, X = [X, Xexo(pmax+1:T,:)]; end
        B = (X.'*X)\(X.'*Yr);
        U = Yr - X*B;
        Sig = (U.'*U)/Teff;
        nparams = K*(K*p + 1 + qExo);
        aic(p) = log(det(Sig)) + 2*nparams/Teff;
    end
    [~, pstar] = min(aic);
end
