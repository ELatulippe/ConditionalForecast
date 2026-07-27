function out = varMA(model, H, D)
% VARMA  Moving-average / impulse-response coefficients of the VAR.
%
%   out = varMA(model, H)        reduced-form MA  Phi_j = J A^j J'
%   out = varMA(model, H, D)     also structural  Theta_j = Phi_j * D
%
% Implements the objects Phi_j (Eq. 5) and Theta_j = Phi_j D (Eq. 6').
% Returns coefficients for horizons j = 0,1,...,H-1.
%
% OUTPUT struct:
%   .Phi     K x K x H   reduced-form MA (Phi(:,:,1) = Phi_0 = I_K)
%   .Theta   K x K x H   structural MA (only if D supplied)
%   .Mred    KH x KH     block-lower-triangular stacked reduced-form map
%                        block (h,s) = Phi_{h-s} for s<=h, 0 otherwise
%   .Mstr    KH x KH     same, structural (only if D supplied)
%
% Mred maps a stacked innovation vector u=[u_{t+1};...;u_{t+H}] into the
% stacked path deviation dy=[dy_{t+1};...;dy_{t+H}] (see note, Eq. 5).

    K = model.K; A = model.A; J = model.J;

    Phi = zeros(K,K,H);
    Aj  = eye(size(A));                 % A^0
    for j = 0:H-1
        Phi(:,:,j+1) = J*Aj*J.';        % Phi_j = J A^j J'
        Aj = A*Aj;                      % advance to A^(j+1)
    end

    Mred = blockLowerTri(Phi, K, H);

    out = struct('Phi',Phi,'Mred',Mred);

    if nargin >= 3 && ~isempty(D)
        Theta = zeros(K,K,H);
        for j = 0:H-1
            Theta(:,:,j+1) = Phi(:,:,j+1)*D;   % Theta_j = Phi_j D
        end
        out.Theta = Theta;
        out.Mstr  = blockLowerTri(Theta, K, H);
    end
end

function M = blockLowerTri(C, K, H)
% Assemble KH x KH block-lower-triangular matrix with block(h,s)=C(:,:,h-s+1).
    M = zeros(K*H, K*H);
    for h = 1:H
        for s = 1:h
            M((h-1)*K+1:h*K, (s-1)*K+1:s*K) = C(:,:,h-s+1);
        end
    end
end
