function [Fcond, info] = bkConditional(model, Fbase, Escen, H, D, MA)
% BKCONDITIONAL  Baumeister-Kilian (2014) structural conditional forecast.
%
%   [Fcond, info] = bkConditional(model, Fbase, Escen, H, D)
%
% Implements Eq. (6'): condition on a path of STRUCTURAL shocks rather than
% on variable values,
%       Fcond = Fbase + sum_{j=0}^{h-1} Theta_j * eps^scenario_{t+h-j}.
%
% Identification of D is short-run recursive (Cholesky of Sigma_u), matching
% the note; supply your own D to use a different scheme.
%
% INPUTS
%   model   struct from estimateVAR
%   Fbase   H x K baseline forecast
%   Escen   K x H matrix of structural-shock scenario values; column h holds
%           eps_{t+h} (in std-dev units, since eps has identity covariance).
%           Rows for shocks you do not want to perturb should be zero.
%   H       horizon
%   D       (optional) K x K impact matrix; default chol(Sigma_u,'lower').
%   MA      (optional) precomputed varMA(model,H,D).
%
% OUTPUTS
%   Fcond   H x K structural conditional forecast
%   info    struct with D and the structural IRFs Theta

    K = model.K;
    if nargin < 5 || isempty(D)
        D = chol(model.Sigmau, 'lower');     % recursive identification
    end
    if nargin < 6 || isempty(MA), MA = varMA(model, H, D); end
    Mstr = MA.Mstr;                          % KH x KH structural map

    epsStack = reshape(Escen, K*H, 1);       % [eps_{t+1};...;eps_{t+H}]
    dfull    = Mstr * epsStack;              % stacked deviation
    Fcond    = Fbase + reshape(dfull, K, H).';

    info = struct('D',D,'Theta',MA.Theta);
end
