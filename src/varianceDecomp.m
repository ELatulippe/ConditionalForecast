function [shares, MSPE, tabHoriz] = varianceDecomp(model, Hmax, D, hReport, varNames, shockNames)
% VARIANCEDECOMP  Forecast-error variance decomposition (Section 2.4).
%
%   [shares, MSPE] = varianceDecomp(model, Hmax, D)
%
% Implements Eqs. (7')-(7''):
%   MSPE(h)        = sum_{j=0}^{h-1} Theta_j Theta_j'
%   MSPE_{k,l}(h)  = sum_{j=0}^{h-1} theta_{k,l,j}^2
%   share_{k,l}(h) = MSPE_{k,l}(h) / sum_l MSPE_{k,l}(h)
%
% INPUTS
%   model      struct from estimateVAR
%   Hmax       maximum horizon
%   D          (optional) impact matrix; default chol(Sigma_u,'lower')
%   hReport    (optional) horizons at which to print a table, e.g. [3 12 24 48]
%   varNames   (optional) 1xK cellstr of variable names (rows of the table)
%   shockNames (optional) 1xK cellstr of shock names (columns)
%
% OUTPUTS
%   shares    K x K x Hmax   shares(k,l,h) = share of h-step FEV of variable k
%                            due to shock l  (each row k sums to 1 over l)
%   MSPE      K x Hmax       total h-step FEV of each variable
%   tabHoriz  struct array   printed decomposition at hReport horizons

    K = model.K;
    if nargin < 3 || isempty(D), D = chol(model.Sigmau,'lower'); end
    MA    = varMA(model, Hmax, D);
    Theta = MA.Theta;                         % K x K x Hmax

    shares = zeros(K,K,Hmax);
    MSPE   = zeros(K,Hmax);
    cumContrib = zeros(K,K);                   % running sum_j theta^2
    for h = 1:Hmax
        cumContrib = cumContrib + Theta(:,:,h).^2;   % add theta_{.,.,h-1}^2
        tot = sum(cumContrib, 2);              % K x 1, sum over shocks l
        MSPE(:,h) = tot;
        % bsxfun (not ./) so this works before R2016b and under Octave.
        shares(:,:,h) = bsxfun(@rdivide, cumContrib, tot);   % normalise each row
    end

    tabHoriz = struct([]);
    if nargin >= 4 && ~isempty(hReport)
        if nargin < 5 || isempty(varNames)
            varNames = arrayfun(@(k)sprintf('y%d',k),1:K,'uni',0);
        end
        if nargin < 6 || isempty(shockNames)
            shockNames = arrayfun(@(k)sprintf('e%d',k),1:K,'uni',0);
        end
        for ii = 1:numel(hReport)
            h = hReport(ii);
            S = 100*shares(:,:,h);
            fprintf('\n=== Variance decomposition, horizon h = %d (%% of FEV) ===\n', h);
            fprintf('%-14s', 'variable\shock');
            fprintf('%12s', shockNames{:}); fprintf('\n');
            for k = 1:K
                fprintf('%-14s', varNames{k});
                fprintf('%12.1f', S(k,:)); fprintf('\n');
            end
            tabHoriz(ii).h = h; tabHoriz(ii).shares = S;
        end
    end
end
