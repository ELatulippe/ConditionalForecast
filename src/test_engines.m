% test_engines.m  -- validate the core VAR/forecast engines on simulated data.
addpath(pwd);
if exist('OCTAVE_VERSION','builtin'), rand('seed',7); randn('seed',7); else, rng(7); end

K = 3; p = 2; T = 400; H = 24;

% --- Build a known, stable VAR(2) and simulate data --------------------
A1 = [ 0.5 0.10 0.00; 0.00 0.40 0.10; 0.05 0.00 0.30];
A2 = [ 0.10 0.00 0.00; 0.00 0.10 0.00; 0.00 0.00 0.10];
nu = [0.2; -0.1; 0.05];
Ptrue = chol([1.0 0.3 0.1; 0.3 1.0 0.2; 0.1 0.2 1.0],'lower');

Y = zeros(T,K); Y(1:2,:) = randn(2,K);
for t=3:T
    Y(t,:) = (nu + A1*Y(t-1,:).' + A2*Y(t-2,:).' + Ptrue*randn(K,1)).';
end

model = estimateVAR(Y, p, true);
fprintf('--- Estimation ---\n');
fprintf('K=%d p=%d Teff=%d  ||Sigmau-Sigmau''||=%.2e\n', ...
       model.K, model.p, model.Teff, norm(model.Sigmau-model.Sigmau.'));

% Companion stability
ev = max(abs(eig(model.A)));
fprintf('max |eig(A)| = %.4f (should be < 1 for stability)\n', ev);

% --- Baseline forecast: must converge to unconditional mean ------------
[Fbase, Ystate] = baselineForecast(model, Y(end-p+1:end,:), H);
% Unconditional mean of stationary VAR: mu_y = (I - sum A_i)^{-1} nu
Asum = zeros(K); for i=1:p, Asum = Asum + model.Acell{i}; end
mu_y = (eye(K) - Asum) \ model.nu;
fprintf('\n--- Baseline forecast ---\n');
fprintf('||Fbase(end,:) - uncond.mean|| = %.4e  (should ->0)\n', ...
        norm(Fbase(end,:).' - mu_y));

% --- MA: Phi_0 should equal I_K ---------------------------------------
D0 = chol(model.Sigmau,'lower'); MA = varMA(model, H, D0);
fprintf('\n--- MA coefficients ---\n');
fprintf('||Phi_0 - I|| = %.2e\n', norm(MA.Phi(:,:,1)-eye(K)));

% --- WZ conditional: must hit the imposed path EXACTLY -----------------
cvar   = 2;                       % constrain variable 2
cpath  = Fbase(:,cvar) + 1.5;     % arbitrary target: baseline + 1.5 each period
[Fcond, uhat, wzinfo] = wzConditional(model, Fbase, cvar, cpath, H, MA);
fprintf('\n--- WZ conditional ---\n');
fprintf('max|Fcond(:,cvar) - target|      = %.3e  (must ~0)\n', ...
        max(abs(Fcond(:,cvar) - cpath)));
fprintf('||R*uhat - r|| (constraint resid)= %.3e\n', wzinfo.constraintResidual);
fprintf('unconstrained vars moved? mean|dF(:,1)|=%.4f  mean|dF(:,3)|=%.4f\n', ...
        mean(abs(Fcond(:,1)-Fbase(:,1))), mean(abs(Fcond(:,3)-Fbase(:,3))));

% Minimum-norm check, in the space the solution actually minimises: any other
% feasible structural sequence must have a weakly larger norm.
N  = null(wzinfo.R);
e1 = wzinfo.eps(:);
e2 = e1 + N*randn(size(N,2),1);                 % feasible perturbation
fprintf('min-norm check: ||eps||=%.4f <= ||eps+null||=%.4f  -> %d\n', ...
        norm(e1), norm(e2), norm(e1) <= norm(e2)+1e-9);

% The structural solution must be the Sigma_u-weighted one: comparing against
% the legacy reduced-form min-norm, both are feasible but the structural
% sequence has the smaller ||eps||.
[~, ~, wzred] = wzConditional(model, Fbase, cvar, cpath, H, MA, 'reduced');
fprintf('||eps||: structural=%.4f  reduced-form=%.4f  -> %d\n', ...
        wzinfo.normEps, wzred.normEps, wzinfo.normEps <= wzred.normEps+1e-9);

% --- BK conditional: single structural shock at t+1 == structural IRF ---
D = chol(model.Sigmau,'lower');
MAd = varMA(model, H, D);
Escen = zeros(K,H); Escen(1,1) = 1;          % one unit shock 1 at t+1
[Fbk, bkinfo] = bkConditional(model, Fbase, Escen, H, D, MAd);
dev = Fbk - Fbase;                            % should equal Theta_{h-1}(:,1)
irf1 = squeeze(MAd.Theta(:,1,:)).';           % H x K : response of all vars to shock1
fprintf('\n--- BK conditional ---\n');
fprintf('max|BK dev - structural IRF|     = %.3e  (must ~0)\n', ...
        max(abs(dev(:) - irf1(:))));

% --- Bands: the constrained variable must have ZERO band width ----------
% Every draw is ehat + P z with P = I - pinv(R) R, and R P = 0, so R(ehat+Pz)
% = r exactly. The imposed path therefore has no dispersion at any horizon.
sc = struct('name','test','cvar',cvar,'path',cpath);
Bs = wzBands(Y, p, [], sc, H, struct('nDraws',150,'parameter',false, ...
                                     'verbose',false,'seed',11));
wCon = max(max(abs(Bs.scen(1).cond(:,cvar,end) - Bs.scen(1).cond(:,cvar,1))));
wOth = max(max(abs(Bs.scen(1).cond(:,1,end)    - Bs.scen(1).cond(:,1,1))));
fprintf('\n--- Bands (shock uncertainty) ---\n');
fprintf('constrained-variable band width  = %.3e  (must ~0)\n', wCon);
fprintf('unconstrained-variable width     = %.4f  (must be > 0)\n', wOth);
fprintf('median tracks the point forecast = %.3e\n', ...
        max(abs(Bs.scen(1).cond(:,1,3) - Fcond(:,1))));

% With the parameters fixed, the scenario-minus-baseline difference is
% deterministic, so its band must collapse; with the bootstrap on it must not.
dFix = max(max(abs(Bs.scen(1).diff(:,1,end) - Bs.scen(1).diff(:,1,1))));
Bp = wzBands(Y, p, [], sc, H, struct('nDraws',150,'parameter',true, ...
                                     'shock',false,'biasCorrect',false, ...
                                     'verbose',false,'seed',12));
dBoot = max(max(abs(Bp.scen(1).diff(:,1,end) - Bp.scen(1).diff(:,1,1))));
fprintf('diff band, parameters fixed      = %.3e  (must ~0)\n', dFix);
fprintf('diff band, bootstrap on          = %.4f  (must be > 0)\n', dBoot);
fprintf('constrained width under bootstrap= %.3e  (must ~0)\n', ...
        max(max(abs(Bp.scen(1).cond(:,cvar,end) - Bp.scen(1).cond(:,cvar,1)))));

% --- Joint and partial-horizon constraints ------------------------------
% wzConditional has always accepted several constrained variables; check that
% both legs are hit, and that NaN horizons are genuinely left to the model
% rather than silently pinned.
cv2 = [cvar, 3];
cp2 = [cpath, Fbase(:,3) + 0.8];
[Fj, ~, ij] = wzConditional(model, Fbase, cv2, cp2, H, MA);
fprintf('\n--- Joint / partial-horizon constraints ---\n');
fprintf('joint: leg 1 err %.2e, leg 2 err %.2e\n', ...
        max(abs(Fj(:,cv2(1))-cp2(:,1))), max(abs(Fj(:,cv2(2))-cp2(:,2))));
fprintf('joint rows in R: %d (must be %d)\n', size(ij.R,1), 2*H);

cpP = cpath;  cpP(7:end) = NaN;              % constrain the first 6 only
[Fp, ~, ip] = wzConditional(model, Fbase, cvar, cpP, H, MA);
fprintf('pulse: constrained err %.2e over h=1..6\n', ...
        max(abs(Fp(1:6,cvar)-cpP(1:6))));
fprintf('pulse: rows in R = %d (must be 6)\n', size(ip.R,1));
fprintf('pulse: free horizons move away from the path: %.4f (must be > 0)\n', ...
        max(abs(Fp(7:end,cvar) - cpath(7:end))));

scP = struct('name','pulse','cvar',cvar,'path',cpP);
Bp2 = wzBands(Y, p, [], scP, H, struct('nDraws',120,'parameter',false, ...
                                       'verbose',false,'seed',13));
wCon2 = max(abs(Bp2.scen(1).cond(1:6,cvar,end)   - Bp2.scen(1).cond(1:6,cvar,1)));
wFre2 = max(abs(Bp2.scen(1).cond(7:end,cvar,end) - Bp2.scen(1).cond(7:end,cvar,1)));
fprintf('pulse bands: constrained width %.2e (must ~0), free width %.4f (must be > 0)\n', ...
        wCon2, wFre2);

% --- Soft conditions ----------------------------------------------------
% Omega = 0 must reproduce the hard solution exactly; a positive sd must give
% the constrained variable a band of its own.
[Fh2,~,~]  = wzConditional(model, Fbase, cvar, cpath, H, MA);
[F02,~,~]  = wzConditional(model, Fbase, cvar, cpath, H, MA, [], 1e-10*ones(H,1));
[Fs2,~,is2]= wzConditional(model, Fbase, cvar, cpath, H, MA, [], 0.5*ones(H,1));
fprintf('\n--- Soft conditions ---\n');
fprintf('sd->0 reproduces the hard solution: %.2e\n', max(max(abs(Fh2-F02))));
fprintf('soft ||eps||=%.4f <= hard ||eps||=%.4f -> %d\n', ...
        is2.normEps, wzinfo.normEps, is2.normEps <= wzinfo.normEps + 1e-9);
soSc = struct('name','soft','cvar',cvar,'path',cpath,'sd',0.5*ones(H,1));
Bso  = wzBands(Y, p, [], soSc, H, struct('nDraws',150,'parameter',false, ...
                                         'verbose',false,'seed',21));
fprintf('constrained band width under a soft condition = %.3f (must be > 0)\n', ...
        max(abs(Bso.scen(1).cond(:,cvar,end) - Bso.scen(1).cond(:,cvar,1))));

% --- In-sample counterfactual -------------------------------------------
% Feeding history back as the target must leave history untouched.
tCF = size(Y,1) - H;
Cid = wzCounterfactual(model, Y, arrayfun(@(i)datenum(1990,i,1),(1:size(Y,1))'), ...
                       tCF, cvar, Y(tCF+1:tCF+H,cvar), H, MA);
Calt = wzCounterfactual(model, Y, arrayfun(@(i)datenum(1990,i,1),(1:size(Y,1))'), ...
                        tCF, cvar, Y(tCF+1:tCF+H,cvar)-1, H, MA);
fprintf('\n--- In-sample counterfactual ---\n');
fprintf('target = history leaves history unchanged: %.2e\n', max(max(abs(Cid.diff))));
fprintf('altered target is followed exactly:        %.2e\n', ...
        max(abs(Calt.cf(:,cvar) - (Y(tCF+1:tCF+H,cvar)-1))));
fprintf('unconstrained variables respond:           %.4f\n', mean(abs(Calt.diff(:,1))));

% --- Variance decomposition: rows sum to 1, shares in [0,1] -------------
[shares, MSPE] = varianceDecomp(model, H, D);
rowsums = sum(shares,2);                       % K x 1 x H
fprintf('\n--- Variance decomposition ---\n');
fprintf('max|row sum - 1| over all h      = %.3e  (must ~0)\n', ...
        max(abs(rowsums(:)-1)));
fprintf('min share = %.4f, max share = %.4f (must be in [0,1])\n', ...
        min(shares(:)), max(shares(:)));

fprintf('\nALL CORE ENGINE CHECKS COMPLETE.\n');
