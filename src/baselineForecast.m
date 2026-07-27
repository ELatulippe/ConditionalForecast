function [Fbase, Ystate] = baselineForecast(model, Ylast, H)
% BASELINEFORECAST  Unconditional VAR forecast E[y_{t+h}|F_t], h=1..H.
%
%   [Fbase, Ystate] = baselineForecast(model, Ylast, H)
%
% Implements Eqs. (3)-(4): iterate the companion recursion with the shock
% expectations set to zero.
%
% INPUTS
%   model   struct from estimateVAR
%   Ylast   p x K matrix of the last p observations (rows newest-last),
%           i.e. Ylast = Y(end-p+1:end,:). Used to build the state Y_t.
%   H       forecast horizon.
%
% OUTPUTS
%   Fbase   H x K   baseline forecasts (transformed units), row h = horizon h
%   Ystate  Kp x 1  companion state Y_t = [y_t; y_{t-1}; ...; y_{t-p+1}]

    K = model.K; p = model.p; A = model.A; mu = model.mu; J = model.J;

    % Build companion state Y_t: newest block on top.
    Ystate = zeros(K*p,1);
    for i = 1:p
        Ystate((i-1)*K+1:i*K) = Ylast(end-i+1, :).';   % y_t, y_{t-1}, ...
    end

    Fbase = zeros(H, K);
    Yh = Ystate;
    for h = 1:H
        Yh = mu + A*Yh;              % E[Y_{t+h}|F_t] = mu + A E[Y_{t+h-1}|F_t]
        Fbase(h, :) = (J*Yh).';      % top block = y_{t+h}
    end
end
