function [g, addM] = monthGrid(y0,m0,y1,m1)
% MONTHGRID  Monthly datenum grid (1st of each month) from y0-m0 to y1-m1.
%   g    = column vector of datenums
%   addM = handle @(dn,k) to add k months to a datenum (1st-of-month)
    n0 = y0*12 + (m0-1);  n1 = y1*12 + (m1-1);
    idx = (n0:n1).';
    yy = floor(idx/12);  mm = mod(idx,12)+1;
    g = datenum(yy, mm, 1);
    addM = @(dn,k) addMonthsLocal(dn,k);
end

function out = addMonthsLocal(dn,k)
    dv = datevec(dn);  out = datenum(dv(1), dv(2)+k, 1);
end
