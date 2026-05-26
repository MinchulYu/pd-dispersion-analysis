function [kappa_discard] = NB_find_discard_kappa(e_save,kappa_set,c0,CFL)
%--------------------------------------------------------------------------
% Find the normalized wavenumber kappa at which the NB discarded-mode
% criterion is reached.
%
% The criterion is based on
%     Delta t / T_period > 0.3,
% which is used in the manuscript as a practical guideline for identifying
% high-frequency modes strongly attenuated by numerical dissipation.
%
% Input:
%   e_save    : relative wave-speed error, e = (c-c0)/c0
%   kappa_set : normalized numerical wavenumber, kappa = k*Delta x
%   c0        : exact wave speed
%   CFL       : CFL number, CFL = c0*Delta t/h
%   p_ind     : NB parameter index, kept for consistency with main scripts
%
% Output:
%   kappa_discard : normalized wavenumber where Delta t/T_period = 0.3
%--------------------------------------------------------------------------

syms kappa real

% Remove invalid high-wavenumber part if NaN appears in the dispersion curve
if any(isnan(e_save))
    a = find(isnan(e_save),1,'first');
    kappa_set = kappa_set(1:a-1);
    e_save = e_save(1:a-1);
end

% Numerical wave speed reconstructed from e
c = c0*(1 + e_save);

% Fit c(kappa) with a low-order polynomial
poly_order = 0;
for j = 1:10
    poly_order = poly_order + 1;
    poly_coef = polyfit(kappa_set,c,poly_order);
    c_fit = polyval(poly_coef,kappa_set);

    if sum(abs(c-c_fit)) < 1e-2
        break
    end
end

x = sym(zeros(poly_order+1,1));
for j = 1:poly_order+1
    x(j) = kappa^(poly_order+1-j);
end

% Practical discarded-mode criterion used in the manuscript
discard_ratio = 0.3;   % Delta t / T_period

% From Delta t/T_period = omega*Delta t/(2*pi)
% with omega = k*c = (kappa/h)*c and CFL = c0*Delta t/h:
%     Delta t/T_period = kappa*c*CFL/(2*pi*c0)
% Thus, kappa*c(kappa) - 2*pi*c0*discard_ratio/CFL = 0.
Eq_sym = poly_coef*x*kappa - 2*pi*c0*discard_ratio/CFL;

TargetEq = matlabFunction(Eq_sym);

options = optimset('Display','off','TolFun',1e-16,'TolX',1e-16);

if isempty(kappa_set)
    kappa_discard = 1000;
elseif TargetEq(kappa_set(1))*TargetEq(kappa_set(end)) < 0
    kappa_discard = fzero(TargetEq,[kappa_set(1),kappa_set(end)],options);
else
    kappa_discard = 1000;
end