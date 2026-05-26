clear
close all
clc

%% =========================================================
% User settings
% =========================================================
h  = 1;        % h = Delta x
c0 = 1;        % exact wave speed

kappa_range = pi;
kappa_num   = 101;

options = optimset('Display','off', ...
                   'TolFun',1e-14, ...
                   'TolX',1e-14);

%% =========================================================
% Normalized wavenumber definition
% =========================================================
% k       : numerical wavenumber
% k0      : exact wavenumber
% c       : numerical wave speed
% c0      : exact wave speed
% e       : relative wave-speed error, e = (c-c0)/c0
% kappa   : normalized numerical wavenumber, kappa = k*h
% kappa0  : normalized exact wavenumber, kappa0 = k0*h

kappa_set = linspace(1e-6,kappa_range,kappa_num);

%% =========================================================
% Spatial dispersion relationship: f_spatial(kappa,kappa0) = 0
% =========================================================
% Available examples:
%   FE
%   PD with m = 2
%   PD with m = 3
%   PD with m = 4
%
% Important:
%   kappa  = k*h
%   kappa0 = k0*h

f_spatial =  @(kappa,kappa0) ...
     -exp(-1i*kappa) - exp(1i*kappa) - kappa0^2 + 2;

% Replace f_spatial with one of the following if needed:
%
% FE:
% f_spatial = @(kappa,kappa0) ...
%     -exp(-1i*kappa) - exp(1i*kappa) - kappa0^2 + 2;
%
% PD with m = 2:
% f_spatial = @(kappa,kappa0) ...
%     -1/4*cos(2*kappa) - cos(kappa) - kappa0^2 + 5/4;
%
% PD with m = 3:
% f_spatial = @(kappa,kappa0) ...
%     -2/9*cos(2*kappa) - 2/27*cos(3*kappa) ...
%     -4/9*cos(kappa) - kappa0^2 + 20/27;
%
% PD with m = 4:
% f_spatial = @(kappa,kappa0) ...
%     -1/8*cos(2*kappa) - 1/12*cos(3*kappa) ...
%     -1/32*cos(4*kappa) - 1/4*cos(kappa) ...
%     - kappa0^2 + 47/96;

%% =========================================================
% Solve spatial dispersion relation for kappa0
% =========================================================
kappa0_set = zeros(size(kappa_set));

for i = 1:length(kappa_set)

    kappa_now = kappa_set(i);

    kappa0_set(i) = fzero(@(kappa0_now) ...
                          f_spatial(kappa_now,kappa0_now), ...
                          kappa_now,options);
end

%% =========================================================
% Spatio-temporal dispersion relationship: f_temporal(CFL,e,kappa,kappa0) = 0
% =========================================================
% Important:
%   e = (c-c0)/c0
%   c = c0*(1+e)

% -------------------------------
% CDM setting
% -------------------------------
CFL_set = [0.4 0.7 1.0];

f_temporal = @(CFL,e,kappa,kappa0) ...
    2*cos(CFL*kappa*(1+e)) + CFL^2*kappa0^2 - 2;

% -------------------------------
% NB method setting
% -------------------------------
% p = 0.54;
% CFL_set = [1 1.5 1.85];
% 
% f_temporal = @(CFL,e,kappa,kappa0) ...
%     exp(-25469*kappa*(e+1)/(36744*kappa0)) ...
%     .* cos(CFL*kappa*(e+1)) ...
%     .* (4 - 2*(p-1/2)*(p-1)^2*CFL^4*kappa0^4) ...
%        .^(((1/2)*kappa*e + (1/2)*kappa)/kappa0) ...
%   + (1 - (1/2)*(p-1/2)*(p-1)^2*CFL^4*kappa0^4) ...
%     .* cos(CFL*kappa*(e+1)) ...
%     .* exp(25469*kappa*(e+1)/(36744*kappa0)) ...
%     .* (4 - 2*(p-1/2)*(p-1)^2*CFL^4*kappa0^4) ...
%        .^((-(1/2)*e - 1/2)*kappa/kappa0) ...
%   - 2 ...
%   + ((1/2)*p^3 - (1/2)*p^2)*CFL^4*kappa0^4 ...
%   + CFL^2*kappa0^2;

%% =========================================================
% Solve spatio-temporal dispersion relation for e
% =========================================================
e_save = zeros(length(CFL_set),length(kappa_set));

for icfl = 1:length(CFL_set)

    CFL = CFL_set(icfl);
    e_guess = 0;

    for i = 1:length(kappa_set)

        kappa_now  = kappa_set(i);
        kappa0_now = kappa0_set(i);

        fun_e = @(e_now) f_temporal(CFL,e_now,kappa_now,kappa0_now);

        try
            e_save(icfl,i) = fzero(fun_e,e_guess,options);
            e_guess = e_save(icfl,i);   % continuation
        catch
            e_save(icfl,i) = NaN;
            e_guess = 0;
        end
    end
end

%% =========================================================
% Plot: spatio-temporal dispersion relation
% =========================================================
figure('Position',[100 100 600 450])
hold on
box on

plot(kappa_set/pi,e_save(1,:),'-k','LineWidth',2)
plot(kappa_set/pi,e_save(2,:),'-b','LineWidth',2)
plot(kappa_set/pi,e_save(3,:),'-r','LineWidth',2)

x1 = xlabel('$$k\Delta x/\pi$$','Interpreter','latex');
y1 = ylabel('$$(c-c_0)/c_0$$','Interpreter','latex');

legend(['CFL = ',num2str(CFL_set(1))], ...
       ['CFL = ',num2str(CFL_set(2))], ...
       ['CFL = ',num2str(CFL_set(3))], ...
       'Location','best')

set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1.0)
set(x1,'FontName','Times New Roman','FontSize',20)
set(y1,'FontName','Times New Roman','FontSize',20)

xlim([0 1])
grid on

%% =========================================================
% Plot: spatial dispersion relation
% =========================================================
figure('Position',[150 150 600 450])
hold on
box on

plot(kappa_set/pi,kappa0_set/pi,'-r','LineWidth',2)
plot([0 1],[0 1],'k','LineWidth',2)

x1 = xlabel('$$k\Delta x/\pi$$','Interpreter','latex');
y1 = ylabel('$$k_0\Delta x/\pi$$','Interpreter','latex');

legend('Numerical','Exact','Location','northwest')

set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1.0)
set(x1,'FontName','Times New Roman','FontSize',20)
set(y1,'FontName','Times New Roman','FontSize',20)

axis equal
axis([0 1 0 1])
grid on