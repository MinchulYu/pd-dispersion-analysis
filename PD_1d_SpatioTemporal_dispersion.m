clear
close all
clc

scrz = get(0,'screensize');

%==========================================================================
% 1D PD spatial and spatio-temporal dispersion analysis
%
% Notation:
%   k       : numerical wavenumber
%   k0      : exact wavenumber
%   omega0  : exact angular frequency obtained from the discrete system
%   c       : numerical wave speed
%   c0      : exact wave speed
%   e       : relative wave-speed error, e = (c-c0)/c0
%   kappa   : normalized numerical wavenumber, k*Delta x
%   kappa0  : normalized exact wavenumber, k0*Delta x
%
% This code plots:
%   1) spatial dispersion relation: kappa/pi vs kappa0/pi
%   2) spatio-temporal dispersion relation: kappa/pi vs e
%
% Required functions:
%   Peri_const_mk_1d.m
%   NB_find_discard_kappa.m
%==========================================================================

%==========================================================================
% number of nodes
%==========================================================================
Nx = 21;     % 51 or 21
m = 2;       % horizon ratio, delta/h

%==========================================================================
% nodal spacing
%==========================================================================
h = 1;       % h = Delta x

%==========================================================================
% domain size
%==========================================================================
x_length = (Nx-1)*h;
delta_x = x_length/(Nx-1);

%==========================================================================
% coordinates
%==========================================================================
coor = zeros(Nx,1);
for i = 1:Nx
    coor(i,1) = delta_x*(i-1);
end

%==========================================================================
% assign horizon radius
%==========================================================================
horizon = m*delta_x;

%==========================================================================
% nodal volume
%==========================================================================
vol = zeros(Nx,1);
for i = 1:Nx
    if (coor(i,1) == 0) || (coor(i,1) == x_length)
        vol(i) = delta_x^2/2;
    else
        vol(i) = delta_x^2;
    end
end

%==========================================================================
% connectivity for PD
%
% lm_pd(i,:) stores the family nodes of center node i within the horizon.
%==========================================================================
lm_pd = 0;
for i = 1:Nx
    num_family = 0;

    for j = 1:i-1
        if sqrt((coor(i,1)-coor(j,1))^2) < 1.0001*horizon
            num_family = num_family + 1;
            lm_pd(i,num_family) = j;
        end
    end

    for j = i+1:Nx
        if sqrt((coor(i,1)-coor(j,1))^2) < 1.0001*horizon
            num_family = num_family + 1;
            lm_pd(i,num_family) = j;
        end
    end
end

%% Spatial dispersion analysis (1D PD)

syms kappa0_sym kappa_sym real

E = 1;
A = 1;
rho = 1;

c0 = 1;                         % exact wave speed
micro_modulus = 2*E/A/horizon^2;

kappa_range = pi;
kappa_num = 51;

% kappa_sym  = k*h
% kappa0_sym = k0*h
k  = kappa_sym/h;       % numerical wavenumber
k0 = kappa0_sym/h;      % exact wavenumber

K = zeros(Nx,Nx);
M = eye(Nx,Nx);
M(1,1) = 1/2;
M(end,end) = 1/2;

% Assemble PD stiffness matrix
for i = 1:Nx
    xi = coor(i,1);

    lm_length = 0;
    for a = 1:size(lm_pd,2)
        if lm_pd(i,a) ~= 0
            lm_length = a;
        end
    end

    for j = 1:lm_length
        xj = coor(lm_pd(i,j),1);

        kloc = Peri_const_mk_1d(xi,xj,horizon, ...
                                vol(i),vol(lm_pd(i,j)), ...
                                micro_modulus,delta_x);

        K(i,[i,lm_pd(i,j)]) = K(i,[i,lm_pd(i,j)]) + kloc;
    end
end

% Time-independent wave equation:
%   (K - omega0^2 M)U = 0
% With c0 = 1, omega0 = k0*c0 = k0, and kappa0 = k0*h.
Mat = K - (k0*h)^2*M;

% Plane-wave vector with numerical wavenumber k
normalized_coord = coor - coor(ceil(end/2));

U = sym(zeros(length(M),1));
for i = 1:length(M)
    U(i) = cos(k*(normalized_coord(i)));
end

MatU = Mat*U;
TargetEq_spatial = MatU(ceil(length(M)/2));

SpatialFun = matlabFunction(TargetEq_spatial, ...
                            'Vars',[kappa0_sym,kappa_sym]);

options = optimset('Display','off','TolFun',1e-16,'TolX',1e-16);

% kappa_set stores kappa = k*Delta x, where k is the numerical wavenumber.
kappa_set = linspace(1e-4,kappa_range,kappa_num);

% Solve the spatial dispersion relation to obtain kappa0 = k0*Delta x.
kappa0 = zeros(1,length(kappa_set));
for i = 1:length(kappa_set)
    kappa_now = kappa_set(i);
    kappa0(i) = fzero(@(kappa0_now) SpatialFun(kappa0_now,kappa_now), ...
                      kappa_now,options);
end

% Spatial dispersion plot: kappa/pi vs kappa0/pi
figure('position',[50 50 0.3*scrz(3) 0.4*scrz(4)])
plot(kappa_set/pi,kappa0/pi,'-b','linewidth',2); hold on
plot([0 1],[0 1],'k','linewidth',2); hold on
axis equal
axis([0 1 0 1])
xticks([0 0.25 0.5 0.75 1])
xlabel('$$k\Delta x/ \pi$$','interpreter','latex')
ylabel('$$k_0\Delta x/\pi$$','interpreter','latex')
legend('Numerical','Exact','location','northwest')
set(gca,'fontname','times new roman','fontsize',18)
grid off

% exportgraphics(gca,'PD_1dwave_spatial.tiff','resolution',1000)

%% Spatio-temporal dispersion analysis
%
% This part combines the spatial dispersion relation k-k0 with
% the selected time integration scheme.
%
% scheme_ind : 1 = CDM / 2 = NB method
% p_ind      : 1 -> p = 0.54 / 2 -> p = 2-sqrt(2) for NB method

for aaa = 3
    if aaa == 1
        scheme_ind = 1;
    elseif aaa == 2
        scheme_ind = 2;
        p_ind = 1;
    elseif aaa == 3
        scheme_ind = 2;
        p_ind = 2;
    end

    kappa_discard = 1000*ones(1,3);
    gh_plot = zeros(3,length(kappa_set));
    e_save_all = zeros(3,length(kappa_set));

    % Convert normalized wavenumbers to dimensional wavenumbers.
    k0 = kappa0/h;       % exact wavenumber
    k  = kappa_set/h;    % numerical wavenumber

    switch scheme_ind

        case 1          % Central difference method
            syms e_sym real

            CFL_set = [0.4 0.7 1];

            for iii = 1:3
                CFL = CFL_set(iii);
                e_search = 0;
                e_save = zeros(1,length(kappa_set));

                for ii = 1:length(kappa_set)
                    dt = CFL*h/c0;
                    c = c0*(1 + e_sym);     % numerical wave speed

                    % Characteristic-polynomial coefficients for CDM
                    A1 = -1/2*(k0(ii)^2*c0^2*dt^2 - 2);
                    A2 = 1;

                    % Numerical damping ratio
                    gh = -1/(2*k0(ii)*c0*dt)*log(A2);

                    % Spatio-temporal dispersion equation
                    x2 = exp(-gh*k(ii)*c*dt)*cos(k(ii)*c*dt);
                    x1 = 1;
                    x0 = exp( gh*k(ii)*c*dt)*cos(-k(ii)*c*dt);

                    TargetEq_temporal = x2 - 2*A1*x1 + A2*x0;
                    TemporalFun = matlabFunction(TargetEq_temporal, ...
                                                  'Vars',e_sym);

                    options = optimset('Display','off','TolFun',1e-16,'TolX',1e-16);

                    if isnan(TemporalFun(e_search)) || ...
                       isinf(TemporalFun(e_search)) || ...
                       ~isreal(TemporalFun(e_search))

                        e_save(1,ii-3:end) = NaN;
                        gh_plot(iii,ii-3:end) = NaN;

                    else
                        e_save(ii) = fzero(TemporalFun,e_search,options);
                        gh_plot(iii,ii) = -1/(2*k0(ii)*c0*dt)*log(A2);
                    end

                    e_search = e_save(ii);
                end

                e_save_all(iii,:) = e_save;
            end

        case 2          % Noh-Bathe method
            syms e_sym real

            if p_ind == 1
                p = 0.54;
                CFL_set = [1 1.5 1.85];
            elseif p_ind == 2
                p = 2-sqrt(2);
                CFL_set = [1 1.5 1.70];
            end

            q1 = (1-2*p)/(2*p*(1-p));

            for iii = 1:3
                CFL = CFL_set(iii);
                e_search = 0;
                e_save = zeros(1,length(kappa_set));

                for ii = 1:length(kappa_set)
                    dt = CFL*h/c0;
                    c = c0*(1 + e_sym);     % numerical wave speed

                    % Characteristic-polynomial coefficients for NB method
                    A1 = 1 ...
                        - 1/2*k0(ii)^2*c0^2*dt^2 ...
                        + 1/4*p*(1-p)*(p^2*q1-p*q1+1/2)*k0(ii)^4*c0^4*dt^4;

                    A2 = 1 ...
                        + 1/2*p*q1*(1-p)^3*k0(ii)^4*c0^4*dt^4;

                    % Numerical damping ratio
                    gh = -1/(2*k0(ii)*c0*dt)*log(A2);

                    % Spatio-temporal dispersion equation
                    x2 = exp(-gh*k(ii)*c*dt)*cos(k(ii)*c*dt);
                    x1 = 1;
                    x0 = exp( gh*k(ii)*c*dt)*cos(-k(ii)*c*dt);

                    TargetEq_temporal = x2 - 2*A1*x1 + A2*x0;
                    TemporalFun = matlabFunction(TargetEq_temporal, ...
                                                  'Vars',e_sym);

                    options = optimset('Display','off','TolFun',1e-16,'TolX',1e-16);

                    if isnan(TemporalFun(e_search)) || ...
                       isinf(TemporalFun(e_search)) || ...
                       ~isreal(TemporalFun(e_search))

                        e_save(1,ii-3:end) = NaN;
                        gh_plot(iii,ii-3:end) = NaN;

                    else
                        e_save(ii) = fzero(TemporalFun,e_search,options);
                        gh_plot(iii,ii) = -1/(2*k0(ii)*c0*dt)*log(A2);
                    end

                    e_search = e_save(ii);
                end

                e_save_all(iii,:) = e_save;

                kappa_discard(1,iii) = ...
                    NB_find_discard_kappa(e_save,kappa_set,c0,CFL);
            end
    end

    file_name1 = ['PD_m',num2str(m),'_1dwave_'];
    file_name2 = ['dispersion_analysis\data\PD_1dwave_m',num2str(m),'_'];

    if scheme_ind == 1
        file_name1 = [file_name1,'CDM'];
        file_name2 = [file_name2,'CDM_dispersion_data'];
        ylim_max1 = Inf;

    elseif scheme_ind == 2 && p_ind == 1
        file_name1 = [file_name1,'NB054'];
        file_name2 = [file_name2,'NB054_dispersion_data'];
        ylim_max1 = max(max(gh_plot)) + 0.1*max(max(gh_plot));

    elseif scheme_ind == 2 && p_ind == 2
        file_name1 = [file_name1,'NB059'];
        file_name2 = [file_name2,'NB059_dispersion_data'];
        ylim_max1 = max(max(gh_plot)) + 0.1*max(max(gh_plot));
    end

    % save(file_name2)

    %% Spatio-temporal dispersion plot
    figure('position',[50 50 0.3*scrz(3) 0.4*scrz(4)])
    colororder({'k','k'})

    yyaxis left
    plot(kappa_set(kappa_set<kappa_discard(1,1))/pi, ...
        e_save_all(1,kappa_set<kappa_discard(1,1)), ...
        '-k','linewidth',2); hold on
    plot(kappa_set(kappa_set<kappa_discard(1,2))/pi, ...
        e_save_all(2,kappa_set<kappa_discard(1,2)), ...
        '-b','linewidth',2); hold on
    plot(kappa_set(kappa_set<kappa_discard(1,3))/pi, ...
        e_save_all(3,kappa_set<kappa_discard(1,3)), ...
        '-r','linewidth',2); hold on

    plot(kappa_set(kappa_set>kappa_discard(1,1))/pi, ...
        e_save_all(1,kappa_set>kappa_discard(1,1)), ...
        '.k','linewidth',2,'markersize',10); hold on
    plot(kappa_set(kappa_set>kappa_discard(1,2))/pi, ...
        e_save_all(2,kappa_set>kappa_discard(1,2)), ...
        '.b','linewidth',2,'markersize',10); hold on
    plot(kappa_set(kappa_set>kappa_discard(1,3))/pi, ...
        e_save_all(3,kappa_set>kappa_discard(1,3)), ...
        '.r','linewidth',2,'markersize',10); hold on

    y1 = ylabel('$$(c-c_0)/c_0$$','interpreter','latex');
    axis([0 1 -0.7 0.3])

    yyaxis right
    plot(kappa_set/pi,gh_plot(2,:),'-ob', ...
        'markerindices',(1:5:kappa_num), ...
        'markersize',5,'markerfacecolor','b','linewidth',2); hold on
    plot(kappa_set/pi,gh_plot(3,:),'-or', ...
        'markerindices',(1:5:kappa_num), ...
        'markersize',5,'markerfacecolor','r','linewidth',2); hold on
    plot(kappa_set/pi,gh_plot(1,:),'-ok', ...
        'markerindices',(1:5:kappa_num), ...
        'markersize',5,'markerfacecolor','k','linewidth',2); hold on

    y2 = ylabel('numerical damping ratio $$\xi$$','interpreter','latex');
    ylim([0 ylim_max1])

    legend(['CFL = ',num2str(CFL_set(1))], ...
           ['CFL = ',num2str(CFL_set(2))], ...
           ['CFL = ',num2str(CFL_set(3))], ...
           'location','northwest')

    x1 = xlabel('$$k\Delta x/ \pi$$','interpreter','latex');

    set(gca,'fontname','times new roman','fontsize',19);
    set(x1,'fontsize',24);
    set(y1,'fontsize',24);
    set(y2,'fontsize',20);

    % exportgraphics(gca,[file_name1,'.tiff'],'resolution',1000)
    % savefig(file_name1)

end

% green color [0 0.6 0.2]