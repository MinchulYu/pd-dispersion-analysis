clear
close all
clc

scrz = get(0,'screensize');

%==========================================================================
% 2D PD spatial dispersion analysis
%
% Notation:
%   k       : numerical wavenumber
%   k0      : exact wavenumber
%   omega0  : exact angular frequency obtained from the discrete system
%   kappa   : normalized numerical wavenumber, k*Delta x
%   kappa0  : normalized exact wavenumber, k0*Delta x
%
% This code plots the spatial dispersion relationship between
% kappa/pi and kappa0/pi.
%
% Required functions:
%   Peri_const_mk_2d.m
%   NB_find_discard_kappa.m
%==========================================================================

%==========================================================================
% number of node in x and y direction
%==========================================================================
Nx = 51;
Ny = Nx;
n_node = Nx*Ny;
n_nodof = 2;
n_dof = n_node*2;
m = 3;

%==========================================================================
% nodal spacing
%==========================================================================
h = 1;                  % h = Delta x = Delta y

%==========================================================================
% domain size
%==========================================================================
x_length = (Nx-1)*h;
y_length = x_length;
delta_x = x_length/(Nx-1);
delta_y = y_length/(Ny-1);

%==========================================================================
% coordinates
%==========================================================================
coor = zeros(n_node,2);       % nodal coordinate array
for j = 1:Ny
    for i = 1:Nx
        coor(i+(j-1)*Nx,1) = delta_x*(i-1);
        coor(i+(j-1)*Nx,2) = delta_y*(j-1);
    end
end

%==========================================================================
% assign horizon radius
%==========================================================================
horizon = m*delta_x;

%==========================================================================
% nodal volume
%
% Boundary and corner volumes are reduced to account for the finite
% computational domain.
%==========================================================================
vol = zeros(n_node,1);
for i = 1:n_node
    if (coor(i,1) == 0) || (coor(i,1) == x_length) || ...
       (coor(i,2) == 0) || (coor(i,2) == y_length)
        vol(i) = delta_x^2/2;
    else
        vol(i) = delta_x^2;
    end
end
vol(1) = delta_x^2/4;
vol(Nx) = delta_x^2/4;
vol(Nx*(Ny-1)+1) = delta_x^2/4;
vol(Nx*Ny) = delta_x^2/4;

%==========================================================================
% connectivity for PD
%
% lm(i,:) stores the family nodes of center node i within the horizon.
% The factor 1.001 is used to include nodes located numerically on the
% horizon boundary.
%==========================================================================
lm = 0;
for i = 1:n_node
    num_family = 0;
    for j = 1:n_node
        dist_ij = sqrt((coor(i,1)-coor(j,1))^2 + ...
                       (coor(i,2)-coor(j,2))^2);

        if dist_ij < 1.001*horizon && dist_ij > 0
            num_family = num_family + 1;
            lm(i,num_family) = j;
        end
    end
end

%% Spatial dispersion analysis (2D PD)

theta_set = [0 pi/4];   % propagation angles: 0 and 45 degrees

E = 1;
rho = 1;
nu = 1/4;

% Exact wave speeds
% vp_exact = sqrt(E/(rho*(1-nu^2)));                  % plane stress
vp_exact = sqrt(E*(1-nu)/(rho*(1-2*nu)*(1+nu)));      % plane strain
vs_exact = sqrt(E/(2*rho*(1+nu)));

% Micro-modulus
% micro_modulus = 24*E/(pi*horizon^3*(1-nu));          % 2D plane stress (YD Ha, Int J Fract 2010)
micro_modulus = 24*E/(pi*horizon^3*(1-nu-2*nu^2));    % 2D plane strain (G Zhang, Eng Fract Mech 2016)
% micro_modulus = 48*E/5/pi/horizon^3;

kappa_range = pi;       % maximum value of kappa = k*Delta x
kappa_num = 51;

% Lumped mass matrix and PD stiffness matrix
M = rho*[diag(vol), zeros(length(vol),length(vol)); ...
         zeros(length(vol),length(vol)), diag(vol)];
K = zeros(n_dof,n_dof);

% Assemble the PD stiffness matrix.
% Peri_const_mk_2d returns the local contribution associated with
% [u_i, u_j, v_i, v_j].
for i = 1:n_node
    xi = coor(i,1);
    yi = coor(i,2);

    lm_length = 0;
    for a = 1:size(lm,2)
        if lm(i,a) ~= 0
            lm_length = a;
        end
    end

    for j = 1:lm_length
        xj = coor(lm(i,j),1);
        yj = coor(lm(i,j),2);

        kloc = Peri_const_mk_2d(xi,yi,xj,yj, ...
                                horizon,vol(i),vol(lm(i,j)), ...
                                micro_modulus,delta_x);

        assemble_ind = [i lm(i,j) n_node+i n_node+lm(i,j)];

        K(i,assemble_ind) = K(i,assemble_ind) + kloc(1,:);
        K(n_node+i,assemble_ind) = K(n_node+i,assemble_ind) + kloc(2,:);
    end
end

M = sparse(M);
K = sparse(K);

% Matrix blocks used in the 2-by-2 eigenvalue reduction
MM = M(1:n_node,1:n_node);
KK1 = K(1:n_node,1:n_node);
KK2 = K(1:n_node,n_node+1:n_dof);
KK3 = K(n_node+1:n_dof,n_node+1:n_dof);

% kappa_set stores kappa = k*Delta x, where k is the numerical wavenumber.
kappa_set = linspace(1e-4,kappa_range,kappa_num);

% omega0 is the exact angular frequency obtained from the discrete system.
omega0p_save = zeros(length(theta_set),length(kappa_set));
omega0s_save = zeros(length(theta_set),length(kappa_set));

for iii = 1:length(theta_set)
    theta = theta_set(iii);

    for jjj = 1:length(kappa_set)

        kappa = kappa_set(jjj);     % kappa = k*Delta x
        k = kappa/h;                % numerical wavenumber

        % Shift the coordinate origin to the domain center.
        % This avoids large phase values and does not affect the dispersion result.
        normalized_coor = coor - [x_length/2, y_length/2];

        % Plane-wave vector for the selected propagation angle
        U = zeros(n_node,1);
        for i = 1:n_node
            U(i,1) = exp(1i*k*(normalized_coor(i,1)*cos(theta) + ...
                               normalized_coor(i,2)*sin(theta)));
        end

        % Rayleigh quotient terms for the reduced 2-by-2 eigenvalue problem
        d1 = real((U'*KK1*U)/(U'*MM*U));
        d2 = real((U'*KK2*U)/(U'*MM*U));
        d3 = real((U'*KK3*U)/(U'*MM*U));

        % The larger and smaller eigenvalues correspond to P- and S-waves.
        omega0p_save(iii,jjj) = sqrt((d1+d3)/2 + ...
                                sqrt(((d1-d3)/2)^2 + d2^2));

        omega0s_save(iii,jjj) = sqrt((d1+d3)/2 - ...
                                sqrt(((d1-d3)/2)^2 + d2^2));
    end
end

% Convert omega0 to exact wavenumber k0.
k0p_save = omega0p_save/vp_exact;
k0s_save = omega0s_save/vs_exact;

% Store normalized exact wavenumber kappa0 = k0*Delta x.
% Multiplication by h makes the result invariant with respect to h.
kappa0p_save = k0p_save*h;
kappa0s_save = k0s_save*h;

% Spatial dispersion plot: kappa/pi vs kappa0/pi
figure('position',[50 50 0.3*scrz(3) 0.4*scrz(4)])

plot(kappa_set/pi,kappa0s_save(length(theta_set),:)/pi,'.b', ...
    'linewidth',2,'markersize',10); hold on
plot(kappa_set/pi,kappa0p_save(length(theta_set),:)/pi,'.r', ...
    'linewidth',2,'markersize',10); hold on

h2 = plot(kappa_set/pi,kappa0s_save(1,:)/pi,'-ob', ...
    'linewidth',2,'markerindices',(1:5:kappa_num), ...
    'markersize',5,'markerfacecolor','b'); hold on

h1 = plot(kappa_set/pi,kappa0p_save(1,:)/pi,'-or', ...
    'linewidth',2,'markerindices',(1:5:kappa_num), ...
    'markersize',5,'markerfacecolor','r'); hold on

h3 = plot([0 1],[0 1],'k','linewidth',2); hold on

axis equal
axis([0 1 0 1])
xticks([0 0.25 0.5 0.75 1])

x1 = xlabel('$$k\Delta x/ \pi$$','interpreter','latex');
y1 = ylabel('$$k_0 \Delta x/\pi$$','interpreter','latex');

legend([h1 h2 h3],'P-wave','S-wave','Exact','location','northwest')

set(gca,'fontname','times new roman','fontsize',20)
set(x1,'fontname','times new roman','fontsize',24)
set(y1,'fontname','times new roman','fontsize',24)
grid off

% exportgraphics(gca,['PD_m',num2str(m),'_spatial.tiff'],'resolution',1000)
% savefig(['PD_m',num2str(m),'_spatial'])

%% Spatio-Temporal dispersion analysis
%
% This part combines the spatial dispersion relation k-k0 with
% the selected time integration scheme.
%
% Notation:
%   k      : numerical wavenumber
%   k0     : exact wavenumber
%   c      : numerical wave speed
%   c0     : exact wave speed
%   e      : relative wave-speed error, (c-c0)/c0
%   kappa  : normalized numerical wavenumber, k*Delta x
%   kappa0 : normalized exact wavenumber, k0*Delta x
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

    % Initialize the discarded-mode boundary and result arrays
    kappa_discard{1} = 1000*ones(2,3);
    kappa_discard{2} = 1000*ones(2,3);

    ghp_plot{1} = zeros(3,length(kappa_set));
    ghp_plot{2} = zeros(3,length(kappa_set));
    ghs_plot{1} = zeros(3,length(kappa_set));
    ghs_plot{2} = zeros(3,length(kappa_set));

    ep_save{1} = zeros(3,length(kappa_set));
    ep_save{2} = zeros(3,length(kappa_set));
    es_save{1} = zeros(3,length(kappa_set));
    es_save{2} = zeros(3,length(kappa_set));

    for theta_ind = 1:length(theta_set)
        for wave_ind = 1:2

            if wave_ind == 1
                % P-wave
                k0 = kappa0p_save(theta_ind,:)/h;   % exact wavenumber
                k  = kappa_set/h;                   % numerical wavenumber
                c0 = vp_exact;                      % exact P-wave speed
            else
                % S-wave
                k0 = kappa0s_save(theta_ind,:)/h;   % exact wavenumber
                k  = kappa_set/h;                   % numerical wavenumber
                c0 = vs_exact;                      % exact S-wave speed
            end

            switch scheme_ind

                case 1      % Central difference method
                    syms e_sym real

                    e = zeros(1,length(kappa_set));

                    CFL_set = [0.4 0.7 1];

                    for iii = 1:3
                        CFL = CFL_set(iii);
                        e_search = 0;

                        for ii = 1:length(kappa_set)
                            dt = CFL*h/c0;
                            c = c0 + c0*e_sym;  % numerical wave speed c = c0*(1+e)

                            % Characteristic-polynomial coefficients for CDM
                            A1 = -1/2*(k0(ii)^2*c0^2*dt^2 - 2);
                            A2 = 1;

                            % Numerical damping ratio
                            gh = -1/(2*k0(ii)*c0*dt)*log(A2);

                            % Temporal dispersion equation
                            x2 = exp(-gh*k(ii)*c*dt)*cos(k(ii)*c*dt);
                            x1 = 1;
                            x0 = exp( gh*k(ii)*c*dt)*cos(-k(ii)*c*dt);

                            TargetEq_temporal = x2 - 2*A1*x1 + A2*x0;
                            TemporalFun = matlabFunction(TargetEq_temporal);

                            options = optimset('Display','off','TolFun',1e-16,'TolX',1e-16);

                            if isnan(TemporalFun(e_search)) || ...
                               isinf(TemporalFun(e_search)) || ...
                               ~isreal(TemporalFun(e_search))

                                e(1,ii-3:end) = NaN;

                                if wave_ind == 1
                                    ghp_plot{theta_ind}(iii,ii-3:end) = NaN;
                                else
                                    ghs_plot{theta_ind}(iii,ii-3:end) = NaN;
                                end

                            else
                                e(ii) = fzero(TemporalFun,e_search,options);

                                if wave_ind == 1
                                    ghp_plot{theta_ind}(iii,ii) = -1/(2*k0(ii)*c0*dt)*log(A2);
                                else
                                    ghs_plot{theta_ind}(iii,ii) = -1/(2*k0(ii)*c0*dt)*log(A2);
                                end
                            end

                            e_search = e(ii);
                        end

                        if wave_ind == 1
                            ep_save{theta_ind}(iii,:) = e;
                        else
                            es_save{theta_ind}(iii,:) = e;
                        end
                    end

                case 2      % Noh-Bathe method
                    syms e_sym real

                    e = zeros(1,length(kappa_set));

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

                        for ii = 1:length(kappa_set)
                            dt = CFL*h/c0;
                            c = c0 + c0*e_sym;  % numerical wave speed c = c0*(1+e)

                            % Characteristic-polynomial coefficients for NB method
                            A1 = 1 ...
                                - 1/2*k0(ii)^2*c0^2*dt^2 ...
                                + 1/4*p*(1-p)*(p^2*q1-p*q1+1/2)*k0(ii)^4*c0^4*dt^4;

                            A2 = 1 ...
                                + 1/2*p*q1*(1-p)^3*k0(ii)^4*c0^4*dt^4;

                            % Numerical damping ratio
                            gh = -1/(2*k0(ii)*c0*dt)*log(A2);

                            % Temporal dispersion equation
                            x2 = exp(-gh*k(ii)*c*dt)*cos(k(ii)*c*dt);
                            x1 = 1;
                            x0 = exp( gh*k(ii)*c*dt)*cos(-k(ii)*c*dt);

                            TargetEq_temporal = simplify(x2 - 2*A1*x1 + A2*x0);
                            TemporalFun = matlabFunction(TargetEq_temporal);

                            options = optimset('Display','off','TolFun',1e-16,'TolX',1e-16);

                            if isnan(TemporalFun(e_search)) || ...
                               isinf(TemporalFun(e_search)) || ...
                               ~isreal(TemporalFun(e_search))

                                e(1,ii-3:end) = NaN;

                                if wave_ind == 1
                                    ghp_plot{theta_ind}(iii,ii-3:end) = NaN;
                                else
                                    ghs_plot{theta_ind}(iii,ii-3:end) = NaN;
                                end

                            else
                                e(ii) = fzero(TemporalFun,e_search,options);

                                if wave_ind == 1
                                    ghp_plot{theta_ind}(iii,ii) = -1/(2*k0(ii)*c0*dt)*log(A2);
                                else
                                    ghs_plot{theta_ind}(iii,ii) = -1/(2*k0(ii)*c0*dt)*log(A2);
                                end
                            end

                            e_search = e(ii);
                        end

                        if wave_ind == 1
                            ep_save{theta_ind}(iii,:) = e;
                            kappa_discard{theta_ind}(1,iii) = ...
                                NB_find_discard_kappa(e,kappa_set,c0,CFL);
                        else
                            es_save{theta_ind}(iii,:) = e;
                            kappa_discard{theta_ind}(2,iii) = ...
                                NB_find_discard_kappa(e,kappa_set,c0,CFL);
                        end
                    end
            end
        end
    end

    file_name1 = ['PD_m',num2str(m),'_pwave_'];
    file_name2 = ['PD_m',num2str(m),'_swave_'];
    file_name3 = ['dispersion_analysis\data\PD_elastic_m',num2str(m),'_'];

    if scheme_ind == 1
        file_name1 = [file_name1,'CDM'];
        file_name2 = [file_name2,'CDM'];
        file_name3 = [file_name3,'CDM_dispersion_data'];
        ylim_max1 = Inf;
        ylim_max2 = Inf;

    elseif scheme_ind == 2 && p_ind == 1
        file_name1 = [file_name1,'NB054'];
        file_name2 = [file_name2,'NB054'];
        file_name3 = [file_name3,'NB054_dispersion_data'];
        ylim_max1 = max(max([ghp_plot{1};ghp_plot{2}])) + ...
                    0.1*max(max([ghp_plot{1};ghp_plot{2}]));
        ylim_max2 = max(max([ghs_plot{1};ghs_plot{2}])) + ...
                    0.1*max(max([ghs_plot{1};ghs_plot{2}]));

    elseif scheme_ind == 2 && p_ind == 2
        file_name1 = [file_name1,'NB059'];
        file_name2 = [file_name2,'NB059'];
        file_name3 = [file_name3,'NB059_dispersion_data'];
        ylim_max1 = max(max([ghp_plot{1};ghp_plot{2}])) + ...
                    0.1*max(max([ghp_plot{1};ghp_plot{2}]));
        ylim_max2 = max(max([ghs_plot{1};ghs_plot{2}])) + ...
                    0.1*max(max([ghs_plot{1};ghs_plot{2}]));
    end

    % save(file_name3)

    %% P-wave spatio-temporal dispersion plot
    figure('position',[50 50 0.3*scrz(3) 0.4*scrz(4)])
    colororder({'k','k'})

    yyaxis right
    plot(kappa_set/pi,ghp_plot{length(theta_set)}(2,:),'-o','color',[0.5 0.5 1], ...
        'markerindices',(1:5:kappa_num),'markersize',5, ...
        'markerfacecolor',[0.5 0.5 1],'linewidth',2); hold on
    plot(kappa_set/pi,ghp_plot{length(theta_set)}(3,:),'-o','color',[1 0.5 0.5], ...
        'markerindices',(1:5:kappa_num),'markersize',5, ...
        'markerfacecolor',[1 0.5 0.5],'linewidth',2);
    plot(kappa_set/pi,ghp_plot{length(theta_set)}(1,:),'-o','color',[0.5 0.5 0.5], ...
        'markerindices',(1:5:kappa_num),'markersize',5, ...
        'markerfacecolor',[0.5 0.5 0.5],'linewidth',2);
    plot(kappa_set/pi,ghp_plot{1}(2,:),'-ob','markerindices',(1:5:kappa_num), ...
        'markersize',5,'markerfacecolor','b','linewidth',2);
    plot(kappa_set/pi,ghp_plot{1}(3,:),'-or','markerindices',(1:5:kappa_num), ...
        'markersize',5,'markerfacecolor','r','linewidth',2);
    plot(kappa_set/pi,ghp_plot{1}(1,:),'-ok','markerindices',(1:5:kappa_num), ...
        'markersize',5,'markerfacecolor','k','linewidth',2); hold off

    y2 = ylabel('numerical damping ratio $$\xi$$','interpreter','latex');
    ylim([0 ylim_max1])

    yyaxis left
    plot(kappa_set(kappa_set<kappa_discard{2}(1,1))/pi, ...
        ep_save{2}(1,kappa_set<kappa_discard{2}(1,1)), ...
        '-','color',[0.5 0.5 0.5],'linewidth',2); hold on
    plot(kappa_set(kappa_set<kappa_discard{2}(1,2))/pi, ...
        ep_save{2}(2,kappa_set<kappa_discard{2}(1,2)), ...
        '-','color',[0.5 0.5 1],'linewidth',2);
    plot(kappa_set(kappa_set<kappa_discard{2}(1,3))/pi, ...
        ep_save{2}(3,kappa_set<kappa_discard{2}(1,3)), ...
        '-','color',[1 0.5 0.5],'linewidth',2);

    plot(kappa_set(kappa_set>kappa_discard{2}(1,1))/pi, ...
        ep_save{2}(1,kappa_set>kappa_discard{2}(1,1)), ...
        '.','color',[0.5 0.5 0.5],'linewidth',2,'markersize',10);
    plot(kappa_set(kappa_set>kappa_discard{2}(1,2))/pi, ...
        ep_save{2}(2,kappa_set>kappa_discard{2}(1,2)), ...
        '.','color',[0.5 0.5 1],'linewidth',2,'markersize',10);
    plot(kappa_set(kappa_set>kappa_discard{2}(1,3))/pi, ...
        ep_save{2}(3,kappa_set>kappa_discard{2}(1,3)), ...
        '.','color',[1 0.5 0.5],'linewidth',2,'markersize',10);

    h1 = plot(kappa_set(kappa_set<kappa_discard{1}(1,1))/pi, ...
        ep_save{1}(1,kappa_set<kappa_discard{1}(1,1)), ...
        '-k','linewidth',2);
    h2 = plot(kappa_set(kappa_set<kappa_discard{1}(1,2))/pi, ...
        ep_save{1}(2,kappa_set<kappa_discard{1}(1,2)), ...
        '-b','linewidth',2);
    h3 = plot(kappa_set(kappa_set<kappa_discard{1}(1,3))/pi, ...
        ep_save{1}(3,kappa_set<kappa_discard{1}(1,3)), ...
        '-r','linewidth',2);

    plot(kappa_set(kappa_set>kappa_discard{1}(1,1))/pi, ...
        ep_save{1}(1,kappa_set>kappa_discard{1}(1,1)), ...
        '.k','linewidth',2,'markersize',10);
    plot(kappa_set(kappa_set>kappa_discard{1}(1,2))/pi, ...
        ep_save{1}(2,kappa_set>kappa_discard{1}(1,2)), ...
        '.b','linewidth',2,'markersize',10);
    plot(kappa_set(kappa_set>kappa_discard{1}(1,3))/pi, ...
        ep_save{1}(3,kappa_set>kappa_discard{1}(1,3)), ...
        '.r','linewidth',2,'markersize',10); hold off

    y1 = ylabel('$$(c^P-c_0^P)/c_0^P$$','interpreter','latex');
    x1 = xlabel('$$k\Delta x/ \pi$$','interpreter','latex');

    set(gca,'fontname','times new roman','fontsize',18);
    axis([0 1 -0.7 0.3])

    legend([h1 h2 h3], ...
        ['CFL = ',num2str(CFL_set(1))], ...
        ['CFL = ',num2str(CFL_set(2))], ...
        ['CFL = ',num2str(CFL_set(3))], ...
        'location','northwest');

    set(x1,'fontsize',24);
    set(y1,'fontsize',24);
    set(y2,'fontsize',20);

    % exportgraphics(gca,[file_name1,'.tiff'],'resolution',1000)
    % savefig(file_name1)

    %% S-wave spatio-temporal dispersion plot
    figure('position',[50 50 0.3*scrz(3) 0.4*scrz(4)])
    colororder({'k','k'})

    yyaxis right
    plot(kappa_set/pi,ghs_plot{length(theta_set)}(2,:),'-o','color',[0.5 0.5 1], ...
        'markerindices',(1:5:kappa_num),'markersize',5, ...
        'markerfacecolor',[0.5 0.5 1],'linewidth',2); hold on
    plot(kappa_set/pi,ghs_plot{length(theta_set)}(3,:),'-o','color',[1 0.5 0.5], ...
        'markerindices',(1:5:kappa_num),'markersize',5, ...
        'markerfacecolor',[1 0.5 0.5],'linewidth',2);
    plot(kappa_set/pi,ghs_plot{length(theta_set)}(1,:),'-o','color',[0.5 0.5 0.5], ...
        'markerindices',(1:5:kappa_num),'markersize',5, ...
        'markerfacecolor',[0.5 0.5 0.5],'linewidth',2);
    plot(kappa_set/pi,ghs_plot{1}(2,:),'-ob','markerindices',(1:5:kappa_num), ...
        'markersize',5,'markerfacecolor','b','linewidth',2);
    plot(kappa_set/pi,ghs_plot{1}(3,:),'-or','markerindices',(1:5:kappa_num), ...
        'markersize',5,'markerfacecolor','r','linewidth',2);
    plot(kappa_set/pi,ghs_plot{1}(1,:),'-ok','markerindices',(1:5:kappa_num), ...
        'markersize',5,'markerfacecolor','k','linewidth',2); hold off

    y2 = ylabel('numerical damping ratio $$\xi$$','interpreter','latex');
    ylim([0 ylim_max2])

    yyaxis left
    plot(kappa_set(kappa_set<kappa_discard{2}(2,1))/pi, ...
        es_save{2}(1,kappa_set<kappa_discard{2}(2,1)), ...
        '-','color',[0.5 0.5 0.5],'linewidth',2); hold on
    plot(kappa_set(kappa_set<kappa_discard{2}(2,2))/pi, ...
        es_save{2}(2,kappa_set<kappa_discard{2}(2,2)), ...
        '-','color',[0.5 0.5 1],'linewidth',2);
    plot(kappa_set(kappa_set<kappa_discard{2}(2,3))/pi, ...
        es_save{2}(3,kappa_set<kappa_discard{2}(2,3)), ...
        '-','color',[1 0.5 0.5],'linewidth',2);

    plot(kappa_set(kappa_set>kappa_discard{2}(2,1))/pi, ...
        es_save{2}(1,kappa_set>kappa_discard{2}(2,1)), ...
        '.','color',[0.5 0.5 0.5],'linewidth',2,'markersize',10);
    plot(kappa_set(kappa_set>kappa_discard{2}(2,2))/pi, ...
        es_save{2}(2,kappa_set>kappa_discard{2}(2,2)), ...
        '.','color',[0.5 0.5 1],'linewidth',2,'markersize',10);
    plot(kappa_set(kappa_set>kappa_discard{2}(2,3))/pi, ...
        es_save{2}(3,kappa_set>kappa_discard{2}(2,3)), ...
        '.','color',[1 0.5 0.5],'linewidth',2,'markersize',10);

    h1 = plot(kappa_set(kappa_set<kappa_discard{1}(2,1))/pi, ...
        es_save{1}(1,kappa_set<kappa_discard{1}(2,1)), ...
        '-k','linewidth',2);
    h2 = plot(kappa_set(kappa_set<kappa_discard{1}(2,2))/pi, ...
        es_save{1}(2,kappa_set<kappa_discard{1}(2,2)), ...
        '-b','linewidth',2);
    h3 = plot(kappa_set(kappa_set<kappa_discard{1}(2,3))/pi, ...
        es_save{1}(3,kappa_set<kappa_discard{1}(2,3)), ...
        '-r','linewidth',2);

    plot(kappa_set(kappa_set>kappa_discard{1}(2,1))/pi, ...
        es_save{1}(1,kappa_set>kappa_discard{1}(2,1)), ...
        '.k','linewidth',2,'markersize',10);
    plot(kappa_set(kappa_set>kappa_discard{1}(2,2))/pi, ...
        es_save{1}(2,kappa_set>kappa_discard{1}(2,2)), ...
        '.b','linewidth',2,'markersize',10);
    plot(kappa_set(kappa_set>kappa_discard{1}(2,3))/pi, ...
        es_save{1}(3,kappa_set>kappa_discard{1}(2,3)), ...
        '.r','linewidth',2,'markersize',10); hold off

    y1 = ylabel('$$(c^S-c_0^S)/c_0^S$$','interpreter','latex');
    x1 = xlabel('$$k\Delta x/ \pi$$','interpreter','latex');

    set(gca,'fontname','times new roman','fontsize',18);
    axis([0 1 -0.7 0.3])

    legend([h1 h2 h3], ...
        ['CFL = ',num2str(CFL_set(1))], ...
        ['CFL = ',num2str(CFL_set(2))], ...
        ['CFL = ',num2str(CFL_set(3))], ...
        'location','northwest');

    set(x1,'fontsize',24);
    set(y1,'fontsize',24);
    set(y2,'fontsize',20);

    % exportgraphics(gca,[file_name2,'.tiff'],'resolution',1000)
    % savefig(file_name2)

end